/* imap_codec -- see imap_codec.h and doc/SPEC.md S1.2. */
#include "imap_codec.h"

#include <stdint.h>

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ buffers */

static int buf_reserve(imap_buf *b, size_t extra) {
    if (b->n + extra <= b->cap) return 0;
    size_t cap = b->cap ? b->cap : 64;
    while (cap < b->n + extra) cap *= 2;
    uint8_t *p = (uint8_t *)realloc(b->p, cap);
    if (!p) return -1;
    b->p = p;
    b->cap = cap;
    return 0;
}

void imap_buf_free(imap_buf *b) {
    free(b->p);
    b->p = NULL;
    b->n = b->cap = 0;
}

/* ------------------------------------------------------------------ varints */

size_t imap_put_uvarint(uint8_t *dst, uint64_t v) {
    size_t i = 0;
    while (v >= 0x80) {
        dst[i++] = (uint8_t)(v | 0x80);
        v >>= 7;
    }
    dst[i++] = (uint8_t)v;
    return i;
}

int imap_get_uvarint(const uint8_t **pp, const uint8_t *end, uint64_t *out) {
    const uint8_t *p = *pp;
    uint64_t v = 0;
    int shift = 0;
    while (p < end) {
        uint8_t byte = *p++;
        v |= (uint64_t)(byte & 0x7f) << shift;
        if (!(byte & 0x80)) { *pp = p; *out = v; return 0; }
        shift += 7;
        if (shift > 63) break;          /* over-long: malformed */
    }
    *pp = p;
    return -1;
}

static int buf_put_uvarint(imap_buf *b, uint64_t v) {
    if (buf_reserve(b, 10) != 0) return -1;
    b->n += imap_put_uvarint(b->p + b->n, v);
    return 0;
}

/* ------------------------------------------------------------------- encode */

int imap_encode_chunk(const imap_run *runs, int64_t n, imap_buf streams[IMAP_N_STREAMS]) {
    for (int i = 0; i < IMAP_N_STREAMS; i++) {
        streams[i].p = NULL;
        streams[i].n = streams[i].cap = 0;
    }
    if (n <= 0) return 0;

    /* strand is a bitmap: size it up front and write bits directly */
    size_t nbytes = (size_t)((n + 7) / 8);
    if (buf_reserve(&streams[IMAP_STREAM_STRAND], nbytes) != 0) goto fail;
    memset(streams[IMAP_STREAM_STRAND].p, 0, nbytes);
    streams[IMAP_STREAM_STRAND].n = nbytes;

    int64_t prev_a_exit = runs[0].a;
    int64_t prev_b_exit = imap_b_enter(runs[0].b, runs[0].len, runs[0].strand);

    for (int64_t i = 0; i < n; i++) {
        const imap_run *r = &runs[i];

        int64_t ga = r->a - prev_a_exit;
        if (ga < 0) goto fail;            /* axis a must not go backwards in a chunk */
        if (buf_put_uvarint(&streams[IMAP_STREAM_A], (uint64_t)ga) != 0) goto fail;

        int64_t enter = imap_b_enter(r->b, r->len, r->strand);
        int64_t gb = enter - prev_b_exit;
        if (buf_put_uvarint(&streams[IMAP_STREAM_B], imap_zigzag(gb)) != 0) goto fail;

        if (r->len < 1) goto fail;
        if (buf_put_uvarint(&streams[IMAP_STREAM_LEN], (uint64_t)(r->len - 1)) != 0) goto fail;

        if (r->strand) streams[IMAP_STREAM_STRAND].p[i >> 3] |= (uint8_t)(1u << (i & 7));

        prev_a_exit = r->a + r->len;
        prev_b_exit = imap_b_exit(r->b, r->len, r->strand);
    }
    return 0;

fail:
    for (int i = 0; i < IMAP_N_STREAMS; i++) imap_buf_free(&streams[i]);
    return -1;
}

/* ------------------------------------------------------------------- decode */

int imap_decode_chunk(const imap_buf streams[IMAP_N_STREAMS], int64_t n,
                      int64_t a0, int64_t b_enter0, imap_run *out) {
    if (n <= 0) return 0;
    if ((size_t)((n + 7) / 8) > streams[IMAP_STREAM_STRAND].n) return -1;
    /* The bases come from the chunk directory, i.e. off disk. A negative one makes the
     * overflow guards below overflow themselves, so reject it here. Everything after
     * this point stays non-negative by induction. */
    if (a0 < 0 || b_enter0 < 0) return -1;

    const uint8_t *pa  = streams[IMAP_STREAM_A].p,   *ea  = pa  + streams[IMAP_STREAM_A].n;
    const uint8_t *pb  = streams[IMAP_STREAM_B].p,   *eb  = pb  + streams[IMAP_STREAM_B].n;
    const uint8_t *pl  = streams[IMAP_STREAM_LEN].p, *el  = pl  + streams[IMAP_STREAM_LEN].n;
    const uint8_t *str = streams[IMAP_STREAM_STRAND].p;

    int64_t prev_a_exit = a0;
    int64_t prev_b_exit = b_enter0;

    for (int64_t i = 0; i < n; i++) {
        uint64_t ua, ub, ul;
        if (imap_get_uvarint(&pa, ea, &ua) != 0) return -1;
        if (imap_get_uvarint(&pb, eb, &ub) != 0) return -1;
        if (imap_get_uvarint(&pl, el, &ul) != 0) return -1;

        /* Everything above came off disk and is attacker-controlled. Each step below
         * is checked before it is taken, because signed overflow is undefined and a
         * caller handed a negative len would use it as a size. */
        uint8_t s = (uint8_t)((str[i >> 3] >> (i & 7)) & 1u);
        if (ul > (uint64_t)INT64_MAX - 1) return -1;
        int64_t len = (int64_t)ul + 1;                       /* len >= 1 by construction */

        if (ua > (uint64_t)(INT64_MAX - prev_a_exit)) return -1;
        int64_t a = prev_a_exit + (int64_t)ua;

        int64_t db = imap_unzigzag(ub);
        if (db > 0 ? (prev_b_exit > INT64_MAX - db) : (prev_b_exit < INT64_MIN - db)) return -1;
        int64_t enter = prev_b_exit + db;

        int64_t b = enter;
        if (s) { b = enter - len; }                          /* enter >= len checked below */
        if (a < 0 || b < 0) return -1;
        if (len > INT64_MAX - a || len > INT64_MAX - b) return -1;

        out[i].a = a;
        out[i].b = b;
        out[i].len = len;
        out[i].strand = s;

        prev_a_exit = a + len;
        prev_b_exit = imap_b_exit(b, len, s);
    }
    return 0;
}
