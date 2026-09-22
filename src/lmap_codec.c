/* lmap_codec -- see lmap_codec.h and doc/SPEC.md S1.2. */
#include "lmap_codec.h"

#include <stdint.h>

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ buffers */

static int buf_reserve(lmap_buf *b, size_t extra) {
    if (b->n + extra <= b->cap) return 0;
    size_t cap = b->cap ? b->cap : 64;
    while (cap < b->n + extra) cap *= 2;
    uint8_t *p = (uint8_t *)realloc(b->p, cap);
    if (!p) return -1;
    b->p = p;
    b->cap = cap;
    return 0;
}

void lmap_buf_free(lmap_buf *b) {
    free(b->p);
    b->p = NULL;
    b->n = b->cap = 0;
}

/* ------------------------------------------------------------------ varints */

size_t lmap_put_uvarint(uint8_t *dst, uint64_t v) {
    size_t i = 0;
    while (v >= 0x80) {
        dst[i++] = (uint8_t)(v | 0x80);
        v >>= 7;
    }
    dst[i++] = (uint8_t)v;
    return i;
}

int lmap_get_uvarint(const uint8_t **pp, const uint8_t *end, uint64_t *out) {
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

static int buf_put_uvarint(lmap_buf *b, uint64_t v) {
    if (buf_reserve(b, 10) != 0) return -1;
    b->n += lmap_put_uvarint(b->p + b->n, v);
    return 0;
}

/* ------------------------------------------------------------------- encode */

/* Axis-b reference points for the delta chain (SPEC 1.2).  Under order a the chain
 * follows the traversal, so a reverse run is entered at b+len and left at b; under
 * order b runs are sorted by b-start and the plain start/end is already monotone. */
static inline int64_t b_ref_enter(int order, const lmap_run *r) {
    return order == LMAP_ORDER_A ? lmap_b_enter(r->b, r->len, r->strand) : r->b;
}
static inline int64_t b_ref_exit(int order, const lmap_run *r) {
    return order == LMAP_ORDER_A ? lmap_b_exit(r->b, r->len, r->strand) : r->b + r->len;
}

int lmap_encode_chunk(const lmap_run *runs, int64_t n, int order,
                      lmap_buf streams[LMAP_N_STREAMS]) {
    for (int i = 0; i < LMAP_N_STREAMS; i++) {
        streams[i].p = NULL;
        streams[i].n = streams[i].cap = 0;
    }
    if (order != LMAP_ORDER_A && order != LMAP_ORDER_B) return -1;
    if (n <= 0) return 0;

    /* strand is a bitmap: size it up front and write bits directly */
    size_t nbytes = (size_t)((n + 7) / 8);
    if (buf_reserve(&streams[LMAP_STREAM_STRAND], nbytes) != 0) goto fail;
    memset(streams[LMAP_STREAM_STRAND].p, 0, nbytes);
    streams[LMAP_STREAM_STRAND].n = nbytes;

    int64_t prev_a_exit = runs[0].a;
    int64_t prev_b_exit = b_ref_enter(order, &runs[0]);

    for (int64_t i = 0; i < n; i++) {
        const lmap_run *r = &runs[i];

        int64_t ga = r->a - prev_a_exit;
        if (order == LMAP_ORDER_A) {
            if (ga < 0) goto fail;        /* under order a, axis a never goes backwards */
            if (buf_put_uvarint(&streams[LMAP_STREAM_A], (uint64_t)ga) != 0) goto fail;
        } else {
            if (buf_put_uvarint(&streams[LMAP_STREAM_A], lmap_zigzag(ga)) != 0) goto fail;
        }

        int64_t gb = b_ref_enter(order, r) - prev_b_exit;
        if (buf_put_uvarint(&streams[LMAP_STREAM_B], lmap_zigzag(gb)) != 0) goto fail;

        if (r->len < 1) goto fail;
        if (buf_put_uvarint(&streams[LMAP_STREAM_LEN], (uint64_t)(r->len - 1)) != 0) goto fail;

        if (r->strand) streams[LMAP_STREAM_STRAND].p[i >> 3] |= (uint8_t)(1u << (i & 7));

        prev_a_exit = r->a + r->len;
        prev_b_exit = b_ref_exit(order, r);
    }
    return 0;

fail:
    for (int i = 0; i < LMAP_N_STREAMS; i++) lmap_buf_free(&streams[i]);
    return -1;
}

/* ------------------------------------------------------------------- decode */

/* checked signed add: 0 and *out on success, -1 on overflow */
static inline int add_ok(int64_t x, int64_t d, int64_t *out) {
    if (d > 0 ? x > INT64_MAX - d : x < INT64_MIN - d) return -1;
    *out = x + d;
    return 0;
}

int lmap_decode_chunk(const lmap_buf streams[LMAP_N_STREAMS], int64_t n, int order,
                      int64_t base_a, int64_t base_b, lmap_run *out) {
    if (order != LMAP_ORDER_A && order != LMAP_ORDER_B) return -1;
    if (n <= 0) return 0;
    if ((size_t)((n + 7) / 8) > streams[LMAP_STREAM_STRAND].n) return -1;
    /* The bases come from the chunk directory, i.e. off disk.  A negative one would
     * make every later bound meaningless, so reject it here; everything after this
     * point is checked step by step and stays non-negative. */
    if (base_a < 0 || base_b < 0) return -1;

    const uint8_t *pa  = streams[LMAP_STREAM_A].p,   *ea  = pa  + streams[LMAP_STREAM_A].n;
    const uint8_t *pb  = streams[LMAP_STREAM_B].p,   *eb  = pb  + streams[LMAP_STREAM_B].n;
    const uint8_t *pl  = streams[LMAP_STREAM_LEN].p, *el  = pl  + streams[LMAP_STREAM_LEN].n;
    const uint8_t *str = streams[LMAP_STREAM_STRAND].p;

    int64_t prev_a_exit = base_a;
    int64_t prev_b_exit = base_b;

    for (int64_t i = 0; i < n; i++) {
        uint64_t ua, ub, ul;
        if (lmap_get_uvarint(&pa, ea, &ua) != 0) return -1;
        if (lmap_get_uvarint(&pb, eb, &ub) != 0) return -1;
        if (lmap_get_uvarint(&pl, el, &ul) != 0) return -1;

        /* Everything above is attacker-controlled.  Each step is checked before it is
         * taken: signed overflow is undefined, and a caller handed a negative len
         * would use it as a size. */
        uint8_t s = (uint8_t)((str[i >> 3] >> (i & 7)) & 1u);
        if (ul > (uint64_t)INT64_MAX - 1) return -1;
        int64_t len = (int64_t)ul + 1;                        /* len >= 1 */

        int64_t a;
        if (order == LMAP_ORDER_A) {
            if (ua > (uint64_t)(INT64_MAX - prev_a_exit)) return -1;
            a = prev_a_exit + (int64_t)ua;
        } else if (add_ok(prev_a_exit, lmap_unzigzag(ua), &a) != 0) return -1;

        int64_t enter;
        if (add_ok(prev_b_exit, lmap_unzigzag(ub), &enter) != 0) return -1;
        int64_t b = (order == LMAP_ORDER_A && s) ? enter - len : enter;   /* enter>=0, len>=1 */

        if (a < 0 || b < 0) return -1;
        if (len > INT64_MAX - a || len > INT64_MAX - b) return -1;

        out[i].a = a;
        out[i].b = b;
        out[i].len = len;
        out[i].strand = s;

        prev_a_exit = a + len;
        prev_b_exit = order == LMAP_ORDER_A ? lmap_b_exit(b, len, s) : b + len;
    }
    return 0;
}
