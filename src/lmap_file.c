/* lmap_file -- see lmap_file.h and doc/SPEC.md S2. */
/* POSIX.1-2008 for pread/strdup/fstat, whatever -std= the embedding build uses. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "lmap_file.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

/* ------------------------------------------- little-endian load/store helpers
 * Everything on disk is serialized through these, so the format does not depend
 * on host endianness or struct padding. */

static void st32(uint8_t *p, uint32_t v) {
    p[0]=(uint8_t)v; p[1]=(uint8_t)(v>>8); p[2]=(uint8_t)(v>>16); p[3]=(uint8_t)(v>>24);
}
static void st64(uint8_t *p, uint64_t v) {
    st32(p, (uint32_t)v); st32(p+4, (uint32_t)(v>>32));
}
static uint32_t ld32(const uint8_t *p) {
    return (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
}
static uint64_t ld64(const uint8_t *p) {
    return (uint64_t)ld32(p) | ((uint64_t)ld32(p+4) << 32);
}

/* zlib's crc32 takes a uInt (32-bit), so a single call silently checksums only
 * len mod 2^32 bytes.  Feed it in bounded steps instead. */
static uint32_t crc_all(const uint8_t *p, uint64_t n) {
    uLong c = crc32(0L, Z_NULL, 0);
    while (n) {
        uInt k = n > 0x40000000u ? 0x40000000u : (uInt)n;
        c = crc32(c, p ? p : (const Bytef *)"", k);
        if (p) p += k;
        n -= k;
    }
    return (uint32_t)c;
}

/* growable byte vector */
typedef struct { uint8_t *p; size_t n, cap; } vec;

static int vec_need(vec *v, size_t extra) {
    if (v->n + extra <= v->cap) return 0;
    size_t cap = v->cap ? v->cap : 256;
    while (cap < v->n + extra) cap *= 2;
    uint8_t *q = realloc(v->p, cap);
    if (!q) return -1;
    v->p = q; v->cap = cap; return 0;
}
static int vec_put(vec *v, const void *src, size_t n) {
    if (vec_need(v, n) != 0) return -1;
    memcpy(v->p + v->n, src, n); v->n += n; return 0;
}
static int vec_put32(vec *v, uint32_t x) { uint8_t b[4]; st32(b,x); return vec_put(v,b,4); }
static int vec_put64(vec *v, uint64_t x) { uint8_t b[8]; st64(b,x); return vec_put(v,b,8); }
static void vec_free(vec *v) { free(v->p); v->p=NULL; v->n=v->cap=0; }

/* ------------------------------------------------ per-stream raw deflate
 * Raw deflate (no zlib wrapper): each stream would otherwise pay a 2-byte header and a
 * 4-byte Adler checksum, and a file has four streams per chunk and ~10^4 chunks.  The
 * section CRC already covers integrity. */
static int raw_deflate(const uint8_t *in, size_t n, uint8_t **out, size_t *outn) {
    z_stream z; memset(&z, 0, sizeof z);
    if (deflateInit2(&z, 9, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY) != Z_OK) return -1;
    uLong cap = deflateBound(&z, (uLong)n);
    uint8_t *buf = malloc(cap ? cap : 1);
    if (!buf) { deflateEnd(&z); return -1; }
    z.next_in = (Bytef *)in; z.avail_in = (uInt)n;
    z.next_out = buf; z.avail_out = (uInt)cap;
    int rc = deflate(&z, Z_FINISH);
    size_t got = cap - z.avail_out;
    deflateEnd(&z);
    if (rc != Z_STREAM_END) { free(buf); return -1; }
    *out = buf; *outn = got;
    return 0;
}

/* Inflate exactly `want` bytes from exactly `n` bytes of input, or fail. */
static int raw_inflate(const uint8_t *in, size_t n, uint8_t *out, size_t want) {
    z_stream z; memset(&z, 0, sizeof z);
    if (inflateInit2(&z, -15) != Z_OK) return -1;
    z.next_in = (Bytef *)in; z.avail_in = (uInt)n;
    z.next_out = out; z.avail_out = (uInt)want;
    int rc = inflate(&z, Z_FINISH);
    int ok = rc == Z_STREAM_END && z.avail_out == 0 && z.avail_in == 0;
    inflateEnd(&z);
    return ok ? 0 : -1;
}

/* ================================================================== WRITER */

typedef struct {
    char *name; uint64_t length, amin, amax; uint32_t first_chunk, n_chunks; int seen;
    int64_t last_a_exit; int has_runs;       /* axis-a order is per member, across chunks */
} wmember;

struct lmap_writer {
    FILE    *fp;
    char    *path, *profile, *schema;
    uint32_t count;  uint64_t bspan;  int codec;  int order;  int a_overlap;

    wmember *mem[2];  uint32_t n_mem[2], cap_mem[2];

    lmap_run *pend;  uint32_t n_pend, cap_pend;      /* runs in the open chunk */
    uint32_t  pend_amem, pend_bmem;
    uint64_t  pend_bmin, pend_bmax;

    lmap_chunk *chunks; uint32_t n_chunks, cap_chunks;

    /* application sections (SPEC 2.5): stored and checksummed, never interpreted */
    struct xsec { char *id; uint8_t *p; size_t n; } *xs;  uint32_t n_xs, cap_xs;

    /* order b: one whole axis-a member buffered, with each run's axis-b member */
    struct brun { lmap_run r; uint32_t bm; } *bbuf;  uint32_t n_bbuf, cap_bbuf, bbuf_amem;
    vec runs;                                        /* the runs section body */
    int failed;
};

lmap_writer *lmap_writer_open(const char *path, const char *profile, const char *schema) {
    lmap_writer *w = calloc(1, sizeof *w);
    if (!w) return NULL;
    w->fp = fopen(path, "wb");
    if (!w->fp) { free(w); return NULL; }
    w->path = strdup(path);
    w->profile = strdup(profile ? profile : "");
    w->schema = strdup(schema ? schema : "");
    w->count = 8192; w->bspan = 1000000; w->codec = LMAP_CODEC_DEFLATE; w->order = LMAP_ORDER_A;
    if (!w->path || !w->profile || !w->schema) { lmap_writer_abort(w); return NULL; }
    return w;
}

static int writer_has_runs(const lmap_writer *w);

int lmap_writer_set_params(lmap_writer *w, uint32_t count, uint64_t bspan, int codec) {
    if (!w || writer_has_runs(w)) return -1;            /* only before any run */
    if (count == 0 || count > LMAP_MAX_CHUNK_RUNS) return -1;
    if (codec != LMAP_CODEC_NONE && codec != LMAP_CODEC_DEFLATE) return -1;
    w->count = count; w->bspan = bspan; w->codec = codec; return 0;
}

static int writer_has_runs(const lmap_writer *w) { return w->n_chunks || w->n_pend || w->n_bbuf; }

int lmap_writer_set_order(lmap_writer *w, int order) {
    if (!w || writer_has_runs(w)) return -1;            /* only before any run */
    if (order != LMAP_ORDER_A && order != LMAP_ORDER_B) return -1;
    if (w->a_overlap && order != LMAP_ORDER_B) return -1;
    w->order = order; return 0;
}

int lmap_writer_set_a_overlap(lmap_writer *w) {
    /* Order a stores a as an unsigned delta from the previous run's end, which cannot
     * go backwards; order b stores it zigzag against the previous run's start. */
    if (!w || writer_has_runs(w) || w->order != LMAP_ORDER_B) return -1;
    w->a_overlap = 1; return 0;
}

int lmap_writer_add_section(lmap_writer *w, const char *id, const void *data, size_t n) {
    if (!w || w->failed || !id || strncmp(id, "x.", 2) != 0 || strlen(id) > 255) return -1;
    if (n && !data) return -1;
    for (uint32_t i = 0; i < w->n_xs; i++) if (!strcmp(w->xs[i].id, id)) return -1;  /* unique */
    if (w->n_xs == w->cap_xs) {
        uint32_t cap = w->cap_xs ? w->cap_xs * 2 : 4;
        struct xsec *nx = realloc(w->xs, cap * sizeof *nx);
        if (!nx) return -1;
        w->xs = nx; w->cap_xs = cap;
    }
    struct xsec *x = &w->xs[w->n_xs];
    x->id = strdup(id);
    x->p = malloc(n ? n : 1);
    if (!x->id || !x->p) { free(x->id); free(x->p); return -1; }
    if (n) memcpy(x->p, data, n);
    x->n = n;
    w->n_xs++;
    return 0;
}

int32_t lmap_writer_add_member(lmap_writer *w, int axis, const char *name, uint64_t length) {
    if (!w || (axis != 0 && axis != 1) || !name) return -1;
    if (w->n_mem[axis] == w->cap_mem[axis]) {
        uint32_t cap = w->cap_mem[axis] ? w->cap_mem[axis]*2 : 16;
        wmember *m = realloc(w->mem[axis], cap * sizeof *m);
        if (!m) return -1;
        w->mem[axis] = m; w->cap_mem[axis] = cap;
    }
    wmember *m = &w->mem[axis][w->n_mem[axis]];
    memset(m, 0, sizeof *m);
    m->name = strdup(name);
    if (!m->name) return -1;
    m->length = length;
    return (int32_t)w->n_mem[axis]++;
}

/* Encode and append the pending chunk. */
/* Encode n runs, already in stored order, as one chunk of members (am, bm). */
static int emit_chunk(lmap_writer *w, const lmap_run *runs, uint32_t n, uint32_t am, uint32_t bm) {
    if (n == 0) return 0;
    lmap_buf st[LMAP_N_STREAMS];
    if (lmap_encode_chunk(runs, n, w->order, st) != 0) return -1;

    /* Chunk blob (SPEC 2.2): u32 raw_len[4], u32 stored_len[4], then the four stored
     * streams.  A stream is deflated only where that makes it smaller; otherwise it is
     * stored verbatim, and stored_len == raw_len says so -- there is no separate flag. */
    uint8_t *packed[LMAP_N_STREAMS] = {0};
    size_t   stored[LMAP_N_STREAMS];
    uint64_t rawsum = 0;
    int bad = 0;
    for (int i = 0; i < LMAP_N_STREAMS; i++) {
        rawsum += st[i].n;
        stored[i] = st[i].n;
        if (w->codec == LMAP_CODEC_DEFLATE && st[i].n > 0) {
            uint8_t *z = NULL; size_t zn = 0;
            if (raw_deflate(st[i].p, st[i].n, &z, &zn) != 0) { bad = 1; break; }
            if (zn < st[i].n) { packed[i] = z; stored[i] = zn; } else free(z);
        }
    }
    vec blob = {0};
    if (!bad && 16 + rawsum > LMAP_MAX_CHUNK_RAW) bad = 1;
    for (int i = 0; !bad && i < LMAP_N_STREAMS; i++) bad |= vec_put32(&blob, (uint32_t)st[i].n);
    for (int i = 0; !bad && i < LMAP_N_STREAMS; i++) bad |= vec_put32(&blob, (uint32_t)stored[i]);
    for (int i = 0; !bad && i < LMAP_N_STREAMS; i++)
        bad |= vec_put(&blob, packed[i] ? packed[i] : st[i].p, stored[i]);
    for (int i = 0; i < LMAP_N_STREAMS; i++) { free(packed[i]); lmap_buf_free(&st[i]); }
    if (bad) { vec_free(&blob); return -1; }
    const uint8_t *body = blob.p; size_t bodylen = blob.n;
    uint32_t rawlen32 = (uint32_t)(16 + rawsum);

    if (w->n_chunks == w->cap_chunks) {
        uint32_t cap = w->cap_chunks ? w->cap_chunks*2 : 256;
        lmap_chunk *c = realloc(w->chunks, cap * sizeof *c);
        if (!c) { vec_free(&blob); return -1; }
        w->chunks = c; w->cap_chunks = cap;
    }
    /* extents from min/max: the stored order is only sorted on its own axis */
    int64_t amin = runs[0].a, amax = runs[0].a + runs[0].len;
    int64_t bmin = runs[0].b, bmax = runs[0].b + runs[0].len;
    for (uint32_t k = 1; k < n; k++) {
        if (runs[k].a < amin) amin = runs[k].a;
        if (runs[k].a + runs[k].len > amax) amax = runs[k].a + runs[k].len;
        if (runs[k].b < bmin) bmin = runs[k].b;
        if (runs[k].b + runs[k].len > bmax) bmax = runs[k].b + runs[k].len;
    }
    uint64_t a_span64 = (uint64_t)(amax - amin);
    uint64_t b_span64 = (uint64_t)(bmax - bmin);
    if (a_span64 > 0xFFFFFFFFull || b_span64 > 0xFFFFFFFFull) { vec_free(&blob); return -1; }

    lmap_chunk *c = &w->chunks[w->n_chunks];
    memset(c, 0, sizeof *c);
    c->a_member = am; c->b_member = bm;
    c->a_min = (uint64_t)amin;
    c->a_span = (uint32_t)a_span64;
    c->b_min = (uint64_t)bmin;
    c->b_span = (uint32_t)b_span64;
    /* the first run in stored order is the min on the ordering axis, so only the other
     * axis's base needs storing */
    c->base = w->order == LMAP_ORDER_A
            ? (uint64_t)lmap_b_enter(runs[0].b, runs[0].len, runs[0].strand)
            : (uint64_t)runs[0].a;
    c->off = w->runs.n;
    c->clen = (uint32_t)bodylen;
    c->rawlen = rawlen32;
    c->n_runs = n;
    c->codec = (uint32_t)w->codec;
    c->crc = crc_all(body, bodylen);

    int rc = vec_put(&w->runs, body, bodylen);
    vec_free(&blob);
    if (rc != 0) return -1;

    w->n_chunks++;
    return 0;
}


/* Order a: the pending list, cut in axis-a order as runs arrived, is one chunk. */
static int flush_chunk(lmap_writer *w) {
    if (w->n_pend == 0) return 0;
    if (emit_chunk(w, w->pend, w->n_pend, w->pend_amem, w->pend_bmem) != 0) return -1;
    w->n_pend = 0;
    return 0;
}

static int cmp_brun(const void *x, const void *y) {
    const struct brun *p = (const struct brun *)x, *q = (const struct brun *)y;
    if (p->bm != q->bm) return p->bm < q->bm ? -1 : 1;
    if (p->r.b != q->r.b) return p->r.b < q->r.b ? -1 : 1;
    if (p->r.a != q->r.a) return p->r.a < q->r.a ? -1 : 1;
    /* (b, a) is unique within a member unless axis a may overlap; keep that case
     * deterministic */
    if (p->r.len != q->r.len) return p->r.len < q->r.len ? -1 : 1;
    return (int)p->r.strand - (int)q->r.strand;
}

/* Order b (SPEC 1.3): a whole axis-a member is buffered, sorted by (b_member, b, a), and
 * cut into chunks ALONG AXIS B -- tui's layout.  Cutting in axis-a order and sorting only
 * inside each chunk is not equivalent: on a universal column axis the span cap then fires
 * at nearly every rearrangement, and a real .tui came out 22x larger that way. */
static int flush_member_b(lmap_writer *w) {
    uint32_t n = w->n_bbuf;
    if (n == 0) return 0;
    qsort(w->bbuf, n, sizeof *w->bbuf, cmp_brun);
    lmap_run *tmp = malloc((size_t)(n < w->count ? n : w->count) * sizeof *tmp);
    if (!tmp) return -1;
    uint32_t i = 0;
    int rc = 0;
    while (i < n && rc == 0) {
        uint32_t bm = w->bbuf[i].bm;
        int64_t blo = w->bbuf[i].r.b, bhi = blo + w->bbuf[i].r.len;
        int64_t alo = w->bbuf[i].r.a, ahi = alo + w->bbuf[i].r.len;
        uint32_t j = i + 1;
        while (j < n && j - i < w->count && w->bbuf[j].bm == bm) {
            const lmap_run *r = &w->bbuf[j].r;
            int64_t nbhi = r->b + r->len > bhi ? r->b + r->len : bhi;     /* sorted: blo fixed */
            int64_t nalo = r->a < alo ? r->a : alo;
            int64_t nahi = r->a + r->len > ahi ? r->a + r->len : ahi;
            if (w->bspan && (uint64_t)(nbhi - blo) > w->bspan) break;
            if ((uint64_t)(nbhi - blo) > 0xFFFFFFFFull || (uint64_t)(nahi - nalo) > 0xFFFFFFFFull) break;
            bhi = nbhi; alo = nalo; ahi = nahi; j++;
        }
        for (uint32_t k = i; k < j; k++) tmp[k - i] = w->bbuf[k].r;
        rc = emit_chunk(w, tmp, j - i, w->bbuf_amem, bm);
        i = j;
    }
    free(tmp);
    w->n_bbuf = 0;
    return rc;
}

int lmap_writer_add_run(lmap_writer *w, uint32_t a_member, uint32_t b_member,
                        int64_t a, int64_t b, int64_t len, uint8_t strand) {
    if (!w || w->failed) return -1;
    if (a < 0 || b < 0 || len < 1) return -1;
    if (a_member >= w->n_mem[0] || b_member >= w->n_mem[1]) return -1;
    if (len > INT64_MAX - a || len > INT64_MAX - b) return -1;
    /* A run must lie inside both members (SPEC 1.2).  The reader rejects a chunk
     * whose extent passes its member's end, so accepting one here would write a
     * file that cannot be opened. */
    if ((uint64_t)(a + len) > w->mem[0][a_member].length ||
        (uint64_t)(b + len) > w->mem[1][b_member].length) return -1;

    /* Axis a is a partial function (SPEC 1.2): a member's runs must not overlap or
     * go backwards, and that has to hold across chunks, not just inside one --
     * otherwise two chunks of the same member can claim the same bases. */
    wmember *am = &w->mem[0][a_member];
    if (am->has_runs && a < am->last_a_exit && !w->a_overlap) { w->failed = 1; return -1; }

    if (w->order == LMAP_ORDER_B) {
        if (w->n_bbuf && a_member != w->bbuf_amem) {
            if (flush_member_b(w) != 0) { w->failed = 1; return -1; }
        }
        if (w->n_bbuf == w->cap_bbuf) {
            uint32_t cap = w->cap_bbuf ? w->cap_bbuf * 2 : 4096;
            struct brun *nb = realloc(w->bbuf, (size_t)cap * sizeof *nb);
            if (!nb) { w->failed = 1; return -1; }
            w->bbuf = nb; w->cap_bbuf = cap;
        }
        struct brun *br = &w->bbuf[w->n_bbuf++];
        br->r.a = a; br->r.b = b; br->r.len = len; br->r.strand = strand ? 1 : 0;
        br->bm = b_member;
        w->bbuf_amem = a_member;
        am->last_a_exit = a + len; am->has_runs = 1;
        return 0;
    }

    uint64_t nlo = (uint64_t)b, nhi = (uint64_t)(b + len);
    if (w->n_pend) {
        if (nlo > w->pend_bmin) nlo = w->pend_bmin;
        if (nhi < w->pend_bmax) nhi = w->pend_bmax;
    }
    /* a_span and b_span are u32 on disk (SPEC 2.2).  A chunk whose extent would not
     * fit is closed early rather than silently wrapping, because a loose or wrapped
     * extent makes every reverse-lookup skip wrong. */
    uint64_t na_end = (uint64_t)(a + len);
    uint64_t na_span = w->n_pend ? na_end - (uint64_t)w->pend[0].a : (uint64_t)len;
    int cut = w->n_pend &&
              (w->n_pend >= w->count ||
               a_member != w->pend_amem || b_member != w->pend_bmem ||
               (w->bspan && (nhi - nlo) > w->bspan) ||
               na_span > 0xFFFFFFFFull || (nhi - nlo) > 0xFFFFFFFFull);
    if (cut) { if (flush_chunk(w) != 0) { w->failed = 1; return -1; } }

    if (w->n_pend == 0) {
        w->pend_amem = a_member; w->pend_bmem = b_member;
        w->pend_bmin = (uint64_t)b; w->pend_bmax = (uint64_t)(b + len);
    } else {
        if ((uint64_t)b < w->pend_bmin) w->pend_bmin = (uint64_t)b;
        if ((uint64_t)(b + len) > w->pend_bmax) w->pend_bmax = (uint64_t)(b + len);
    }
    if (w->n_pend == w->cap_pend) {
        uint32_t cap = w->cap_pend ? w->cap_pend*2 : 1024;
        lmap_run *r = realloc(w->pend, cap * sizeof *r);
        if (!r) { w->failed = 1; return -1; }
        w->pend = r; w->cap_pend = cap;
    }
    lmap_run *r = &w->pend[w->n_pend++];
    r->a = a; r->b = b; r->len = len; r->strand = strand ? 1 : 0;
    am->last_a_exit = a + len; am->has_runs = 1;
    return 0;
}

/* name-sorted permutation over a member array.  The key travels with the element
 * rather than through a global, so two writers can run concurrently. */
typedef struct { const char *name; uint32_t id; } namekey;
static int cmp_by_name(const void *x, const void *y) {
    const namekey *a = (const namekey *)x, *b = (const namekey *)y;
    int c = strcmp(a->name, b->name);
    return c ? c : (a->id < b->id ? -1 : a->id > b->id);
}

static int build_memdir(lmap_writer *w, int axis, vec *out) {
    uint32_t n = w->n_mem[axis];
    vec names = {0};
    uint32_t *noff = calloc(n ? n : 1, sizeof *noff);
    uint32_t *nlen = calloc(n ? n : 1, sizeof *nlen);
    namekey *perm = malloc((n ? n : 1) * sizeof *perm);
    if (!noff || !nlen || !perm) { free(noff); free(nlen); free(perm); return -1; }
    for (uint32_t i = 0; i < n; i++) {
        size_t L = strlen(w->mem[axis][i].name);
        noff[i] = (uint32_t)names.n; nlen[i] = (uint32_t)L;
        if (vec_put(&names, w->mem[axis][i].name, L) != 0) goto fail;
        perm[i].name = w->mem[axis][i].name; perm[i].id = i;
    }
    qsort(perm, n, sizeof *perm, cmp_by_name);

    if (vec_put32(out, n) != 0) goto fail;
    if (vec_put32(out, (uint32_t)names.n) != 0) goto fail;
    for (uint32_t i = 0; i < n; i++) {
        wmember *m = &w->mem[axis][i];
        if (vec_put32(out, noff[i]) || vec_put32(out, nlen[i]) ||
            vec_put64(out, m->length) || vec_put64(out, m->seen ? m->amin : 0) ||
            vec_put64(out, m->seen ? m->amax : 0) ||
            vec_put32(out, m->first_chunk) || vec_put32(out, m->n_chunks)) goto fail;
    }
    if (vec_put(out, names.p, names.n) != 0) goto fail;
    for (uint32_t i = 0; i < n; i++) if (vec_put32(out, perm[i].id) != 0) goto fail;
    vec_free(&names); free(noff); free(nlen); free(perm);
    return 0;
fail:
    vec_free(&names); free(noff); free(nlen); free(perm);
    return -1;
}

/* dir.b: chunk ids ordered by (b_member, b_min), each with the running max of
 * b_min+b_span, so a reverse walk can stop early. */
static int cmp_chunk_a(const void *x, const void *y) {
    const lmap_chunk *p = (const lmap_chunk *)x, *q = (const lmap_chunk *)y;
    if (p->a_member != q->a_member) return p->a_member < q->a_member ? -1 : 1;
    if (p->a_min != q->a_min) return p->a_min < q->a_min ? -1 : 1;
    return 0;
}

/* Record [first, first+n) of a member's entries in a directory, plus its extent. */
static void note_range(wmember *m, uint32_t pos, uint64_t lo, uint64_t hi) {
    if (!m->seen) { m->seen = 1; m->first_chunk = pos; m->amin = lo; m->amax = hi; }
    m->n_chunks++;
    if (lo < m->amin) m->amin = lo;
    if (hi > m->amax) m->amax = hi;
}

typedef struct { uint32_t b_member; uint64_t b_min; uint32_t id; } bkey;
static int cmp_by_b(const void *x, const void *y) {
    const bkey *p = (const bkey *)x, *q = (const bkey *)y;
    if (p->b_member != q->b_member) return p->b_member < q->b_member ? -1 : 1;
    if (p->b_min != q->b_min) return p->b_min < q->b_min ? -1 : 1;
    return p->id < q->id ? -1 : p->id > q->id;
}

static int write_all(FILE *fp, const void *p, size_t n) {
    return (n == 0 || fwrite(p, 1, n, fp) == n) ? 0 : -1;
}

int lmap_writer_close(lmap_writer *w) {
    if (!w) return -1;
    int rc = -1;
    vec foot = {0}, dira = {0}, dirax = {0}, dirb = {0}, ma = {0}, mb = {0};
    if (w->failed) goto done;
    if ((w->order == LMAP_ORDER_B ? flush_member_b(w) : flush_chunk(w)) != 0) goto done;

    /* dir.a is sorted by (a_member, a_min).  Chunks are addressed by payload offset,
     * so reordering the directory moves nothing on disk.  Per-member order across
     * chunks was enforced at add_run, so within a member this is also a_min order
     * with disjoint extents -- which is what makes a binary search over it valid. */
    qsort(w->chunks, w->n_chunks, sizeof *w->chunks, cmp_chunk_a);
    for (int ax = 0; ax < 2; ax++)
        for (uint32_t i = 0; i < w->n_mem[ax]; i++) {
            w->mem[ax][i].seen = 0; w->mem[ax][i].first_chunk = 0; w->mem[ax][i].n_chunks = 0;
        }
    uint64_t a_run_max = 0;
    for (uint32_t i = 0; i < w->n_chunks; i++) {
        lmap_chunk *c = &w->chunks[i];
        note_range(&w->mem[0][c->a_member], i, c->a_min, c->a_min + c->a_span);
        /* dir.a.max: running max of a_end, reset per member -- the axis-a twin of
         * dir.b's.  Under order b a member's chunks overlap on axis a. */
        uint64_t a_end = c->a_min + c->a_span;
        if (i == 0 || c->a_member != w->chunks[i-1].a_member || a_end > a_run_max) a_run_max = a_end;
        if (vec_put64(&dirax, a_run_max) != 0) goto done;
        if (vec_put32(&dira, c->a_member) || vec_put32(&dira, c->b_member) ||
            vec_put64(&dira, c->a_min) || vec_put32(&dira, c->a_span) ||
            vec_put32(&dira, c->b_span) || vec_put64(&dira, c->b_min) ||
            vec_put64(&dira, c->base) || vec_put64(&dira, c->off) ||
            vec_put32(&dira, c->clen) || vec_put32(&dira, c->rawlen) ||
            vec_put32(&dira, c->n_runs) || vec_put32(&dira, c->codec) ||
            vec_put32(&dira, c->crc) || vec_put32(&dira, 0)) goto done;
    }
    {
        /* dir.b: chunk ids sorted by (b_member, b_min), with a running max of b_end
         * that RESETS at each member.  Members are separate coordinate spaces, so a
         * bound carried over from the previous scaffold would be meaningless -- and
         * too large, which silently weakens the early exit for the next one.
         * mem.b ranges index dir.b positions, where a member's chunks are contiguous;
         * they are not contiguous in dir.a. */
        bkey *perm = malloc((w->n_chunks ? w->n_chunks : 1) * sizeof *perm);
        if (!perm) goto done;
        for (uint32_t i = 0; i < w->n_chunks; i++) {
            perm[i].b_member = w->chunks[i].b_member;
            perm[i].b_min = w->chunks[i].b_min;
            perm[i].id = i;
        }
        qsort(perm, w->n_chunks, sizeof *perm, cmp_by_b);
        uint64_t run_max = 0;
        int bad = 0;
        for (uint32_t i = 0; i < w->n_chunks; i++) {
            lmap_chunk *c = &w->chunks[perm[i].id];
            uint64_t end = c->b_min + c->b_span;
            int new_member = (i == 0) || c->b_member != w->chunks[perm[i-1].id].b_member;
            if (new_member || end > run_max) run_max = end;
            note_range(&w->mem[1][c->b_member], i, c->b_min, end);
            bad |= vec_put32(&dirb, perm[i].id);
            bad |= vec_put64(&dirb, run_max);
        }
        free(perm);
        if (bad) goto done;
    }
    if (build_memdir(w, 0, &ma) != 0) goto done;
    if (build_memdir(w, 1, &mb) != 0) goto done;

    /* --- header --- */
    uint32_t hdr_crc;
    {
        uint8_t hdr[LMAP_HEADER_SIZE];
        memset(hdr, 0, sizeof hdr);
        memcpy(hdr, LMAP_MAGIC, 8);
        st32(hdr + 8, LMAP_FORMAT_MAJOR); st32(hdr + 12, 0);
        st64(hdr + 16, (w->order == LMAP_ORDER_B ? LMAP_FEAT_ORDER_B : 0) |
                       (w->a_overlap ? LMAP_FEAT_A_OVERLAP : 0));
        strncpy((char *)hdr + 24, w->profile, 31);
        hdr_crc = crc_all(hdr, sizeof hdr);
        if (write_all(w->fp, hdr, sizeof hdr) != 0) goto done;
    }

    /* --- sections --- */
    typedef struct { const char *id; const uint8_t *p; size_t n; } secdesc;
    const secdesc builtin[] = {
        { "schema",    (const uint8_t *)w->schema, strlen(w->schema) },
        { "mem.a",     ma.p,      ma.n      },
        { "mem.b",     mb.p,      mb.n      },
        { "runs",      w->runs.p, w->runs.n },
        { "dir.a",     dira.p,    dira.n    },
        { "dir.a.max", dirax.p,   dirax.n   },
        { "dir.b",     dirb.p,    dirb.n    },
    };
    const uint32_t NB = (uint32_t)(sizeof builtin / sizeof builtin[0]);
    const uint32_t NSEC = NB + w->n_xs;
    secdesc  *sec = malloc(NSEC * sizeof *sec);
    uint64_t *off = malloc(NSEC * sizeof *off), *len = malloc(NSEC * sizeof *len);
    uint32_t *crc = malloc(NSEC * sizeof *crc);
    if (!sec || !off || !len || !crc) { free(sec); free(off); free(len); free(crc); goto done; }
    for (uint32_t i = 0; i < NB; i++) sec[i] = builtin[i];
    for (uint32_t i = 0; i < w->n_xs; i++) {
        sec[NB + i].id = w->xs[i].id; sec[NB + i].p = w->xs[i].p; sec[NB + i].n = w->xs[i].n;
    }
    uint64_t cur = LMAP_HEADER_SIZE;
    int wbad = 0;
    for (uint32_t i = 0; i < NSEC && !wbad; i++) {
        uint8_t pad[8] = {0};
        size_t padn = (size_t)((8 - (cur % 8)) % 8);
        if (write_all(w->fp, pad, padn) != 0) { wbad = 1; break; }
        cur += padn;
        off[i] = cur; len[i] = sec[i].n;
        crc[i] = crc_all(sec[i].p, (uint64_t)sec[i].n);
        if (write_all(w->fp, sec[i].p, sec[i].n) != 0) { wbad = 1; break; }
        cur += sec[i].n;
    }

    /* --- footer --- */
    if (!wbad && vec_put32(&foot, NSEC) != 0) wbad = 1;
    for (uint32_t i = 0; i < NSEC && !wbad; i++) {
        uint8_t idlen = (uint8_t)strlen(sec[i].id);
        uint8_t cb = LMAP_CODEC_NONE, fl = 0;
        if (vec_put(&foot, &idlen, 1) || vec_put(&foot, sec[i].id, idlen) ||
            vec_put64(&foot, off[i]) || vec_put64(&foot, len[i]) ||
            vec_put64(&foot, (uint64_t)sec[i].n) ||
            vec_put(&foot, &cb, 1) || vec_put(&foot, &fl, 1) ||
            vec_put32(&foot, crc[i])) wbad = 1;
    }
    free(sec); free(off); free(len); free(crc);
    if (wbad) goto done;
    {
        uint8_t pad[8] = {0};
        size_t padn = (size_t)((8 - (cur % 8)) % 8);
        if (write_all(w->fp, pad, padn) != 0) goto done;
        cur += padn;
        uint64_t foot_off = cur;
        if (write_all(w->fp, foot.p, foot.n) != 0) goto done;
        cur += foot.n;

        uint8_t tr[LMAP_TRAILER_SIZE];
        memset(tr, 0, sizeof tr);
        st64(tr, foot_off);
        st64(tr + 8, (uint64_t)foot.n);
        st32(tr + 16, crc_all(foot.p, (uint64_t)foot.n));
        st32(tr + 20, hdr_crc);       /* the header changes how every chunk decodes */
        memcpy(tr + 24, LMAP_MAGIC, 8);
        if (write_all(w->fp, tr, sizeof tr) != 0) goto done;
    }
    if (fflush(w->fp) != 0) goto done;
    rc = 0;
done:
    vec_free(&foot); vec_free(&dira); vec_free(&dirax); vec_free(&dirb); vec_free(&ma); vec_free(&mb);
    if (w->fp) { if (fclose(w->fp) != 0) rc = -1; w->fp = NULL; }
    if (rc != 0 && w->path) remove(w->path);   /* no half-written file left behind */
    lmap_writer_abort(w);
    return rc;
}

void lmap_writer_abort(lmap_writer *w) {
    if (!w) return;
    if (w->fp) fclose(w->fp);
    for (int ax = 0; ax < 2; ax++) {
        for (uint32_t i = 0; i < w->n_mem[ax]; i++) free(w->mem[ax][i].name);
        free(w->mem[ax]);
    }
    free(w->pend); free(w->chunks); free(w->bbuf);
    for (uint32_t i = 0; i < w->n_xs; i++) { free(w->xs[i].id); free(w->xs[i].p); }
    free(w->xs);
    vec_free(&w->runs);
    free(w->path); free(w->profile); free(w->schema);
    free(w);
}

/* ================================================================== READER */

struct lmap_file {
    lmap_io *io; int own_io;
    char profile[32];
    char *schema;
    int order, a_overlap;
    lmap_chunk *chunks; uint32_t n_chunks;
    uint32_t *dirb_id; uint64_t *dirb_max;
    uint64_t *dira_max;
    lmap_member *mem[2]; uint32_t n_mem[2];
    uint32_t *by_name[2];
    uint64_t runs_off, runs_len; uint32_t runs_crc;
    struct xent { char *id; uint64_t off, len; uint32_t crc; } *xs; uint32_t n_xs;
};

static int read_section_body(lmap_io *io, uint64_t off, uint64_t len, uint8_t **out) {
    uint8_t *p = malloc(len ? (size_t)len : 1);
    if (!p) return -1;
    if (lmap_pread(io, p, (int64_t)off, (int64_t)len) != 0) { free(p); return -1; }
    *out = p; return 0;
}

static int parse_memdir(lmap_file *f, int axis, const uint8_t *p, size_t n) {
    if (n < 8) return -1;
    uint32_t cnt = ld32(p), blob = ld32(p + 4);
    /* 64-bit throughout: on ILP32 a size_t product here wraps, and cnt/blob are
     * both read straight off disk. */
    uint64_t need = 8ull + (uint64_t)cnt * LMAP_MEM_ENTRY + (uint64_t)blob + (uint64_t)cnt * 4ull;
    if (need > (uint64_t)n) return -1;
    const uint8_t *e = p + 8;
    const uint8_t *names = e + (size_t)cnt * LMAP_MEM_ENTRY;
    f->mem[axis] = calloc(cnt ? cnt : 1, sizeof *f->mem[axis]);
    f->by_name[axis] = calloc(cnt ? cnt : 1, sizeof **f->by_name);
    if (!f->mem[axis] || !f->by_name[axis]) return -1;
    f->n_mem[axis] = cnt;   /* set first: lmap_close must free whatever we allocated */
    for (uint32_t i = 0; i < cnt; i++, e += LMAP_MEM_ENTRY) {
        uint32_t noff = ld32(e), nlen = ld32(e + 4);
        if ((uint64_t)noff + (uint64_t)nlen > (uint64_t)blob) return -1;
        lmap_member *m = &f->mem[axis][i];
        m->name = malloc((size_t)nlen + 1);     /* nlen+1u wraps to 0 at UINT32_MAX */
        if (!m->name) return -1;
        memcpy(m->name, names + noff, nlen); m->name[nlen] = 0;
        m->length   = ld64(e + 8);
        m->axis_min = ld64(e + 16);
        m->axis_max = ld64(e + 24);
        m->first_chunk = ld32(e + 32);
        m->n_chunks    = ld32(e + 36);
    }
    const uint8_t *perm = names + blob;
    for (uint32_t i = 0; i < cnt; i++) {
        uint32_t id = ld32(perm + (size_t)i * 4);
        if (id >= cnt) return -1;
        f->by_name[axis][i] = id;
    }
    return 0;
}

lmap_file *lmap_open_io(lmap_io *io, int own_io) {
    lmap_file *f = calloc(1, sizeof *f);
    if (!f) return NULL;
    f->io = io; f->own_io = own_io;

    int64_t sz = lmap_io_size(io);
    uint8_t hdr[LMAP_HEADER_SIZE], tr[LMAP_TRAILER_SIZE];
    if (sz < LMAP_HEADER_SIZE + LMAP_TRAILER_SIZE) goto fail;
    if (lmap_pread(io, hdr, 0, LMAP_HEADER_SIZE) != 0) goto fail;
    if (memcmp(hdr, LMAP_MAGIC, 8) != 0) goto fail;
    if (lmap_pread(io, tr, sz - LMAP_TRAILER_SIZE, LMAP_TRAILER_SIZE) != 0) goto fail;
    if (memcmp(tr + 24, LMAP_MAGIC, 8) != 0) goto fail;
    /* Verify the header before trusting any field in it.  The feature bits decide how
     * every chunk decodes: a single flipped bit used to make some chunks decode, without
     * complaint, into different runs. */
    if (crc_all(hdr, LMAP_HEADER_SIZE) != ld32(tr + 20)) goto fail;
    if (ld32(hdr + 8) != LMAP_FORMAT_MAJOR) goto fail;       /* exact: no older reader */
    uint64_t feat = ld64(hdr + 16);
    if (feat & ~LMAP_FEAT_KNOWN) goto fail;                 /* a feature we don't know */
    f->order = (feat & LMAP_FEAT_ORDER_B) ? LMAP_ORDER_B : LMAP_ORDER_A;
    f->a_overlap = (feat & LMAP_FEAT_A_OVERLAP) ? 1 : 0;
    if (f->a_overlap && f->order != LMAP_ORDER_B) goto fail;   /* not encodable in order a */
    memcpy(f->profile, hdr + 24, 31); f->profile[31] = 0;
    uint64_t foff = ld64(tr), flen = ld64(tr + 8);
    if (foff < LMAP_HEADER_SIZE || flen > (uint64_t)sz || foff > (uint64_t)sz - flen) goto fail;

    uint8_t *foot = NULL;
    if (read_section_body(io, foff, flen, &foot) != 0) goto fail;
    if (crc_all(foot, flen) != ld32(tr + 16)) { free(foot); goto fail; }

    /* section table */
    const uint8_t *p = foot, *end = foot + flen;
    if (end - p < 4) { free(foot); goto fail; }
    uint32_t nsec = ld32(p); p += 4;
    /* Smallest possible entry is 1 id byte + 8 + 8 + 8 + 1 + 1 + 4 = 31 bytes.
     * Reject an inflated count before allocating anything for it. */
    if ((uint64_t)nsec * 31ull > flen) goto footfail;
    uint8_t *dira = NULL, *dirax = NULL, *dirb = NULL, *ma = NULL, *mb = NULL;
    uint64_t dira_n = 0, dirax_n = 0, dirb_n = 0, ma_n = 0, mb_n = 0;
    int seen_runs = 0;
    for (uint32_t i = 0; i < nsec; i++) {
        if (end - p < 1) goto footfail;
        uint8_t idlen = *p++;
        if (end - p < idlen + 8 + 8 + 8 + 1 + 1 + 4) goto footfail;
        char id[256]; memcpy(id, p, idlen); id[idlen] = 0; p += idlen;
        uint64_t so = ld64(p); p += 8;
        uint64_t sl = ld64(p); p += 8;
        uint64_t sraw = ld64(p); p += 8;
        uint8_t scodec = *p++; uint8_t sflags = *p++;
        uint32_t scrc = ld32(p); p += 4;
        if (so > (uint64_t)sz || sl > (uint64_t)sz - so) goto footfail;
        if (scodec != LMAP_CODEC_NONE || sflags != 0 || sraw != sl) goto footfail;

        uint8_t **dst = NULL; uint64_t *dstn = NULL;
        if      (!strcmp(id, "dir.a")) { dst = &dira; dstn = &dira_n; }
        else if (!strcmp(id, "dir.a.max")) { dst = &dirax; dstn = &dirax_n; }
        else if (!strcmp(id, "dir.b")) { dst = &dirb; dstn = &dirb_n; }
        else if (!strcmp(id, "mem.a")) { dst = &ma;   dstn = &ma_n;   }
        else if (!strcmp(id, "mem.b")) { dst = &mb;   dstn = &mb_n;   }
        else if (!strcmp(id, "runs")) {
            if (seen_runs) goto footfail;         /* duplicate ids are not legal */
            seen_runs = 1;
            /* Not read here: the payload is the one section that scales with the data,
             * and an index is opened to touch a few chunks.  Each chunk carries its own
             * CRC, checked on read; lmap_verify checks this section's. */
            f->runs_off = so; f->runs_len = sl; f->runs_crc = scrc;
            continue;
        }
        else if (!strcmp(id, "schema")) {
            if (f->schema) goto footfail;
            uint8_t *sb = NULL;
            if (read_section_body(io, so, sl, &sb) != 0) goto footfail;
            if (crc_all(sb, sl) != scrc) { free(sb); goto footfail; }
            f->schema = malloc((size_t)sl + 1);
            if (!f->schema) { free(sb); goto footfail; }
            memcpy(f->schema, sb, (size_t)sl); f->schema[sl] = 0; free(sb);
            continue;
        } else if (!strncmp(id, "x.", 2)) {       /* application section: remember, load on demand */
            for (uint32_t k = 0; k < f->n_xs; k++) if (!strcmp(f->xs[k].id, id)) goto footfail;
            struct xent *nx = realloc(f->xs, (f->n_xs + 1) * sizeof *nx);
            if (!nx) goto footfail;
            f->xs = nx;
            f->xs[f->n_xs].id = strdup(id);
            if (!f->xs[f->n_xs].id) goto footfail;
            f->xs[f->n_xs].off = so; f->xs[f->n_xs].len = sl; f->xs[f->n_xs].crc = scrc;
            f->n_xs++;
            continue;
        } else continue;                          /* unknown section: skip, not an error */

        if (*dst) goto footfail;                  /* duplicate id: would leak the first */
        if (read_section_body(io, so, sl, dst) != 0) goto footfail;
        *dstn = sl;
        if (crc_all(*dst, sl) != scrc) goto footfail;
    }
    free(foot); foot = NULL;

    if (!dira || dira_n % LMAP_DIRA_ENTRY) goto tablefail;
    if (!dirax || dirax_n != (dira_n / LMAP_DIRA_ENTRY) * 8) goto tablefail;
    f->n_chunks = (uint32_t)(dira_n / LMAP_DIRA_ENTRY);
    f->chunks = calloc(f->n_chunks ? f->n_chunks : 1, sizeof *f->chunks);
    if (!f->chunks) goto tablefail;
    for (uint32_t i = 0; i < f->n_chunks; i++) {
        const uint8_t *e = dira + (size_t)i * LMAP_DIRA_ENTRY;
        lmap_chunk *c = &f->chunks[i];
        c->a_member = ld32(e);      c->b_member = ld32(e + 4);
        c->a_min    = ld64(e + 8);  c->a_span   = ld32(e + 16);
        c->b_span   = ld32(e + 20); c->b_min    = ld64(e + 24);
        c->base     = ld64(e + 32); c->off      = ld64(e + 40);
        c->clen     = ld32(e + 48); c->rawlen   = ld32(e + 52);
        c->n_runs   = ld32(e + 56); c->codec = ld32(e + 60);
        c->crc      = ld32(e + 64);
        if (c->codec != LMAP_CODEC_NONE && c->codec != LMAP_CODEC_DEFLATE) goto tablefail;
        if (ld32(e + 68) != 0) goto tablefail;              /* reserved */
        /* n_runs sizes the decode buffer, so it must be answerable to the bytes that
         * are actually there: every run costs at least one varint byte in each of the
         * a, b and len streams plus a strand bit, after the 16-byte stream header.
         * rawlen in turn must be reachable from clen. */
        if (c->n_runs == 0 || c->n_runs > LMAP_MAX_CHUNK_RUNS) goto tablefail;
        if (c->rawlen > LMAP_MAX_CHUNK_RAW) goto tablefail;
        if ((uint64_t)c->rawlen < 16ull + 3ull * c->n_runs + (c->n_runs + 7ull) / 8) goto tablefail;
        /* clen = 32 + stored bytes, rawlen = 16 + raw bytes, and stored <= raw per stream */
        if (c->clen < 32) goto tablefail;
        if (c->codec == LMAP_CODEC_NONE ? (uint64_t)c->clen != (uint64_t)c->rawlen + 16
            : ((uint64_t)c->clen > (uint64_t)c->rawlen + 16 ||
               (uint64_t)c->rawlen > (uint64_t)LMAP_ZLIB_MAX_RATIO * c->clen + 64)) goto tablefail;
        if (c->off > f->runs_len || c->clen > f->runs_len - c->off) goto tablefail;
        /* a_min, b_min and base become int64 decode bases. */
        if (c->a_min > (uint64_t)INT64_MAX || c->b_min > (uint64_t)INT64_MAX ||
            c->base > (uint64_t)INT64_MAX) goto tablefail;
    }
    f->dira_max = malloc((f->n_chunks ? f->n_chunks : 1) * sizeof *f->dira_max);
    if (!f->dira_max) goto tablefail;
    for (uint32_t i = 0; i < f->n_chunks; i++) f->dira_max[i] = ld64(dirax + (size_t)i * 8);
    if (dirb && dirb_n == (uint64_t)f->n_chunks * LMAP_DIRB_ENTRY) {
        f->dirb_id  = malloc((f->n_chunks ? f->n_chunks : 1) * sizeof *f->dirb_id);
        f->dirb_max = malloc((f->n_chunks ? f->n_chunks : 1) * sizeof *f->dirb_max);
        if (!f->dirb_id || !f->dirb_max) goto tablefail;
        for (uint32_t i = 0; i < f->n_chunks; i++) {
            const uint8_t *e = dirb + (size_t)i * LMAP_DIRB_ENTRY;
            uint32_t id = ld32(e);
            if (id >= f->n_chunks) goto tablefail;
            f->dirb_id[i] = id; f->dirb_max[i] = ld64(e + 4);
        }
    }
    if (ma && parse_memdir(f, 0, ma, (size_t)ma_n) != 0) goto tablefail;
    if (mb && parse_memdir(f, 1, mb, (size_t)mb_n) != 0) goto tablefail;
    /* Every member id we hand to a caller must actually name a member. */
    for (uint32_t i = 0; i < f->n_chunks; i++) {
        if (f->chunks[i].a_member >= f->n_mem[0] ||
            f->chunks[i].b_member >= f->n_mem[1]) goto tablefail;
    }
    /* The invariants the queries rely on, checked rather than trusted.  A file that
     * violates any of them would make a binary search or an early exit silently
     * return the wrong answer, which is worse than refusing to open it. */
    for (int ax = 0; ax < 2; ax++) {
        uint64_t covered = 0;
        for (uint32_t i = 0; i < f->n_mem[ax]; i++) {
            lmap_member *m = &f->mem[ax][i];
            if (m->n_chunks > f->n_chunks ||
                m->first_chunk > f->n_chunks - m->n_chunks) goto tablefail;
            covered += m->n_chunks;
        }
        if (f->n_chunks && covered != f->n_chunks) goto tablefail;
    }
    /* mem.a ranges index dir.a: one member each, a_min increasing, extents disjoint */
    for (uint32_t mi = 0; mi < f->n_mem[0]; mi++) {
        const lmap_member *m = &f->mem[0][mi];
        uint64_t run_max = 0;
        for (uint32_t k = 0; k < m->n_chunks; k++) {
            uint32_t pos = m->first_chunk + k;
            const lmap_chunk *c = &f->chunks[pos];
            uint64_t end = c->a_min + c->a_span;
            if (c->a_member != mi) goto tablefail;
            if (k) {
                const lmap_chunk *pc = &f->chunks[pos - 1];
                if (c->a_min < pc->a_min) goto tablefail;
                /* order a cuts along a, so a member's chunks cannot overlap there */
                if (f->order == LMAP_ORDER_A && c->a_min < pc->a_min + pc->a_span) goto tablefail;
            }
            if (end > m->length) goto tablefail;
            if (k == 0 || end > run_max) run_max = end;
            if (f->dira_max[pos] != run_max) goto tablefail;
        }
    }
    if (f->dirb_id) {
        /* dir.b must be a permutation of the chunk ids */
        uint8_t *hit = calloc(f->n_chunks ? f->n_chunks : 1, 1);
        if (!hit) goto tablefail;
        for (uint32_t i = 0; i < f->n_chunks; i++) {
            if (hit[f->dirb_id[i]]) { free(hit); goto tablefail; }
            hit[f->dirb_id[i]] = 1;
        }
        free(hit);
        /* mem.b ranges index dir.b: one member each, b_min non-decreasing, and the
         * stored prefix max must equal the recomputed one, reset per member */
        for (uint32_t mi = 0; mi < f->n_mem[1]; mi++) {
            const lmap_member *m = &f->mem[1][mi];
            uint64_t run_max = 0;
            for (uint32_t k = 0; k < m->n_chunks; k++) {
                uint32_t pos = m->first_chunk + k;
                const lmap_chunk *c = &f->chunks[f->dirb_id[pos]];
                if (c->b_member != mi) goto tablefail;
                uint64_t end = c->b_min + c->b_span;
                if (end > m->length) goto tablefail;
                if (k && c->b_min < f->chunks[f->dirb_id[pos - 1]].b_min) goto tablefail;
                if (k == 0 || end > run_max) run_max = end;
                if (f->dirb_max[pos] != run_max) goto tablefail;
            }
        }
    }
    free(dira); free(dirax); free(dirb); free(ma); free(mb);
    return f;

footfail:
    free(foot);
tablefail:
    free(dira); free(dirax); free(dirb); free(ma); free(mb);
fail:
    lmap_close(f);      /* honours own_io, so a failed lmap_open_io(io,1) frees io */
    return NULL;
}

lmap_file *lmap_open(const char *path) {
    lmap_io *io = lmap_io_open_file(path);
    if (!io) return NULL;
    return lmap_open_io(io, 1);   /* takes ownership on success and on failure */
}

void lmap_close(lmap_file *f) {
    if (!f) return;
    for (int ax = 0; ax < 2; ax++) {
        for (uint32_t i = 0; i < f->n_mem[ax]; i++) free(f->mem[ax][i].name);
        free(f->mem[ax]); free(f->by_name[ax]);
    }
    free(f->chunks); free(f->dirb_id); free(f->dirb_max); free(f->dira_max); free(f->schema);
    for (uint32_t i = 0; i < f->n_xs; i++) free(f->xs[i].id);
    free(f->xs);
    if (f->own_io) lmap_io_close(f->io);
    free(f);
}

const char *lmap_profile(const lmap_file *f)     { return f ? f->profile : NULL; }
const char *lmap_schema_text(const lmap_file *f) { return f ? f->schema : NULL; }
int         lmap_order(const lmap_file *f)       { return f ? f->order : -1; }
int         lmap_a_overlap(const lmap_file *f)   { return f ? f->a_overlap : 0; }
uint32_t    lmap_n_chunks(const lmap_file *f)    { return f ? f->n_chunks : 0; }
uint32_t    lmap_n_members(const lmap_file *f, int axis) {
    return (f && (axis == 0 || axis == 1)) ? f->n_mem[axis] : 0;
}
const lmap_member *lmap_member_at(const lmap_file *f, int axis, uint32_t i) {
    if (!f || (axis != 0 && axis != 1) || i >= f->n_mem[axis]) return NULL;
    return &f->mem[axis][i];
}
const lmap_chunk *lmap_chunk_at(const lmap_file *f, uint32_t i) {
    return (f && i < f->n_chunks) ? &f->chunks[i] : NULL;
}

int32_t lmap_member_by_name(const lmap_file *f, int axis, const char *name) {
    if (!f || (axis != 0 && axis != 1) || !name || !f->by_name[axis]) return -1;
    uint32_t lo = 0, hi = f->n_mem[axis];
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        uint32_t id = f->by_name[axis][mid];
        int c = strcmp(f->mem[axis][id].name, name);
        if (c == 0) return (int32_t)id;
        if (c < 0) lo = mid + 1; else hi = mid;
    }
    return -1;
}

/* Decode a chunk's streams and check every run lies inside the chunk's declared
 * extents.  Containment is what makes skipping a chunk on its extent sound: a run
 * outside the declared range could be missed by a query that skipped the chunk. */
static int decode_checked(const lmap_file *f, const lmap_chunk *c,
                          const lmap_buf st[LMAP_N_STREAMS], lmap_run *out) {
    int64_t base_a = f->order == LMAP_ORDER_A ? (int64_t)c->a_min : (int64_t)c->base;
    int64_t base_b = f->order == LMAP_ORDER_A ? (int64_t)c->base  : (int64_t)c->b_min;
    if (lmap_decode_chunk(st, c->n_runs, f->order, base_a, base_b, out) != 0) return -1;
    uint64_t a_hi = c->a_min + c->a_span, b_hi = c->b_min + c->b_span;
    for (uint32_t k = 0; k < c->n_runs; k++) {
        const lmap_run *r = &out[k];
        if ((uint64_t)r->a < c->a_min || (uint64_t)(r->a + r->len) > a_hi ||
            (uint64_t)r->b < c->b_min || (uint64_t)(r->b + r->len) > b_hi) return -1;
    }
    return 0;
}

int lmap_read_chunk(lmap_file *f, uint32_t i, lmap_run *out) {
    const lmap_chunk *c = lmap_chunk_at(f, i);
    if (!c || !out || c->clen < 32) return -1;
    uint8_t *raw = malloc(c->clen);
    if (!raw) return -1;
    if (lmap_pread(f->io, raw, (int64_t)(f->runs_off + c->off), (int64_t)c->clen) != 0 ||
        crc_all(raw, c->clen) != c->crc) {
        free(raw); return -1;
    }
    /* Blob (SPEC 2.2): u32 raw_len[4], u32 stored_len[4], the four stored streams.
     * Every length is off disk; each is checked against the directory entry, which was
     * itself bounded at open, before anything is allocated from it. */
    uint32_t rl[LMAP_N_STREAMS], sl[LMAP_N_STREAMS];
    uint64_t rsum = 0, ssum = 0;
    int ok = 1;
    for (int k = 0; k < LMAP_N_STREAMS; k++) {
        rl[k] = ld32(raw + 4*k);
        sl[k] = ld32(raw + 16 + 4*k);
        rsum += rl[k]; ssum += sl[k];
        if (sl[k] > rl[k]) ok = 0;                            /* never stored larger */
        if (c->codec == LMAP_CODEC_NONE && sl[k] != rl[k]) ok = 0;
    }
    if (!ok || rsum + 16 != c->rawlen || ssum + 32 != c->clen) { free(raw); return -1; }

    uint8_t *inf = malloc(rsum ? (size_t)rsum : 1);
    if (!inf) { free(raw); return -1; }
    lmap_buf st[LMAP_N_STREAMS];
    const uint8_t *q = raw + 32;
    uint8_t *o = inf;
    for (int k = 0; ok && k < LMAP_N_STREAMS; k++) {
        if (sl[k] == rl[k]) memcpy(o, q, rl[k]);              /* stored verbatim */
        else if (raw_inflate(q, sl[k], o, rl[k]) != 0) ok = 0;
        st[k].p = o; st[k].n = rl[k]; st[k].cap = rl[k];
        q += sl[k]; o += rl[k];
    }
    int rc = ok ? decode_checked(f, c, st, out) : -1;
    free(raw); free(inf);
    return rc;
}

/* ================================================================== QUERIES */

/* Clip run r to [lo,hi) on `axis`, mapping the other axis through the run's strand.
 * Reverse strand pairs a+i with b+len-1-i, so a window [o1,o2) of offsets on one
 * axis is [len-o2, len-o1) on the other. */
static int clip_run(const lmap_run *r, int axis, int64_t lo, int64_t hi, lmap_hit *h) {
    int64_t s0 = axis == 0 ? r->a : r->b;
    int64_t s = s0 > lo ? s0 : lo;
    int64_t e = (s0 + r->len) < hi ? (s0 + r->len) : hi;
    if (s >= e) return 0;
    int64_t o1 = s - s0, o2 = e - s0;
    int64_t other0 = axis == 0 ? r->b : r->a;
    int64_t other = r->strand ? other0 + r->len - o2 : other0 + o1;
    if (axis == 0) { h->a = s; h->b = other; } else { h->b = s; h->a = other; }
    h->len = e - s;
    h->strand = r->strand;
    return 1;
}

typedef struct { lmap_hit *p; size_t n, cap; } hitvec;
static int hit_push(hitvec *v, const lmap_hit *h) {
    if (v->n == v->cap) {
        size_t cap = v->cap ? v->cap * 2 : 64;
        lmap_hit *q = realloc(v->p, cap * sizeof *q);
        if (!q) return -1;
        v->p = q; v->cap = cap;
    }
    v->p[v->n++] = *h;
    return 0;
}

/* Decode chunk ci and append its runs' overlaps with [lo,hi) on `axis`. */
static int cmp_hit_a(const void *x, const void *y) {
    const lmap_hit *p = (const lmap_hit *)x, *q = (const lmap_hit *)y;
    /* a alone is unique unless the file declares A_OVERLAP; the rest keeps that case
     * deterministic */
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    if (p->len != q->len) return p->len < q->len ? -1 : 1;
    return (int)p->strand - (int)q->strand;
}

/* The one-shot queries: a streaming cursor over the window, clipped, in their
 * documented orders. */
static int query_via_cursor(lmap_file *f, int axis, uint32_t m, int64_t lo, int64_t hi,
                            lmap_hit **out, size_t *n, lmap_query_stats *st,
                            int (*cmp)(const void *, const void *)) {
    if (!f || !out || !n || m >= f->n_mem[axis] || lo < 0 || hi < lo) return -1;
    *out = NULL; *n = 0;
    if (st) memset(st, 0, sizeof *st);
    lmap_cursor *c = lmap_cursor_open(f, axis, m, NULL, 0, 0);
    if (!c) return -1;
    hitvec v = {0};
    lmap_hit h, t;
    int rc = lmap_cursor_seek(c, lo, hi);
    while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
        rc = 0;
        if (lmap_hit_clip(&h, axis, lo, hi, &t) && hit_push(&v, &t) != 0) rc = -1;
    }
    if (st) {
        lmap_cursor_stats cs; lmap_cursor_get_stats(c, &cs);
        st->chunks_examined = (uint32_t)cs.chunks_examined;
        st->chunks_decoded = (uint32_t)cs.chunks_decoded;
    }
    lmap_cursor_close(c);
    if (rc != 0) { free(v.p); return -1; }
    if (v.n > 1) qsort(v.p, v.n, sizeof *v.p, cmp);
    *out = v.p; *n = v.n;
    return 0;
}

int lmap_query_a(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *st) {
    return query_via_cursor(f, 0, m, lo, hi, out, n, st, cmp_hit_a);
}

static int cmp_hit_b(const void *x, const void *y) {
    const lmap_hit *p = (const lmap_hit *)x, *q = (const lmap_hit *)y;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    if (p->a_member != q->a_member) return p->a_member < q->a_member ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return 0;
}

int lmap_query_b(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *st) {
    if (f && !f->dirb_id) return -1;                    /* no reverse index in this file */
    return query_via_cursor(f, 1, m, lo, hi, out, n, st, cmp_hit_b);
}

/* ================================================================ CURSOR
 *
 * A view is the key member's chunks in key_min order with a running max of key_end --
 * borrowed straight from the file's directory, or, when the other axis is filtered to a
 * set of members, gathered from those members' directory ranges and sorted.  Iteration
 * merges the runs of every chunk overlapping the window with a heap, opening a chunk
 * only once the merge reaches its key_min, so runs come out in key order across chunks
 * that overlap, and at most the overlapping chunks are held at once.  Decoded chunks are
 * cached under a byte budget, least recently used first out; a chunk in use by the
 * iteration is pinned and never evicted. */

typedef struct dchunk {
    uint32_t pos, other_member, key_member;
    lmap_run *runs;          /* sorted by (key, other, len, strand) */
    int64_t  *mep;           /* running max of key end */
    uint32_t  n;
    size_t    bytes;
    int       pins;
    int       partial;       /* holds only the runs of one window: never reused */
    uint32_t  hint;          /* last point lookup's first index, for sweeps */
    struct dchunk *prev, *next, *hnext;
} dchunk;

typedef struct { dchunk *c; uint32_t i; } chead;

struct lmap_cursor {
    lmap_file *f;
    int K;
    uint32_t npos, id_base;
    const uint32_t *ids;     /* chunk id per position; NULL: id_base + pos */
    const uint64_t *pmax;    /* running max of key end per position */
    uint32_t *own_ids; uint64_t *own_pmax;
    dchunk **htab; uint32_t hcap, hn;
    dchunk *lru_head, *lru_tail;
    size_t cached, budget;
    int64_t lo, hi;
    uint32_t next_pos;
    uint32_t point_hint;     /* last point lookup's first position, for sweeps */
    chead *heap; uint32_t nheap, capheap;
    lmap_cursor_stats st;
};

static inline uint32_t cur_cid(const lmap_cursor *c, uint32_t pos) {
    return c->ids ? c->ids[pos] : c->id_base + pos;
}
static inline int64_t cur_kmin(const lmap_cursor *c, uint32_t pos) {
    const lmap_chunk *k = &c->f->chunks[cur_cid(c, pos)];
    return (int64_t)(c->K ? k->b_min : k->a_min);
}
static inline int64_t cur_kend(const lmap_cursor *c, uint32_t pos) {
    const lmap_chunk *k = &c->f->chunks[cur_cid(c, pos)];
    return (int64_t)(c->K ? k->b_min + k->b_span : k->a_min + k->a_span);
}
static inline int64_t rkey(const lmap_run *r, int K)   { return K ? r->b : r->a; }
static inline int64_t rother(const lmap_run *r, int K) { return K ? r->a : r->b; }

/* qsort takes no context, and a global would race between threads: one comparator per
 * key axis instead. */
static int cmp_run_key(const lmap_run *p, const lmap_run *q, int K) {
    if (rkey(p, K) != rkey(q, K)) return rkey(p, K) < rkey(q, K) ? -1 : 1;
    if (rother(p, K) != rother(q, K)) return rother(p, K) < rother(q, K) ? -1 : 1;
    if (p->len != q->len) return p->len < q->len ? -1 : 1;
    return (int)p->strand - (int)q->strand;
}
static int cmp_run_key0(const void *x, const void *y) { return cmp_run_key(x, y, 0); }
static int cmp_run_key1(const void *x, const void *y) { return cmp_run_key(x, y, 1); }

typedef struct { uint64_t key; uint32_t cid; } keyed_id;
static int cmp_keyed_id(const void *x, const void *y) {
    const keyed_id *p = x, *q = y;
    if (p->key != q->key) return p->key < q->key ? -1 : 1;
    return p->cid < q->cid ? -1 : p->cid > q->cid;
}

lmap_cursor *lmap_cursor_open(lmap_file *f, int key_axis, uint32_t key_member,
                              const uint32_t *other, uint32_t n_other, size_t cache_bytes) {
    if (!f || (key_axis != 0 && key_axis != 1) || key_member >= f->n_mem[key_axis]) return NULL;
    if (key_axis == 1 && !f->dirb_id) return NULL;
    lmap_cursor *c = calloc(1, sizeof *c);
    if (!c) return NULL;
    c->f = f; c->K = key_axis; c->budget = cache_bytes;
    const lmap_member *km = &f->mem[key_axis][key_member];
    if (!other) {
        c->npos = km->n_chunks;
        if (key_axis == 0) { c->id_base = km->first_chunk; c->pmax = f->dira_max + km->first_chunk; }
        else { c->ids = f->dirb_id + km->first_chunk; c->pmax = f->dirb_max + km->first_chunk; }
    } else {
        int O = 1 - key_axis;
        uint64_t total = 0;
        for (uint32_t i = 0; i < n_other; i++) {
            if (other[i] >= f->n_mem[O]) { free(c); return NULL; }
            total += f->mem[O][other[i]].n_chunks;
        }
        c->own_ids = malloc((total ? total : 1) * sizeof *c->own_ids);
        c->own_pmax = malloc((total ? total : 1) * sizeof *c->own_pmax);
        if (!c->own_ids || !c->own_pmax) { lmap_cursor_close(c); return NULL; }
        keyed_id *kid = malloc((total ? total : 1) * sizeof *kid);
        if (!kid) { lmap_cursor_close(c); return NULL; }
        uint32_t n = 0;
        for (uint32_t i = 0; i < n_other; i++) {
            const lmap_member *om = &f->mem[O][other[i]];
            for (uint32_t k = 0; k < om->n_chunks; k++) {
                uint32_t pos = om->first_chunk + k;
                uint32_t cid = O == 0 ? pos : f->dirb_id[pos];
                const lmap_chunk *ch = &f->chunks[cid];
                if ((key_axis ? ch->b_member : ch->a_member) != key_member) continue;
                kid[n].key = key_axis ? ch->b_min : ch->a_min; kid[n].cid = cid; n++;
            }
        }
        if (n > 1) qsort(kid, n, sizeof *kid, cmp_keyed_id);
        for (uint32_t i = 0; i < n; i++) c->own_ids[i] = kid[i].cid;
        free(kid);
        c->ids = c->own_ids; c->npos = n;
        uint64_t run = 0;
        for (uint32_t pos = 0; pos < n; pos++) {
            uint64_t e = (uint64_t)cur_kend(c, pos);
            if (pos == 0 || e > run) run = e;
            c->own_pmax[pos] = run;
        }
        c->pmax = c->own_pmax;
    }
    c->hcap = 64;
    c->htab = calloc(c->hcap, sizeof *c->htab);
    if (!c->htab) { lmap_cursor_close(c); return NULL; }
    c->lo = c->hi = 0; c->next_pos = c->npos;
    return c;
}

static void lru_remove(lmap_cursor *c, dchunk *d) {
    if (d->prev) d->prev->next = d->next; else if (c->lru_head == d) c->lru_head = d->next;
    if (d->next) d->next->prev = d->prev; else if (c->lru_tail == d) c->lru_tail = d->prev;
    d->prev = d->next = NULL;
}

static void dchunk_free_unhash(lmap_cursor *c, dchunk *d) {
    dchunk **pp = &c->htab[d->pos & (c->hcap - 1)];
    while (*pp && *pp != d) pp = &(*pp)->hnext;
    if (*pp) *pp = d->hnext;
    c->hn--; c->cached -= d->bytes;
    free(d->runs); free(d->mep); free(d);
}

static void evict(lmap_cursor *c) {
    while (c->cached > c->budget && c->lru_tail) {
        dchunk *d = c->lru_tail;
        lru_remove(c, d);
        dchunk_free_unhash(c, d);
    }
}

static void release(lmap_cursor *c, dchunk *d) {
    if (--d->pins > 0) return;
    d->next = c->lru_head; d->prev = NULL;
    if (c->lru_head) c->lru_head->prev = d; else c->lru_tail = d;
    c->lru_head = d;
    evict(c);
}

static int hash_grow(lmap_cursor *c) {
    uint32_t nc = c->hcap * 2;
    dchunk **nt = calloc(nc, sizeof *nt);
    if (!nt) return -1;
    for (uint32_t i = 0; i < c->hcap; i++)
        for (dchunk *d = c->htab[i], *nx; d; d = nx) { nx = d->hnext; d->hnext = nt[d->pos & (nc - 1)]; nt[d->pos & (nc - 1)] = d; }
    free(c->htab); c->htab = nt; c->hcap = nc;
    return 0;
}

/* The decoded chunk at view position pos, pinned, holding at least the runs overlapping
 * [lo, hi) on the key axis.  NULL on error. */
static dchunk *acquire(lmap_cursor *c, uint32_t pos, int64_t lo, int64_t hi) {
    for (dchunk *d = c->htab[pos & (c->hcap - 1)]; d; d = d->hnext)
        if (d->pos == pos && !d->partial) {
            if (d->pins++ == 0) lru_remove(c, d);
            c->st.chunk_reuses++;
            return d;
        }
    uint32_t cid = cur_cid(c, pos);
    const lmap_chunk *ch = &c->f->chunks[cid];
    dchunk *d = calloc(1, sizeof *d);
    if (!d) return NULL;
    d->pos = pos; d->n = ch->n_runs;
    d->key_member = c->K ? ch->b_member : ch->a_member;
    d->other_member = c->K ? ch->a_member : ch->b_member;
    d->runs = malloc((size_t)d->n * sizeof *d->runs);
    d->mep = malloc((size_t)d->n * sizeof *d->mep);
    if (!d->runs || !d->mep || lmap_read_chunk(c->f, cid, d->runs) != 0) {
        free(d->runs); free(d->mep); free(d); return NULL;
    }
    /* Stored order is sorted on the file's order axis -- by (a) under order a, by
     * (b, a, len, strand) under order b, which is already key order when the key is the
     * order axis.  Only the other case needs sorting. */
    if (c->f->order != c->K) {
        /* A cursor that keeps nothing needs only this window's runs, and sorting a whole
         * chunk to use a few of them dominated tui's forward queries. */
        if (c->budget == 0) {
            uint32_t m = 0;
            for (uint32_t i = 0; i < d->n; i++) {
                int64_t k = rkey(&d->runs[i], c->K);
                if (k < hi && k + d->runs[i].len > lo) d->runs[m++] = d->runs[i];
            }
            d->n = m; d->partial = 1;
        }
        qsort(d->runs, d->n, sizeof *d->runs, c->K ? cmp_run_key1 : cmp_run_key0);
    }
    int64_t run = 0;
    for (uint32_t i = 0; i < d->n; i++) {
        int64_t e = rkey(&d->runs[i], c->K) + d->runs[i].len;
        if (i == 0 || e > run) run = e;
        d->mep[i] = run;
    }
    d->bytes = sizeof *d + (size_t)d->n * (sizeof *d->runs + sizeof *d->mep);
    d->pins = 1;
    if ((c->hn + 1) * 2 > c->hcap && hash_grow(c) != 0) { free(d->runs); free(d->mep); free(d); return NULL; }
    d->hnext = c->htab[pos & (c->hcap - 1)]; c->htab[pos & (c->hcap - 1)] = d;
    c->hn++; c->cached += d->bytes;
    if (c->cached > c->st.peak_cached_bytes) c->st.peak_cached_bytes = c->cached;
    c->st.chunks_decoded++;
    return d;
}

/* first index whose running max exceeds v (both arrays are non-decreasing) */
static uint32_t first_above_u64(const uint64_t *a, uint32_t n, int64_t v) {
    uint32_t lo = 0, hi = n;
    while (lo < hi) { uint32_t m = lo + (hi - lo) / 2; if ((int64_t)a[m] > v) hi = m; else lo = m + 1; }
    return lo;
}
static uint32_t first_above_i64(const int64_t *a, uint32_t n, int64_t v) {
    uint32_t lo = 0, hi = n;
    while (lo < hi) { uint32_t m = lo + (hi - lo) / 2; if (a[m] > v) hi = m; else lo = m + 1; }
    return lo;
}

/* Same, searching forward from `hint` when it is a valid lower bound (a[hint-1] <= v):
 * a sweep of increasing v then finds each answer in O(log distance), not O(log n). */
static uint32_t gallop_u64(const uint64_t *a, uint32_t n, int64_t v, uint32_t hint) {
    if (hint > n || (hint > 0 && (int64_t)a[hint - 1] > v)) return first_above_u64(a, n, v);
    uint32_t lo = hint, step = 1;
    while (lo + step <= n && (int64_t)a[lo + step - 1] <= v) { lo += step; step *= 2; }
    uint32_t hi = lo + step <= n ? lo + step : n;
    return lo + first_above_u64(a + lo, hi - lo, v);
}
static uint32_t gallop_i64(const int64_t *a, uint32_t n, int64_t v, uint32_t hint) {
    if (hint > n || (hint > 0 && a[hint - 1] > v)) return first_above_i64(a, n, v);
    uint32_t lo = hint, step = 1;
    while (lo + step <= n && a[lo + step - 1] <= v) { lo += step; step *= 2; }
    uint32_t hi = lo + step <= n ? lo + step : n;
    return lo + first_above_i64(a + lo, hi - lo, v);
}

/* heap order = output order: key, other member, other position, len, strand, chunk */
static int head_less(const lmap_cursor *c, const chead *x, const chead *y) {
    const lmap_run *p = &x->c->runs[x->i], *q = &y->c->runs[y->i];
    int K = c->K;
    if (rkey(p, K) != rkey(q, K)) return rkey(p, K) < rkey(q, K);
    if (x->c->other_member != y->c->other_member) return x->c->other_member < y->c->other_member;
    if (rother(p, K) != rother(q, K)) return rother(p, K) < rother(q, K);
    if (p->len != q->len) return p->len < q->len;
    if (p->strand != q->strand) return p->strand < q->strand;
    return x->c->pos < y->c->pos;
}
static void heap_up(lmap_cursor *c, uint32_t i) {
    while (i > 0) {
        uint32_t p = (i - 1) / 2;
        if (!head_less(c, &c->heap[i], &c->heap[p])) break;
        chead t = c->heap[i]; c->heap[i] = c->heap[p]; c->heap[p] = t; i = p;
    }
}
static void heap_down(lmap_cursor *c, uint32_t i) {
    for (;;) {
        uint32_t l = 2 * i + 1, r = l + 1, m = i;
        if (l < c->nheap && head_less(c, &c->heap[l], &c->heap[m])) m = l;
        if (r < c->nheap && head_less(c, &c->heap[r], &c->heap[m])) m = r;
        if (m == i) break;
        chead t = c->heap[i]; c->heap[i] = c->heap[m]; c->heap[m] = t; i = m;
    }
}

/* Move h->i to the next run overlapping [lo, hi); 0 if the chunk has none left. */
static int head_valid(const lmap_cursor *c, chead *h) {
    for (; h->i < h->c->n; h->i++) {
        const lmap_run *r = &h->c->runs[h->i];
        int64_t k = rkey(r, c->K);
        if (k >= c->hi) return 0;
        if (k + r->len > c->lo) return 1;
    }
    return 0;
}

static void heap_clear(lmap_cursor *c) {
    for (uint32_t i = 0; i < c->nheap; i++) release(c, c->heap[i].c);
    c->nheap = 0;
}

int lmap_cursor_seek(lmap_cursor *c, int64_t lo, int64_t hi) {
    if (!c || lo < 0 || hi < lo) return -1;
    heap_clear(c);
    c->lo = lo; c->hi = hi;
    c->next_pos = first_above_u64(c->pmax, c->npos, lo);   /* earlier chunks all end <= lo */
    return 0;
}

int lmap_cursor_next(lmap_cursor *c, lmap_hit *out) {
    if (!c || !out) return -1;
    for (;;) {
        /* open every chunk that could hold a run at or before the current minimum */
        while (c->next_pos < c->npos) {
            int64_t km = cur_kmin(c, c->next_pos);
            if (km >= c->hi) { c->next_pos = c->npos; break; }
            if (c->nheap && km > rkey(&c->heap[0].c->runs[c->heap[0].i], c->K)) break;
            uint32_t pos = c->next_pos++;
            c->st.chunks_examined++;
            if (cur_kend(c, pos) <= c->lo) continue;
            dchunk *d = acquire(c, pos, c->lo, c->hi);
            if (!d) return -1;
            chead h = { d, first_above_i64(d->mep, d->n, c->lo) };
            if (!head_valid(c, &h)) { release(c, d); continue; }
            if (c->nheap == c->capheap) {
                uint32_t cap = c->capheap ? c->capheap * 2 : 16;
                chead *nh = realloc(c->heap, cap * sizeof *nh);
                if (!nh) { release(c, d); return -1; }
                c->heap = nh; c->capheap = cap;
            }
            c->heap[c->nheap++] = h;
            heap_up(c, c->nheap - 1);
        }
        if (c->nheap == 0) return 0;
        chead *h = &c->heap[0];
        const lmap_run *r = &h->c->runs[h->i];
        out->a = r->a; out->b = r->b; out->len = r->len; out->strand = r->strand;
        out->a_member = c->K ? h->c->other_member : h->c->key_member;
        out->b_member = c->K ? h->c->key_member : h->c->other_member;
        h->i++;
        if (head_valid(c, h)) heap_down(c, 0);
        else {
            dchunk *d = h->c;
            c->heap[0] = c->heap[--c->nheap];
            if (c->nheap) heap_down(c, 0);
            release(c, d);
        }
        return 1;
    }
}

static int cmp_hit_key(const lmap_hit *p, const lmap_hit *q, int K) {
    int64_t kp = K ? p->b : p->a, kq = K ? q->b : q->a;
    if (kp != kq) return kp < kq ? -1 : 1;
    uint32_t mp = K ? p->a_member : p->b_member, mq = K ? q->a_member : q->b_member;
    if (mp != mq) return mp < mq ? -1 : 1;
    int64_t op = K ? p->a : p->b, oq = K ? q->a : q->b;
    if (op != oq) return op < oq ? -1 : 1;
    if (p->len != q->len) return p->len < q->len ? -1 : 1;
    return (int)p->strand - (int)q->strand;
}
static int cmp_hit_key0(const void *x, const void *y) { return cmp_hit_key(x, y, 0); }
static int cmp_hit_key1(const void *x, const void *y) { return cmp_hit_key(x, y, 1); }

int lmap_cursor_point(lmap_cursor *c, int64_t pos, lmap_hit *out, int cap) {
    if (!c || pos < 0 || cap < 0 || (cap > 0 && !out)) return -1;
    int count = 0;
    uint32_t p0 = gallop_u64(c->pmax, c->npos, pos, c->point_hint);
    c->point_hint = p0;
    for (uint32_t p = p0; p < c->npos; p++) {
        if (cur_kmin(c, p) > pos) break;
        c->st.chunks_examined++;
        if (cur_kend(c, p) <= pos) continue;
        dchunk *d = acquire(c, p, pos, pos + 1);
        if (!d) return -1;
        uint32_t i0 = gallop_i64(d->mep, d->n, pos, d->hint);
        d->hint = i0;
        for (uint32_t i = i0; i < d->n; i++) {
            const lmap_run *r = &d->runs[i];
            if (rkey(r, c->K) > pos) break;
            if (rkey(r, c->K) + r->len <= pos) continue;
            if (count < cap) {
                lmap_hit *h = &out[count];
                h->a = r->a; h->b = r->b; h->len = r->len; h->strand = r->strand;
                h->a_member = c->K ? d->other_member : d->key_member;
                h->b_member = c->K ? d->key_member : d->other_member;
            }
            count++;
        }
        release(c, d);
    }
    if (count > 1 && cap > 0) {
        qsort(out, (size_t)(count < cap ? count : cap), sizeof *out, c->K ? cmp_hit_key1 : cmp_hit_key0);
    }
    return count;
}

int lmap_cursor_extent(const lmap_cursor *c, int64_t *lo, int64_t *hi) {
    if (lo) *lo = 0;
    if (hi) *hi = 0;
    if (!c || c->npos == 0) return 0;
    if (lo) *lo = cur_kmin(c, 0);            /* positions are in key_min order */
    if (hi) *hi = (int64_t)c->pmax[c->npos - 1];
    return 1;
}

uint32_t lmap_cursor_n_chunks(const lmap_cursor *c) { return c ? c->npos : 0; }

void lmap_cursor_get_stats(const lmap_cursor *c, lmap_cursor_stats *st) {
    if (c && st) { *st = c->st; st->cached_bytes = c->cached; }
}

void lmap_cursor_close(lmap_cursor *c) {
    if (!c) return;
    heap_clear(c);
    if (c->htab)
        for (uint32_t i = 0; i < c->hcap; i++)
            for (dchunk *d = c->htab[i], *nx; d; d = nx) { nx = d->hnext; free(d->runs); free(d->mep); free(d); }
    free(c->htab); free(c->heap); free(c->own_ids); free(c->own_pmax);
    free(c);
}

int lmap_hit_clip(const lmap_hit *h, int axis, int64_t lo, int64_t hi, lmap_hit *out) {
    if (!h || !out || (axis != 0 && axis != 1)) return 0;
    lmap_run r = { h->a, h->b, h->len, h->strand };
    lmap_hit t;
    if (!clip_run(&r, axis, lo, hi, &t)) return 0;
    t.a_member = h->a_member; t.b_member = h->b_member;
    *out = t;
    return 1;
}

/* ============================================= APPLICATION SECTIONS, NAME RANGES */

int lmap_verify(lmap_file *f) {
    if (!f) return -1;
    enum { STEP = 1 << 20 };
    uint8_t *buf = malloc(STEP);
    if (!buf) return -1;
    uLong crc = crc32(0L, Z_NULL, 0);
    for (uint64_t done = 0; done < f->runs_len; ) {
        uint64_t k = f->runs_len - done < STEP ? f->runs_len - done : STEP;
        if (lmap_pread(f->io, buf, (int64_t)(f->runs_off + done), (int64_t)k) != 0) { free(buf); return -1; }
        crc = crc32(crc, buf, (uInt)k);
        done += k;
    }
    free(buf);
    if ((uint32_t)crc != f->runs_crc) return -1;
    for (uint32_t i = 0; i < f->n_chunks; i++) {
        const lmap_chunk *c = &f->chunks[i];
        uint8_t *raw = malloc(c->clen);
        if (!raw) return -1;
        int bad = lmap_pread(f->io, raw, (int64_t)(f->runs_off + c->off), (int64_t)c->clen) != 0 ||
                  crc_all(raw, c->clen) != c->crc;
        free(raw);
        if (bad) return -1;
    }
    return 0;
}

int lmap_read_section(lmap_file *f, const char *id, uint8_t **out, size_t *n) {
    if (!f || !id || !out || !n) return -1;
    for (uint32_t i = 0; i < f->n_xs; i++) {
        if (strcmp(f->xs[i].id, id) != 0) continue;
        uint8_t *p = NULL;
        if (read_section_body(f->io, f->xs[i].off, f->xs[i].len, &p) != 0) return -1;
        if (crc_all(p, f->xs[i].len) != f->xs[i].crc) { free(p); return -1; }
        *out = p; *n = (size_t)f->xs[i].len;
        return 0;
    }
    return 1;                                         /* absent */
}

int lmap_member_prefix(const lmap_file *f, int axis, const char *prefix,
                       uint32_t *first_rank, uint32_t *count) {
    if (!f || (axis != 0 && axis != 1) || !prefix || !first_rank || !count) return -1;
    size_t pl = strlen(prefix);
    uint32_t lo = 0, hi = f->n_mem[axis];
    while (lo < hi) {                                 /* first name >= prefix */
        uint32_t mid = lo + (hi - lo) / 2;
        if (strcmp(f->mem[axis][f->by_name[axis][mid]].name, prefix) < 0) lo = mid + 1; else hi = mid;
    }
    uint32_t k = lo;
    while (k < f->n_mem[axis] && strncmp(f->mem[axis][f->by_name[axis][k]].name, prefix, pl) == 0) k++;
    *first_rank = lo; *count = k - lo;
    return 0;
}

int32_t lmap_member_by_rank(const lmap_file *f, int axis, uint32_t rank) {
    if (!f || (axis != 0 && axis != 1) || rank >= f->n_mem[axis]) return -1;
    return (int32_t)f->by_name[axis][rank];
}
