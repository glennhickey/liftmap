/* imap_file -- see imap_file.h and doc/SPEC.md S2. */
#include "imap_file.h"

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

struct imap_writer {
    FILE    *fp;
    char    *path, *profile, *schema;
    uint32_t count;  uint64_t bspan;  int codec;  int order;

    wmember *mem[2];  uint32_t n_mem[2], cap_mem[2];

    imap_run *pend;  uint32_t n_pend, cap_pend;      /* runs in the open chunk */
    uint32_t  pend_amem, pend_bmem;
    uint64_t  pend_bmin, pend_bmax;

    imap_chunk *chunks; uint32_t n_chunks, cap_chunks;
    vec runs;                                        /* the runs section body */
    int failed;
};

imap_writer *imap_writer_open(const char *path, const char *profile, const char *schema) {
    imap_writer *w = calloc(1, sizeof *w);
    if (!w) return NULL;
    w->fp = fopen(path, "wb");
    if (!w->fp) { free(w); return NULL; }
    w->path = strdup(path);
    w->profile = strdup(profile ? profile : "");
    w->schema = strdup(schema ? schema : "");
    w->count = 8192; w->bspan = 1000000; w->codec = IMAP_CODEC_DEFLATE; w->order = IMAP_ORDER_A;
    if (!w->path || !w->profile || !w->schema) { imap_writer_abort(w); return NULL; }
    return w;
}

int imap_writer_set_params(imap_writer *w, uint32_t count, uint64_t bspan, int codec) {
    if (!w || w->n_chunks || w->n_pend) return -1;      /* only before any run */
    if (count == 0 || count > IMAP_MAX_CHUNK_RUNS) return -1;
    if (codec != IMAP_CODEC_NONE && codec != IMAP_CODEC_DEFLATE) return -1;
    w->count = count; w->bspan = bspan; w->codec = codec; return 0;
}

int imap_writer_set_order(imap_writer *w, int order) {
    if (!w || w->n_chunks || w->n_pend) return -1;      /* only before any run */
    if (order != IMAP_ORDER_A && order != IMAP_ORDER_B) return -1;
    w->order = order; return 0;
}

int32_t imap_writer_add_member(imap_writer *w, int axis, const char *name, uint64_t length) {
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
static int cmp_run_b(const void *x, const void *y) {
    const imap_run *p = (const imap_run *)x, *q = (const imap_run *)y;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    return (p->a > q->a) - (p->a < q->a);
}

static int flush_chunk(imap_writer *w) {
    if (w->n_pend == 0) return 0;
    /* The chunk was cut in axis-a order; under order b it is stored sorted by (b, a).
     * Axis a is a partial function, so (b, a) is a total order with no ties. */
    if (w->order == IMAP_ORDER_B) qsort(w->pend, w->n_pend, sizeof *w->pend, cmp_run_b);
    imap_buf st[IMAP_N_STREAMS];
    if (imap_encode_chunk(w->pend, w->n_pend, w->order, st) != 0) return -1;

    /* Chunk blob (SPEC 2.2): u32 raw_len[4], u32 stored_len[4], then the four stored
     * streams.  A stream is deflated only where that makes it smaller; otherwise it is
     * stored verbatim, and stored_len == raw_len says so -- there is no separate flag. */
    uint8_t *packed[IMAP_N_STREAMS] = {0};
    size_t   stored[IMAP_N_STREAMS];
    uint64_t rawsum = 0;
    int bad = 0;
    for (int i = 0; i < IMAP_N_STREAMS; i++) {
        rawsum += st[i].n;
        stored[i] = st[i].n;
        if (w->codec == IMAP_CODEC_DEFLATE && st[i].n > 0) {
            uint8_t *z = NULL; size_t zn = 0;
            if (raw_deflate(st[i].p, st[i].n, &z, &zn) != 0) { bad = 1; break; }
            if (zn < st[i].n) { packed[i] = z; stored[i] = zn; } else free(z);
        }
    }
    vec blob = {0};
    if (!bad && 16 + rawsum > IMAP_MAX_CHUNK_RAW) bad = 1;
    for (int i = 0; !bad && i < IMAP_N_STREAMS; i++) bad |= vec_put32(&blob, (uint32_t)st[i].n);
    for (int i = 0; !bad && i < IMAP_N_STREAMS; i++) bad |= vec_put32(&blob, (uint32_t)stored[i]);
    for (int i = 0; !bad && i < IMAP_N_STREAMS; i++)
        bad |= vec_put(&blob, packed[i] ? packed[i] : st[i].p, stored[i]);
    for (int i = 0; i < IMAP_N_STREAMS; i++) { free(packed[i]); imap_buf_free(&st[i]); }
    if (bad) { vec_free(&blob); return -1; }
    const uint8_t *body = blob.p; size_t bodylen = blob.n;
    uint32_t rawlen32 = (uint32_t)(16 + rawsum);

    if (w->n_chunks == w->cap_chunks) {
        uint32_t cap = w->cap_chunks ? w->cap_chunks*2 : 256;
        imap_chunk *c = realloc(w->chunks, cap * sizeof *c);
        if (!c) { vec_free(&blob); return -1; }
        w->chunks = c; w->cap_chunks = cap;
    }
    /* extents from min/max rather than first/last: under order b the runs are no
     * longer in axis-a order by the time we get here */
    int64_t amin = w->pend[0].a, amax = w->pend[0].a + w->pend[0].len;
    for (uint32_t k = 1; k < w->n_pend; k++) {
        if (w->pend[k].a < amin) amin = w->pend[k].a;
        if (w->pend[k].a + w->pend[k].len > amax) amax = w->pend[k].a + w->pend[k].len;
    }
    uint64_t a_span64 = (uint64_t)(amax - amin);
    uint64_t b_span64 = w->pend_bmax - w->pend_bmin;
    if (a_span64 > 0xFFFFFFFFull || b_span64 > 0xFFFFFFFFull) { vec_free(&blob); return -1; }

    imap_chunk *c = &w->chunks[w->n_chunks];
    memset(c, 0, sizeof *c);
    c->a_member = w->pend_amem; c->b_member = w->pend_bmem;
    c->a_min = (uint64_t)amin;
    c->a_span = (uint32_t)a_span64;
    c->b_min = w->pend_bmin;
    c->b_span = (uint32_t)b_span64;
    /* the first run in stored order is the min on the ordering axis, so only the other
     * axis's base needs storing */
    c->base = w->order == IMAP_ORDER_A
            ? (uint64_t)imap_b_enter(w->pend[0].b, w->pend[0].len, w->pend[0].strand)
            : (uint64_t)w->pend[0].a;
    c->off = w->runs.n;
    c->clen = (uint32_t)bodylen;
    c->rawlen = rawlen32;
    c->n_runs = w->n_pend;
    c->codec = (uint32_t)w->codec;

    int rc = vec_put(&w->runs, body, bodylen);
    vec_free(&blob);
    if (rc != 0) return -1;

    w->n_chunks++;
    w->n_pend = 0;
    return 0;
}

int imap_writer_add_run(imap_writer *w, uint32_t a_member, uint32_t b_member,
                        int64_t a, int64_t b, int64_t len, uint8_t strand) {
    if (!w || w->failed) return -1;
    if (a < 0 || b < 0 || len < 1) return -1;
    if (a_member >= w->n_mem[0] || b_member >= w->n_mem[1]) return -1;
    if (len > INT64_MAX - a || len > INT64_MAX - b) return -1;

    /* Axis a is a partial function (SPEC 1.2): a member's runs must not overlap or
     * go backwards, and that has to hold across chunks, not just inside one --
     * otherwise two chunks of the same member can claim the same bases. */
    wmember *am = &w->mem[0][a_member];
    if (am->has_runs && a < am->last_a_exit) { w->failed = 1; return -1; }

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
        imap_run *r = realloc(w->pend, cap * sizeof *r);
        if (!r) { w->failed = 1; return -1; }
        w->pend = r; w->cap_pend = cap;
    }
    imap_run *r = &w->pend[w->n_pend++];
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

static int build_memdir(imap_writer *w, int axis, vec *out) {
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
    const imap_chunk *p = (const imap_chunk *)x, *q = (const imap_chunk *)y;
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

int imap_writer_close(imap_writer *w) {
    if (!w) return -1;
    int rc = -1;
    vec foot = {0}, dira = {0}, dirb = {0}, ma = {0}, mb = {0};
    if (w->failed) goto done;
    if (flush_chunk(w) != 0) goto done;

    /* dir.a is sorted by (a_member, a_min).  Chunks are addressed by payload offset,
     * so reordering the directory moves nothing on disk.  Per-member order across
     * chunks was enforced at add_run, so within a member this is also a_min order
     * with disjoint extents -- which is what makes a binary search over it valid. */
    qsort(w->chunks, w->n_chunks, sizeof *w->chunks, cmp_chunk_a);
    for (int ax = 0; ax < 2; ax++)
        for (uint32_t i = 0; i < w->n_mem[ax]; i++) {
            w->mem[ax][i].seen = 0; w->mem[ax][i].first_chunk = 0; w->mem[ax][i].n_chunks = 0;
        }
    for (uint32_t i = 0; i < w->n_chunks; i++) {
        imap_chunk *c = &w->chunks[i];
        note_range(&w->mem[0][c->a_member], i, c->a_min, c->a_min + c->a_span);
        if (vec_put32(&dira, c->a_member) || vec_put32(&dira, c->b_member) ||
            vec_put64(&dira, c->a_min) || vec_put32(&dira, c->a_span) ||
            vec_put32(&dira, c->b_span) || vec_put64(&dira, c->b_min) ||
            vec_put64(&dira, c->base) || vec_put64(&dira, c->off) ||
            vec_put32(&dira, c->clen) || vec_put32(&dira, c->rawlen) ||
            vec_put32(&dira, c->n_runs) || vec_put32(&dira, c->codec)) goto done;
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
            imap_chunk *c = &w->chunks[perm[i].id];
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
        uint8_t hdr[IMAP_HEADER_SIZE];
        memset(hdr, 0, sizeof hdr);
        memcpy(hdr, IMAP_MAGIC, 8);
        st32(hdr + 8, IMAP_FORMAT_MAJOR); st32(hdr + 12, 0);
        st64(hdr + 16, w->order == IMAP_ORDER_B ? IMAP_FEAT_ORDER_B : 0);
        strncpy((char *)hdr + 24, w->profile, 31);
        hdr_crc = crc_all(hdr, sizeof hdr);
        if (write_all(w->fp, hdr, sizeof hdr) != 0) goto done;
    }

    /* --- sections --- */
    struct { const char *id; const uint8_t *p; size_t n; int codec; size_t raw; } sec[] = {
        { "schema", (const uint8_t *)w->schema, strlen(w->schema), IMAP_CODEC_NONE, strlen(w->schema) },
        { "mem.a",  ma.p,       ma.n,       IMAP_CODEC_NONE, ma.n       },
        { "mem.b",  mb.p,       mb.n,       IMAP_CODEC_NONE, mb.n       },
        { "runs",   w->runs.p,  w->runs.n,  IMAP_CODEC_NONE, w->runs.n  },
        { "dir.a",  dira.p,     dira.n,     IMAP_CODEC_NONE, dira.n     },
        { "dir.b",  dirb.p,     dirb.n,     IMAP_CODEC_NONE, dirb.n     },
    };
    const int NSEC = (int)(sizeof sec / sizeof sec[0]);
    uint64_t off[8], len[8], crc[8];
    uint64_t cur = IMAP_HEADER_SIZE;
    for (int i = 0; i < NSEC; i++) {
        uint8_t pad[8] = {0};
        size_t padn = (size_t)((8 - (cur % 8)) % 8);
        if (write_all(w->fp, pad, padn) != 0) goto done;
        cur += padn;
        off[i] = cur; len[i] = sec[i].n;
        crc[i] = crc_all(sec[i].p, (uint64_t)sec[i].n);
        if (write_all(w->fp, sec[i].p, sec[i].n) != 0) goto done;
        cur += sec[i].n;
    }

    /* --- footer --- */
    if (vec_put32(&foot, (uint32_t)NSEC) != 0) goto done;
    for (int i = 0; i < NSEC; i++) {
        uint8_t idlen = (uint8_t)strlen(sec[i].id);
        if (vec_put(&foot, &idlen, 1) || vec_put(&foot, sec[i].id, idlen) ||
            vec_put64(&foot, off[i]) || vec_put64(&foot, len[i]) ||
            vec_put64(&foot, (uint64_t)sec[i].raw)) goto done;
        uint8_t cb = (uint8_t)sec[i].codec, fl = 0;
        if (vec_put(&foot, &cb, 1) || vec_put(&foot, &fl, 1) ||
            vec_put32(&foot, (uint32_t)crc[i])) goto done;
    }
    {
        uint8_t pad[8] = {0};
        size_t padn = (size_t)((8 - (cur % 8)) % 8);
        if (write_all(w->fp, pad, padn) != 0) goto done;
        cur += padn;
        uint64_t foot_off = cur;
        if (write_all(w->fp, foot.p, foot.n) != 0) goto done;
        cur += foot.n;

        uint8_t tr[IMAP_TRAILER_SIZE];
        memset(tr, 0, sizeof tr);
        st64(tr, foot_off);
        st64(tr + 8, (uint64_t)foot.n);
        st32(tr + 16, crc_all(foot.p, (uint64_t)foot.n));
        st32(tr + 20, hdr_crc);       /* the header changes how every chunk decodes */
        memcpy(tr + 24, IMAP_MAGIC, 8);
        if (write_all(w->fp, tr, sizeof tr) != 0) goto done;
    }
    if (fflush(w->fp) != 0) goto done;
    rc = 0;
done:
    vec_free(&foot); vec_free(&dira); vec_free(&dirb); vec_free(&ma); vec_free(&mb);
    if (w->fp) { if (fclose(w->fp) != 0) rc = -1; w->fp = NULL; }
    if (rc != 0 && w->path) remove(w->path);   /* no half-written file left behind */
    imap_writer_abort(w);
    return rc;
}

void imap_writer_abort(imap_writer *w) {
    if (!w) return;
    if (w->fp) fclose(w->fp);
    for (int ax = 0; ax < 2; ax++) {
        for (uint32_t i = 0; i < w->n_mem[ax]; i++) free(w->mem[ax][i].name);
        free(w->mem[ax]);
    }
    free(w->pend); free(w->chunks);
    vec_free(&w->runs);
    free(w->path); free(w->profile); free(w->schema);
    free(w);
}

/* ================================================================== READER */

struct imap_file {
    imap_io *io; int own_io;
    char profile[32];
    char *schema;
    int order;
    imap_chunk *chunks; uint32_t n_chunks;
    uint32_t *dirb_id; uint64_t *dirb_max;
    imap_member *mem[2]; uint32_t n_mem[2];
    uint32_t *by_name[2];
    uint64_t runs_off, runs_len;
};

static int read_section_body(imap_io *io, uint64_t off, uint64_t len, uint8_t **out) {
    uint8_t *p = malloc(len ? (size_t)len : 1);
    if (!p) return -1;
    if (imap_pread(io, p, (int64_t)off, (int64_t)len) != 0) { free(p); return -1; }
    *out = p; return 0;
}

static int parse_memdir(imap_file *f, int axis, const uint8_t *p, size_t n) {
    if (n < 8) return -1;
    uint32_t cnt = ld32(p), blob = ld32(p + 4);
    /* 64-bit throughout: on ILP32 a size_t product here wraps, and cnt/blob are
     * both read straight off disk. */
    uint64_t need = 8ull + (uint64_t)cnt * IMAP_MEM_ENTRY + (uint64_t)blob + (uint64_t)cnt * 4ull;
    if (need > (uint64_t)n) return -1;
    const uint8_t *e = p + 8;
    const uint8_t *names = e + (size_t)cnt * IMAP_MEM_ENTRY;
    f->mem[axis] = calloc(cnt ? cnt : 1, sizeof *f->mem[axis]);
    f->by_name[axis] = calloc(cnt ? cnt : 1, sizeof **f->by_name);
    if (!f->mem[axis] || !f->by_name[axis]) return -1;
    f->n_mem[axis] = cnt;   /* set first: imap_close must free whatever we allocated */
    for (uint32_t i = 0; i < cnt; i++, e += IMAP_MEM_ENTRY) {
        uint32_t noff = ld32(e), nlen = ld32(e + 4);
        if ((uint64_t)noff + (uint64_t)nlen > (uint64_t)blob) return -1;
        imap_member *m = &f->mem[axis][i];
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

imap_file *imap_open_io(imap_io *io, int own_io) {
    imap_file *f = calloc(1, sizeof *f);
    if (!f) return NULL;
    f->io = io; f->own_io = own_io;

    int64_t sz = imap_io_size(io);
    uint8_t hdr[IMAP_HEADER_SIZE], tr[IMAP_TRAILER_SIZE];
    if (sz < IMAP_HEADER_SIZE + IMAP_TRAILER_SIZE) goto fail;
    if (imap_pread(io, hdr, 0, IMAP_HEADER_SIZE) != 0) goto fail;
    if (memcmp(hdr, IMAP_MAGIC, 8) != 0) goto fail;
    if (imap_pread(io, tr, sz - IMAP_TRAILER_SIZE, IMAP_TRAILER_SIZE) != 0) goto fail;
    if (memcmp(tr + 24, IMAP_MAGIC, 8) != 0) goto fail;
    /* Verify the header before trusting any field in it.  The feature bits decide how
     * every chunk decodes: a single flipped bit used to make some chunks decode, without
     * complaint, into different runs. */
    if (crc_all(hdr, IMAP_HEADER_SIZE) != ld32(tr + 20)) goto fail;
    if (ld32(hdr + 8) != IMAP_FORMAT_MAJOR) goto fail;       /* exact: no older reader */
    uint64_t feat = ld64(hdr + 16);
    if (feat & ~IMAP_FEAT_KNOWN) goto fail;                 /* a feature we don't know */
    f->order = (feat & IMAP_FEAT_ORDER_B) ? IMAP_ORDER_B : IMAP_ORDER_A;
    memcpy(f->profile, hdr + 24, 31); f->profile[31] = 0;
    uint64_t foff = ld64(tr), flen = ld64(tr + 8);
    if (foff < IMAP_HEADER_SIZE || flen > (uint64_t)sz || foff > (uint64_t)sz - flen) goto fail;

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
    uint8_t *dira = NULL, *dirb = NULL, *ma = NULL, *mb = NULL;
    uint64_t dira_n = 0, dirb_n = 0, ma_n = 0, mb_n = 0;
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
        if (scodec != IMAP_CODEC_NONE || sflags != 0 || sraw != sl) goto footfail;

        uint8_t **dst = NULL; uint64_t *dstn = NULL;
        if      (!strcmp(id, "dir.a")) { dst = &dira; dstn = &dira_n; }
        else if (!strcmp(id, "dir.b")) { dst = &dirb; dstn = &dirb_n; }
        else if (!strcmp(id, "mem.a")) { dst = &ma;   dstn = &ma_n;   }
        else if (!strcmp(id, "mem.b")) { dst = &mb;   dstn = &mb_n;   }
        else if (!strcmp(id, "runs")) {
            if (seen_runs) goto footfail;         /* duplicate ids are not legal */
            seen_runs = 1;
            uint8_t *body = NULL;                 /* verify the payload's checksum too */
            if (read_section_body(io, so, sl, &body) != 0) goto footfail;
            int bad = (crc_all(body, sl) != scrc);
            free(body);
            if (bad) goto footfail;
            f->runs_off = so; f->runs_len = sl;
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
        } else continue;                          /* unknown section: skip, not an error */

        if (*dst) goto footfail;                  /* duplicate id: would leak the first */
        if (read_section_body(io, so, sl, dst) != 0) goto footfail;
        *dstn = sl;
        if (crc_all(*dst, sl) != scrc) goto footfail;
    }
    free(foot); foot = NULL;

    if (!dira || dira_n % IMAP_DIRA_ENTRY) goto tablefail;
    f->n_chunks = (uint32_t)(dira_n / IMAP_DIRA_ENTRY);
    f->chunks = calloc(f->n_chunks ? f->n_chunks : 1, sizeof *f->chunks);
    if (!f->chunks) goto tablefail;
    for (uint32_t i = 0; i < f->n_chunks; i++) {
        const uint8_t *e = dira + (size_t)i * IMAP_DIRA_ENTRY;
        imap_chunk *c = &f->chunks[i];
        c->a_member = ld32(e);      c->b_member = ld32(e + 4);
        c->a_min    = ld64(e + 8);  c->a_span   = ld32(e + 16);
        c->b_span   = ld32(e + 20); c->b_min    = ld64(e + 24);
        c->base     = ld64(e + 32); c->off      = ld64(e + 40);
        c->clen     = ld32(e + 48); c->rawlen   = ld32(e + 52);
        c->n_runs   = ld32(e + 56); c->codec = ld32(e + 60);
        if (c->codec != IMAP_CODEC_NONE && c->codec != IMAP_CODEC_DEFLATE) goto tablefail;
        /* n_runs sizes the decode buffer, so it must be answerable to the bytes that
         * are actually there: every run costs at least one varint byte in each of the
         * a, b and len streams plus a strand bit, after the 16-byte stream header.
         * rawlen in turn must be reachable from clen. */
        if (c->n_runs == 0 || c->n_runs > IMAP_MAX_CHUNK_RUNS) goto tablefail;
        if (c->rawlen > IMAP_MAX_CHUNK_RAW) goto tablefail;
        if ((uint64_t)c->rawlen < 16ull + 3ull * c->n_runs + (c->n_runs + 7ull) / 8) goto tablefail;
        /* clen = 32 + stored bytes, rawlen = 16 + raw bytes, and stored <= raw per stream */
        if (c->clen < 32) goto tablefail;
        if (c->codec == IMAP_CODEC_NONE ? (uint64_t)c->clen != (uint64_t)c->rawlen + 16
            : ((uint64_t)c->clen > (uint64_t)c->rawlen + 16 ||
               (uint64_t)c->rawlen > (uint64_t)IMAP_ZLIB_MAX_RATIO * c->clen + 64)) goto tablefail;
        if (c->off > f->runs_len || c->clen > f->runs_len - c->off) goto tablefail;
        /* a_min, b_min and base become int64 decode bases. */
        if (c->a_min > (uint64_t)INT64_MAX || c->b_min > (uint64_t)INT64_MAX ||
            c->base > (uint64_t)INT64_MAX) goto tablefail;
    }
    if (dirb && dirb_n == (uint64_t)f->n_chunks * IMAP_DIRB_ENTRY) {
        f->dirb_id  = malloc((f->n_chunks ? f->n_chunks : 1) * sizeof *f->dirb_id);
        f->dirb_max = malloc((f->n_chunks ? f->n_chunks : 1) * sizeof *f->dirb_max);
        if (!f->dirb_id || !f->dirb_max) goto tablefail;
        for (uint32_t i = 0; i < f->n_chunks; i++) {
            const uint8_t *e = dirb + (size_t)i * IMAP_DIRB_ENTRY;
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
            imap_member *m = &f->mem[ax][i];
            if (m->n_chunks > f->n_chunks ||
                m->first_chunk > f->n_chunks - m->n_chunks) goto tablefail;
            covered += m->n_chunks;
        }
        if (f->n_chunks && covered != f->n_chunks) goto tablefail;
    }
    /* mem.a ranges index dir.a: one member each, a_min increasing, extents disjoint */
    for (uint32_t mi = 0; mi < f->n_mem[0]; mi++) {
        const imap_member *m = &f->mem[0][mi];
        for (uint32_t k = 0; k < m->n_chunks; k++) {
            const imap_chunk *c = &f->chunks[m->first_chunk + k];
            if (c->a_member != mi) goto tablefail;
            if (k) {
                const imap_chunk *pc = &f->chunks[m->first_chunk + k - 1];
                if (c->a_min < pc->a_min + pc->a_span) goto tablefail;
            }
            if (c->a_min + c->a_span > m->length) goto tablefail;
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
            const imap_member *m = &f->mem[1][mi];
            uint64_t run_max = 0;
            for (uint32_t k = 0; k < m->n_chunks; k++) {
                uint32_t pos = m->first_chunk + k;
                const imap_chunk *c = &f->chunks[f->dirb_id[pos]];
                if (c->b_member != mi) goto tablefail;
                uint64_t end = c->b_min + c->b_span;
                if (end > m->length) goto tablefail;
                if (k && c->b_min < f->chunks[f->dirb_id[pos - 1]].b_min) goto tablefail;
                if (k == 0 || end > run_max) run_max = end;
                if (f->dirb_max[pos] != run_max) goto tablefail;
            }
        }
    }
    free(dira); free(dirb); free(ma); free(mb);
    return f;

footfail:
    free(foot);
tablefail:
    free(dira); free(dirb); free(ma); free(mb);
fail:
    imap_close(f);      /* honours own_io, so a failed imap_open_io(io,1) frees io */
    return NULL;
}

imap_file *imap_open(const char *path) {
    imap_io *io = imap_io_open_file(path);
    if (!io) return NULL;
    return imap_open_io(io, 1);   /* takes ownership on success and on failure */
}

void imap_close(imap_file *f) {
    if (!f) return;
    for (int ax = 0; ax < 2; ax++) {
        for (uint32_t i = 0; i < f->n_mem[ax]; i++) free(f->mem[ax][i].name);
        free(f->mem[ax]); free(f->by_name[ax]);
    }
    free(f->chunks); free(f->dirb_id); free(f->dirb_max); free(f->schema);
    if (f->own_io) imap_io_close(f->io);
    free(f);
}

const char *imap_profile(const imap_file *f)     { return f ? f->profile : NULL; }
const char *imap_schema_text(const imap_file *f) { return f ? f->schema : NULL; }
int         imap_order(const imap_file *f)       { return f ? f->order : -1; }
uint32_t    imap_n_chunks(const imap_file *f)    { return f ? f->n_chunks : 0; }
uint32_t    imap_n_members(const imap_file *f, int axis) {
    return (f && (axis == 0 || axis == 1)) ? f->n_mem[axis] : 0;
}
const imap_member *imap_member_at(const imap_file *f, int axis, uint32_t i) {
    if (!f || (axis != 0 && axis != 1) || i >= f->n_mem[axis]) return NULL;
    return &f->mem[axis][i];
}
const imap_chunk *imap_chunk_at(const imap_file *f, uint32_t i) {
    return (f && i < f->n_chunks) ? &f->chunks[i] : NULL;
}

int32_t imap_member_by_name(const imap_file *f, int axis, const char *name) {
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
static int decode_checked(const imap_file *f, const imap_chunk *c,
                          const imap_buf st[IMAP_N_STREAMS], imap_run *out) {
    int64_t base_a = f->order == IMAP_ORDER_A ? (int64_t)c->a_min : (int64_t)c->base;
    int64_t base_b = f->order == IMAP_ORDER_A ? (int64_t)c->base  : (int64_t)c->b_min;
    if (imap_decode_chunk(st, c->n_runs, f->order, base_a, base_b, out) != 0) return -1;
    uint64_t a_hi = c->a_min + c->a_span, b_hi = c->b_min + c->b_span;
    for (uint32_t k = 0; k < c->n_runs; k++) {
        const imap_run *r = &out[k];
        if ((uint64_t)r->a < c->a_min || (uint64_t)(r->a + r->len) > a_hi ||
            (uint64_t)r->b < c->b_min || (uint64_t)(r->b + r->len) > b_hi) return -1;
    }
    return 0;
}

int imap_read_chunk(imap_file *f, uint32_t i, imap_run *out) {
    const imap_chunk *c = imap_chunk_at(f, i);
    if (!c || !out || c->clen < 32) return -1;
    uint8_t *raw = malloc(c->clen);
    if (!raw) return -1;
    if (imap_pread(f->io, raw, (int64_t)(f->runs_off + c->off), (int64_t)c->clen) != 0) {
        free(raw); return -1;
    }
    /* Blob (SPEC 2.2): u32 raw_len[4], u32 stored_len[4], the four stored streams.
     * Every length is off disk; each is checked against the directory entry, which was
     * itself bounded at open, before anything is allocated from it. */
    uint32_t rl[IMAP_N_STREAMS], sl[IMAP_N_STREAMS];
    uint64_t rsum = 0, ssum = 0;
    int ok = 1;
    for (int k = 0; k < IMAP_N_STREAMS; k++) {
        rl[k] = ld32(raw + 4*k);
        sl[k] = ld32(raw + 16 + 4*k);
        rsum += rl[k]; ssum += sl[k];
        if (sl[k] > rl[k]) ok = 0;                            /* never stored larger */
        if (c->codec == IMAP_CODEC_NONE && sl[k] != rl[k]) ok = 0;
    }
    if (!ok || rsum + 16 != c->rawlen || ssum + 32 != c->clen) { free(raw); return -1; }

    uint8_t *inf = malloc(rsum ? (size_t)rsum : 1);
    if (!inf) { free(raw); return -1; }
    imap_buf st[IMAP_N_STREAMS];
    const uint8_t *q = raw + 32;
    uint8_t *o = inf;
    for (int k = 0; ok && k < IMAP_N_STREAMS; k++) {
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
static int clip_run(const imap_run *r, int axis, int64_t lo, int64_t hi, imap_hit *h) {
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

typedef struct { imap_hit *p; size_t n, cap; } hitvec;
static int hit_push(hitvec *v, const imap_hit *h) {
    if (v->n == v->cap) {
        size_t cap = v->cap ? v->cap * 2 : 64;
        imap_hit *q = realloc(v->p, cap * sizeof *q);
        if (!q) return -1;
        v->p = q; v->cap = cap;
    }
    v->p[v->n++] = *h;
    return 0;
}

/* Decode chunk ci and append its runs' overlaps with [lo,hi) on `axis`. */
static int scan_chunk(imap_file *f, uint32_t ci, int axis, int64_t lo, int64_t hi,
                      imap_run **buf, uint32_t *bufn, hitvec *out, imap_query_stats *st) {
    const imap_chunk *c = &f->chunks[ci];
    if (c->n_runs > *bufn) {
        imap_run *nb = realloc(*buf, (size_t)c->n_runs * sizeof *nb);
        if (!nb) return -1;
        *buf = nb; *bufn = c->n_runs;
    }
    if (imap_read_chunk(f, ci, *buf) != 0) return -1;
    if (st) st->chunks_decoded++;
    for (uint32_t k = 0; k < c->n_runs; k++) {
        imap_hit h;
        if (!clip_run(&(*buf)[k], axis, lo, hi, &h)) continue;
        h.a_member = c->a_member; h.b_member = c->b_member;
        if (hit_push(out, &h) != 0) return -1;
    }
    return 0;
}

static int cmp_hit_a(const void *x, const void *y) {
    const imap_hit *p = (const imap_hit *)x, *q = (const imap_hit *)y;
    return (p->a > q->a) - (p->a < q->a);
}

int imap_query_a(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *st) {
    if (!f || !out || !n || m >= f->n_mem[0] || lo < 0 || hi < lo) return -1;
    *out = NULL; *n = 0;
    if (st) memset(st, 0, sizeof *st);
    const imap_member *mem = &f->mem[0][m];
    /* first chunk of the member whose extent ends after lo: extents are disjoint and
     * increasing within a member (validated at open), so this is a binary search */
    uint32_t lo_i = mem->first_chunk, hi_i = mem->first_chunk + mem->n_chunks;
    while (lo_i < hi_i) {
        uint32_t mid = lo_i + (hi_i - lo_i) / 2;
        const imap_chunk *c = &f->chunks[mid];
        if ((int64_t)(c->a_min + c->a_span) <= lo) lo_i = mid + 1; else hi_i = mid;
    }
    hitvec v = {0}; imap_run *buf = NULL; uint32_t bufn = 0; int rc = 0;
    for (uint32_t ci = lo_i; ci < mem->first_chunk + mem->n_chunks; ci++) {
        if ((int64_t)f->chunks[ci].a_min >= hi) break;
        if (st) st->chunks_examined++;
        if (scan_chunk(f, ci, 0, lo, hi, &buf, &bufn, &v, st) != 0) { rc = -1; break; }
    }
    free(buf);
    if (rc != 0) { free(v.p); return -1; }
    /* Chunks are cut in axis-a order, so results are a-ordered across chunks; under
     * order b they are b-ordered within one, so sort.  a is unique within a member. */
    if (f->order == IMAP_ORDER_B && v.n > 1) qsort(v.p, v.n, sizeof *v.p, cmp_hit_a);
    *out = v.p; *n = v.n;
    return 0;
}

static int cmp_hit_b(const void *x, const void *y) {
    const imap_hit *p = (const imap_hit *)x, *q = (const imap_hit *)y;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    if (p->a_member != q->a_member) return p->a_member < q->a_member ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return 0;
}

int imap_query_b(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *st) {
    if (!f || !out || !n || m >= f->n_mem[1] || lo < 0 || hi < lo) return -1;
    if (!f->dirb_id) return -1;                         /* no reverse index in this file */
    *out = NULL; *n = 0;
    if (st) memset(st, 0, sizeof *st);
    const imap_member *mem = &f->mem[1][m];
    uint32_t first = mem->first_chunk, end = mem->first_chunk + mem->n_chunks;
    /* last dir.b position in this member with b_min < hi */
    uint32_t lo_i = first, hi_i = end;
    while (lo_i < hi_i) {
        uint32_t mid = lo_i + (hi_i - lo_i) / 2;
        if ((int64_t)f->chunks[f->dirb_id[mid]].b_min < hi) lo_i = mid + 1; else hi_i = mid;
    }
    hitvec v = {0}; imap_run *buf = NULL; uint32_t bufn = 0; int rc = 0;
    /* Walk backwards. prefix_max[j] bounds the b_end of every chunk at or before j in
     * this member, so once it is <= lo nothing earlier can overlap: stop. */
    for (uint32_t j = lo_i; j > first; j--) {
        uint32_t pos = j - 1;
        if ((int64_t)f->dirb_max[pos] <= lo) break;
        uint32_t ci = f->dirb_id[pos];
        const imap_chunk *c = &f->chunks[ci];
        if (st) st->chunks_examined++;
        if ((int64_t)(c->b_min + c->b_span) <= lo) continue;
        if (scan_chunk(f, ci, 1, lo, hi, &buf, &bufn, &v, st) != 0) { rc = -1; break; }
    }
    free(buf);
    if (rc != 0) { free(v.p); return -1; }
    if (v.n > 1) qsort(v.p, v.n, sizeof *v.p, cmp_hit_b);   /* qsort(NULL, 0) is UB */
    *out = v.p; *n = v.n;
    return 0;
}
