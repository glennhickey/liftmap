/* liftmap -- command-line converters and lookups for liftmap files.
 *
 *   liftmap from-paf   in.paf[.gz]   out.lmap   [--swap] [--allow-overlap]
 *   liftmap from-chain in.chain[.gz] out.lmap   [--swap] [--allow-overlap]
 *   liftmap to-paf     in.lmap [out.paf]        [--max-gap N]
 *   liftmap to-chain   in.lmap [out.chain]      [--max-gap N]
 *   liftmap lift       in.lmap in.bed [out.bed] [--from a|b]
 *   liftmap dump       in.lmap [out.tsv]
 *   liftmap info       in.lmap
 *   liftmap verify     in.lmap
 *
 * Axis convention: axis a is the sequence named FIRST in each record -- the PAF query,
 * the chain target (tName).  So a -> b is each format's own lift direction: PAF query ->
 * target, UCSC liftOver's t -> q.  --swap exchanges the axes on import.
 *
 * Import keeps only what liftmap represents: which bases align to which, and on which
 * strand.  Record boundaries, scores, mapping qualities and match/mismatch counts are not
 * kept, so export regroups the runs into records (--max-gap) and writes placeholder
 * scores.  A round trip preserves every aligned base pair exactly.
 */
#include "../src/lmap_file.h"

#include <errno.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

static void die(const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    fprintf(stderr, "liftmap: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *p, size_t n) { p = realloc(p, n ? n : 1); if (!p) die("out of memory"); return p; }
static char *xstrdup(const char *s) { char *p = strdup(s); if (!p) die("out of memory"); return p; }

/* ------------------------------------------------------------------ line reading */

typedef struct { gzFile gz; char *buf; size_t cap; int64_t lineno; const char *path; } reader;

static void reader_open(reader *r, const char *path) {
    memset(r, 0, sizeof *r);
    r->path = path;
    r->gz = strcmp(path, "-") ? gzopen(path, "rb") : gzdopen(0, "rb");
    if (!r->gz) die("cannot open %s: %s", path, strerror(errno));
    r->cap = 1 << 16; r->buf = xmalloc(r->cap);
}

/* Next line without its newline, or NULL at EOF.  Lines of any length. */
static char *reader_line(reader *r) {
    size_t n = 0;
    for (;;) {
        if (!gzgets(r->gz, r->buf + n, (int)(r->cap - n))) {
            if (n == 0) return NULL;
            break;
        }
        n += strlen(r->buf + n);
        if (n > 0 && r->buf[n-1] == '\n') break;
        if (n + 1 < r->cap) continue;               /* short read, not a full buffer */
        r->cap *= 2; r->buf = xrealloc(r->buf, r->cap);
    }
    while (n > 0 && (r->buf[n-1] == '\n' || r->buf[n-1] == '\r')) r->buf[--n] = 0;
    r->lineno++;
    return r->buf;
}

static void reader_close(reader *r) { gzclose(r->gz); free(r->buf); }

/* ------------------------------------------------------------------ names -> ids */

typedef struct { char *name; int64_t length; } seqinfo;

typedef struct {
    seqinfo *seq; uint32_t n, cap;
    uint32_t *slot; uint32_t nslot;                 /* open addressing, id+1, 0 = empty */
} seqtab;

static uint64_t hash_str(const char *s) {
    uint64_t h = 1469598103934665603ull;
    for (; *s; s++) { h ^= (uint8_t)*s; h *= 1099511628211ull; }
    return h;
}

static void seqtab_grow(seqtab *t) {
    uint32_t ns = t->nslot ? t->nslot * 2 : 1024;
    uint32_t *sl = calloc(ns, sizeof *sl);
    if (!sl) die("out of memory");
    for (uint32_t i = 0; i < t->n; i++) {
        uint64_t h = hash_str(t->seq[i].name) & (ns - 1);
        while (sl[h]) h = (h + 1) & (ns - 1);
        sl[h] = i + 1;
    }
    free(t->slot); t->slot = sl; t->nslot = ns;
}

/* Id of `name`, adding it if new.  A sequence seen twice must have one length. */
static uint32_t seqtab_id(seqtab *t, const char *name, int64_t length, const reader *r) {
    if ((t->n + 1) * 2 > t->nslot) seqtab_grow(t);
    uint64_t h = hash_str(name) & (t->nslot - 1);
    while (t->slot[h]) {
        seqinfo *s = &t->seq[t->slot[h] - 1];
        if (!strcmp(s->name, name)) {
            if (s->length != length)
                die("%s:%" PRId64 ": %s has length %" PRId64 " here but %" PRId64 " earlier",
                    r->path, r->lineno, name, length, s->length);
            return t->slot[h] - 1;
        }
        h = (h + 1) & (t->nslot - 1);
    }
    if (t->n == UINT32_MAX - 1) die("too many sequences");
    if (t->n == t->cap) { t->cap = t->cap ? t->cap * 2 : 1024; t->seq = xrealloc(t->seq, t->cap * sizeof *t->seq); }
    t->seq[t->n].name = xstrdup(name); t->seq[t->n].length = length;
    t->slot[h] = t->n + 1;
    return t->n++;
}

/* ------------------------------------------------------------------ run buffer */

typedef struct { uint32_t am, bm; int64_t a, b, len; uint8_t strand; } run;

typedef struct {
    seqtab ax[2];
    run *r; size_t n, cap;
    int swap;
} importer;

/* One aligned block: [q, q+len) on the first-named sequence pairs with [t, t+len) on the
 * second, reversed if strand.  Both are forward-strand coordinates. */
static void add_block(importer *im, uint32_t first, uint32_t second,
                      int64_t p1, int64_t p2, int64_t len, int strand) {
    if (len <= 0) return;
    if (im->n == im->cap) { im->cap = im->cap ? im->cap * 2 : 1 << 16; im->r = xrealloc(im->r, im->cap * sizeof *im->r); }
    run *x = &im->r[im->n++];
    if (!im->swap) { x->am = first;  x->bm = second; x->a = p1; x->b = p2; }
    else           { x->am = second; x->bm = first;  x->a = p2; x->b = p1; }
    x->len = len; x->strand = (uint8_t)(strand != 0);
}

/* Intern the pair of sequence names on the right axes (first-named on axis a unless
 * --swap) and return their ids through *first / *second. */
static void intern_pair(importer *im, const reader *r, const char *n1, int64_t l1,
                        const char *n2, int64_t l2, uint32_t *first, uint32_t *second) {
    *first  = seqtab_id(&im->ax[im->swap ? 1 : 0], n1, l1, r);
    *second = seqtab_id(&im->ax[im->swap ? 0 : 1], n2, l2, r);
}

static int64_t parse_i64(const char *s, const reader *r, const char *what) {
    char *e; errno = 0;
    long long v = strtoll(s, &e, 10);
    if (errno || e == s || *e) die("%s:%" PRId64 ": bad %s '%s'", r->path, r->lineno, what, s);
    return (int64_t)v;
}

/* Split on tabs (PAF) or runs of whitespace (chain) into at most max fields. */
static int split(char *s, char **f, int max, int ws) {
    int n = 0;
    if (ws) {
        char *tok, *save = NULL;
        for (tok = strtok_r(s, " \t", &save); tok && n < max; tok = strtok_r(NULL, " \t", &save)) f[n++] = tok;
        return n;
    }
    f[n++] = s;
    for (char *p = s; *p && n < max; p++) if (*p == '\t') { *p = 0; f[n++] = p + 1; }
    return n;
}

static void import_paf(importer *im, const char *path) {
    reader r; reader_open(&r, path);
    char *line, *f[4096];
    int64_t nrec = 0;
    while ((line = reader_line(&r))) {
        if (!*line || *line == '#') continue;
        int nf = split(line, f, 4096, 0);
        if (nf < 12) die("%s:%" PRId64 ": PAF needs 12 columns, found %d", path, r.lineno, nf);
        int64_t qlen = parse_i64(f[1], &r, "query length"), qs = parse_i64(f[2], &r, "query start"),
                qe = parse_i64(f[3], &r, "query end");
        int64_t tlen = parse_i64(f[6], &r, "target length"), ts = parse_i64(f[7], &r, "target start"),
                te = parse_i64(f[8], &r, "target end");
        if (strcmp(f[4], "+") && strcmp(f[4], "-")) die("%s:%" PRId64 ": bad strand '%s'", path, r.lineno, f[4]);
        int rev = f[4][0] == '-';
        if (qs < 0 || qs > qe || qe > qlen || ts < 0 || ts > te || te > tlen)
            die("%s:%" PRId64 ": coordinates outside the sequence", path, r.lineno);
        const char *cg = NULL;
        for (int i = 12; i < nf; i++) if (!strncmp(f[i], "cg:Z:", 5)) { cg = f[i] + 5; break; }
        if (!cg) die("%s:%" PRId64 ": no cg:Z: CIGAR -- liftmap needs base-level alignments", path, r.lineno);
        uint32_t q, t;
        intern_pair(im, &r, f[0], qlen, f[5], tlen, &q, &t);
        /* The CIGAR walks the target forward; the query forward on '+', backward on '-'. */
        int64_t tp = ts, qp = rev ? qe : qs;
        for (const char *c = cg; *c; ) {
            char *e; errno = 0;
            long long L = strtoll(c, &e, 10);
            if (e == c || errno || L < 0 || !*e) die("%s:%" PRId64 ": bad CIGAR near '%.20s'", path, r.lineno, c);
            char op = *e; c = e + 1;
            if (L == 0) continue;                   /* some aligners emit 0-length ops */
            switch (op) {
            case 'M': case '=': case 'X':
                if (rev) { add_block(im, q, t, qp - L, tp, L, 1); qp -= L; }
                else     { add_block(im, q, t, qp, tp, L, 0);     qp += L; }
                tp += L; break;
            case 'I': qp += rev ? -L : L; break;
            case 'D': case 'N': tp += L; break;
            default: die("%s:%" PRId64 ": CIGAR op '%c' not supported", path, r.lineno, op);
            }
        }
        if (tp != te || qp != (rev ? qs : qe))
            die("%s:%" PRId64 ": CIGAR does not span the record's coordinates", path, r.lineno);
        nrec++;
    }
    reader_close(&r);
    fprintf(stderr, "liftmap: %" PRId64 " PAF records\n", nrec);
}

static void import_chain(importer *im, const char *path) {
    reader r; reader_open(&r, path);
    char *line, *f[16];
    int64_t nrec = 0;
    int in_chain = 0;
    uint32_t tid = 0, qid = 0;
    int64_t tsize = 0, qsize = 0, tp = 0, qp = 0, tend = 0, qend = 0;
    int trev = 0, qrev = 0;
    while ((line = reader_line(&r))) {
        if (*line == '#') continue;
        int nf = split(line, f, 16, 1);
        if (nf == 0) { if (in_chain) die("%s:%" PRId64 ": chain ends without a final block", path, r.lineno); continue; }
        if (!strcmp(f[0], "chain")) {
            if (in_chain) die("%s:%" PRId64 ": new chain before the last one ended", path, r.lineno);
            if (nf < 12) die("%s:%" PRId64 ": chain header needs 12 fields", path, r.lineno);
            tsize = parse_i64(f[3], &r, "tSize"); qsize = parse_i64(f[8], &r, "qSize");
            trev = f[4][0] == '-'; qrev = f[9][0] == '-';
            tp = parse_i64(f[5], &r, "tStart"); tend = parse_i64(f[6], &r, "tEnd");
            qp = parse_i64(f[10], &r, "qStart"); qend = parse_i64(f[11], &r, "qEnd");
            if (tp < 0 || tp > tend || tend > tsize || qp < 0 || qp > qend || qend > qsize)
                die("%s:%" PRId64 ": coordinates outside the sequence", path, r.lineno);
            intern_pair(im, &r, f[2], tsize, f[7], qsize, &tid, &qid);
            in_chain = 1; nrec++;
            continue;
        }
        if (!in_chain) die("%s:%" PRId64 ": alignment data outside a chain", path, r.lineno);
        if (nf != 1 && nf != 3) die("%s:%" PRId64 ": expected 'size dt dq' or a final 'size'", path, r.lineno);
        int64_t size = parse_i64(f[0], &r, "block size");
        /* Chain coordinates are on each side's own strand; convert to forward. */
        int64_t a = trev ? tsize - (tp + size) : tp;
        int64_t b = qrev ? qsize - (qp + size) : qp;
        add_block(im, tid, qid, a, b, size, trev != qrev);
        tp += size; qp += size;
        if (nf == 3) {
            int64_t dt = parse_i64(f[1], &r, "dt"), dq = parse_i64(f[2], &r, "dq");
            if (dt < 0 || dq < 0) die("%s:%" PRId64 ": negative gap", path, r.lineno);
            tp += dt; qp += dq;
        } else {
            if (tp != tend || qp != qend) die("%s:%" PRId64 ": blocks do not reach the chain's end", path, r.lineno);
            in_chain = 0;
        }
    }
    if (in_chain) die("%s: file ends inside a chain", path);
    reader_close(&r);
    fprintf(stderr, "liftmap: %" PRId64 " chains\n", nrec);
}

static int cmp_run(const void *x, const void *y) {
    const run *p = x, *q = y;
    if (p->am != q->am) return p->am < q->am ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->bm != q->bm) return p->bm < q->bm ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    if (p->len != q->len) return p->len < q->len ? -1 : 1;
    return (int)p->strand - (int)q->strand;
}

/* A run's diagonal: constant along a colinear stretch.  Forward pairs a+i with b+i, so
 * b - a; reverse pairs a+i with b+len-1-i, so a + b + len. */
static int64_t diag(const run *x) { return x->strand ? x->a + x->b + x->len : x->b - x->a; }

static int cmp_diag(const void *x, const void *y) {
    const run *p = x, *q = y;
    if (p->am != q->am) return p->am < q->am ? -1 : 1;
    if (p->bm != q->bm) return p->bm < q->bm ? -1 : 1;
    if (p->strand != q->strand) return (int)p->strand - (int)q->strand;
    int64_t dp = diag(p), dq = diag(q);
    if (dp != dq) return dp < dq ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return 0;
}

/* Reduce the runs to the aligned base-pair SET, as maximal runs: on each diagonal, union
 * runs that overlap or abut.  This is canonical -- the same base pairs always give the
 * same runs, however the input split them into records or blocks, and duplicate records
 * collapse. */
static void canonicalize(importer *im) {
    qsort(im->r, im->n, sizeof *im->r, cmp_diag);
    size_t m = 0;
    for (size_t i = 0; i < im->n; i++) {
        run *x = &im->r[i];
        if (m > 0) {
            run *p = &im->r[m-1];
            if (p->am == x->am && p->bm == x->bm && p->strand == x->strand &&
                diag(p) == diag(x) && x->a <= p->a + p->len) {
                int64_t d = diag(p);
                int64_t end = x->a + x->len > p->a + p->len ? x->a + x->len : p->a + p->len;
                p->len = end - p->a;
                p->b = p->strand ? d - p->a - p->len : d + p->a;
                continue;
            }
        }
        im->r[m++] = *x;
    }
    im->n = m;
    qsort(im->r, im->n, sizeof *im->r, cmp_run);
}

static void write_import(importer *im, const char *out, const char *profile,
                         const char *source, int allow_overlap) {
    canonicalize(im);
    int64_t overlaps = 0; size_t first_ov = 0;
    for (size_t i = 1; i < im->n; i++)
        if (im->r[i].am == im->r[i-1].am && im->r[i].a < im->r[i-1].a + im->r[i-1].len)
            if (overlaps++ == 0) first_ov = i;
    if (overlaps && !allow_overlap) {
        run *x = &im->r[first_ov];
        die("%" PRId64 " runs overlap an earlier one on axis a (first at %s:%" PRId64 "): axis a "
            "maps each base at most once unless --allow-overlap, which stores the file in "
            "order b.  --swap to index the other side as axis a.", overlaps,
            im->ax[0].seq[x->am].name, x->a);
    }

    lmap_writer *w = lmap_writer_open(out, profile, NULL);
    if (!w) die("cannot write %s: %s", out, strerror(errno));
    if (overlaps && (lmap_writer_set_order(w, LMAP_ORDER_B) || lmap_writer_set_a_overlap(w)))
        die("cannot set order b / overlap on %s", out);
    for (int ax = 0; ax < 2; ax++)
        for (uint32_t i = 0; i < im->ax[ax].n; i++)
            if (lmap_writer_add_member(w, ax, im->ax[ax].seq[i].name, (uint64_t)im->ax[ax].seq[i].length) != (int32_t)i)
                die("cannot add sequence %s", im->ax[ax].seq[i].name);
    for (size_t i = 0; i < im->n; i++) {
        run *x = &im->r[i];
        if (lmap_writer_add_run(w, x->am, x->bm, x->a, x->b, x->len, x->strand))
            die("writer refused run %s:%" PRId64 " -> %s:%" PRId64 " len %" PRId64,
                im->ax[0].seq[x->am].name, x->a, im->ax[1].seq[x->bm].name, x->b, x->len);
    }
    char meta[1024];
    snprintf(meta, sizeof meta, "source\t%s\naxis_a\t%s\naxis_b\t%s\n", source,
             !strcmp(profile, "paf") ? (im->swap ? "target" : "query") : (im->swap ? "q" : "t"),
             !strcmp(profile, "paf") ? (im->swap ? "query" : "target") : (im->swap ? "t" : "q"));
    if (lmap_writer_add_section(w, "x.liftmap.meta", meta, strlen(meta)) || lmap_writer_close(w))
        die("failed writing %s", out);
    fprintf(stderr, "liftmap: %zu runs, %u + %u sequences%s -> %s\n", im->n, im->ax[0].n, im->ax[1].n,
            overlaps ? " (axis a overlaps: order b)" : "", out);
}

/* ------------------------------------------------------------------ export */

typedef struct { lmap_hit *h; size_t n; } hitlist;

/* Every run of axis-a member m, in (a, b) order. */
static hitlist member_runs(lmap_file *f, uint32_t m) {
    hitlist hl = {0};
    const lmap_member *mem = lmap_member_at(f, 0, m);
    if (lmap_query_a(f, m, 0, (int64_t)mem->length, &hl.h, &hl.n, NULL))
        die("query failed on %s", mem->name);
    return hl;
}

/* Can hit y extend a record ending with hit x?  Same target, same strand, both axes
 * advancing, gaps within max_gap. */
static int extends(const lmap_hit *x, const lmap_hit *y, int64_t max_gap) {
    if (x->b_member != y->b_member || x->strand != y->strand) return 0;
    int64_t da = y->a - (x->a + x->len);
    int64_t db = x->strand ? x->b - (y->b + y->len) : y->b - (x->b + x->len);
    return da >= 0 && db >= 0 && da <= max_gap && db <= max_gap;
}

static FILE *open_out(const char *path) {
    if (!path || !strcmp(path, "-")) return stdout;
    FILE *o = fopen(path, "w");
    if (!o) die("cannot write %s: %s", path, strerror(errno));
    return o;
}

static void close_out(FILE *o, const char *path) {
    if (fflush(o) || (o != stdout && fclose(o))) die("error writing %s", path ? path : "stdout");
}

static void export_records(lmap_file *f, FILE *o, int chain, int64_t max_gap) {
    uint32_t na = lmap_n_members(f, 0);
    int64_t id = 0;
    for (uint32_t rank = 0; rank < na; rank++) {
        uint32_t m = (uint32_t)lmap_member_by_rank(f, 0, rank);
        const lmap_member *am = lmap_member_at(f, 0, m);
        hitlist hl = member_runs(f, m);
        size_t i = 0;
        while (i < hl.n) {
            size_t j = i + 1;
            while (j < hl.n && extends(&hl.h[j-1], &hl.h[j], max_gap)) j++;
            const lmap_hit *h = hl.h;
            const lmap_member *bm = lmap_member_at(f, 1, h[i].b_member);
            int s = h[i].strand;
            int64_t alen = 0;
            for (size_t k = i; k < j; k++) alen += h[k].len;
            int64_t a0 = h[i].a, a1 = h[j-1].a + h[j-1].len;
            int64_t b0 = s ? h[j-1].b : h[i].b, b1 = s ? h[i].b + h[i].len : h[j-1].b + h[j-1].len;
            if (chain) {
                /* t = axis a (forward), q = axis b on the record's strand */
                int64_t qs = s ? (int64_t)bm->length - b1 : b0, qe = s ? (int64_t)bm->length - b0 : b1;
                fprintf(o, "chain %" PRId64 " %s %" PRIu64 " + %" PRId64 " %" PRId64 " %s %" PRIu64 " %c %"
                        PRId64 " %" PRId64 " %" PRId64 "\n", alen, am->name, am->length, a0, a1,
                        bm->name, bm->length, s ? '-' : '+', qs, qe, ++id);
                for (size_t k = i; k < j; k++) {
                    if (k + 1 < j) {
                        int64_t dt = h[k+1].a - (h[k].a + h[k].len);
                        int64_t dq = s ? h[k].b - (h[k+1].b + h[k+1].len) : h[k+1].b - (h[k].b + h[k].len);
                        fprintf(o, "%" PRId64 "\t%" PRId64 "\t%" PRId64 "\n", h[k].len, dt, dq);
                    } else fprintf(o, "%" PRId64 "\n\n", h[k].len);
                }
            } else {
                /* PAF: query = axis a, target = axis b.  The CIGAR walks the target forward,
                 * which on '-' visits the runs last to first. */
                int64_t blen = 0;
                fprintf(o, "%s\t%" PRIu64 "\t%" PRId64 "\t%" PRId64 "\t%c\t%s\t%" PRIu64 "\t%" PRId64 "\t%" PRId64,
                        am->name, am->length, a0, a1, s ? '-' : '+', bm->name, bm->length, b0, b1);
                /* block length = M + I + D */
                blen = (a1 - a0) + (b1 - b0) - alen;
                fprintf(o, "\t%" PRId64 "\t%" PRId64 "\t255\tcg:Z:", alen, blen);
                for (size_t k = 0; k < j - i; k++) {
                    size_t c = s ? j - 1 - k : i + k, nx = s ? c - 1 : c + 1;
                    fprintf(o, "%" PRId64 "M", h[c].len);
                    if (k + 1 < j - i) {
                        int64_t dq = s ? h[c].a - (h[nx].a + h[nx].len) : h[nx].a - (h[c].a + h[c].len);
                        int64_t dt = h[nx].b - (h[c].b + h[c].len);
                        if (dq) fprintf(o, "%" PRId64 "I", dq);
                        if (dt) fprintf(o, "%" PRId64 "D", dt);
                    }
                }
                fputc('\n', o);
            }
            i = j;
        }
        free(hl.h);
    }
}

/* ------------------------------------------------------------------ lift */

static void lift_bed(lmap_file *f, const char *bed, FILE *o, int from_b) {
    reader r; reader_open(&r, bed);
    int src = from_b ? 1 : 0, dst = 1 - src;
    char *line, *fl[64];
    int64_t nin = 0, nunmapped = 0, nout = 0;
    while ((line = reader_line(&r))) {
        if (!*line || *line == '#' || !strncmp(line, "track", 5) || !strncmp(line, "browser", 7)) continue;
        int nf = split(line, fl, 64, 0);
        if (nf < 3) die("%s:%" PRId64 ": BED needs 3 columns", bed, r.lineno);
        int64_t s = parse_i64(fl[1], &r, "start"), e = parse_i64(fl[2], &r, "end");
        if (s < 0 || e < s) die("%s:%" PRId64 ": bad interval", bed, r.lineno);
        nin++;
        int32_t m = lmap_member_by_name(f, src, fl[0]);
        lmap_hit *h = NULL; size_t n = 0;
        if (m >= 0 && (from_b ? lmap_query_b(f, (uint32_t)m, s, e, &h, &n, NULL)
                              : lmap_query_a(f, (uint32_t)m, s, e, &h, &n, NULL)))
            die("query failed on %s:%" PRId64 "-%" PRId64, fl[0], s, e);
        if (n == 0) nunmapped++;
        for (size_t k = 0; k < n; k++) {
            const lmap_member *dm = lmap_member_at(f, dst, dst ? h[k].b_member : h[k].a_member);
            int64_t ds = dst ? h[k].b : h[k].a;
            fprintf(o, "%s\t%" PRId64 "\t%" PRId64, dm->name, ds, ds + h[k].len);
            for (int c = 3; c < nf; c++) {
                const char *v = fl[c];
                if (c == 5 && h[k].strand && (!strcmp(v, "+") || !strcmp(v, "-"))) v = v[0] == '+' ? "-" : "+";
                fprintf(o, "\t%s", v);
            }
            fputc('\n', o);
            nout++;
        }
        free(h);
    }
    reader_close(&r);
    fprintf(stderr, "liftmap: %" PRId64 " intervals in, %" PRId64 " unmapped, %" PRId64 " lifted pieces out\n",
            nin, nunmapped, nout);
}

/* ------------------------------------------------------------------ main */

static void usage(void) {
    fprintf(stderr,
        "usage:\n"
        "  liftmap from-paf   in.paf[.gz]   out.lmap   [--swap] [--allow-overlap]\n"
        "  liftmap from-chain in.chain[.gz] out.lmap   [--swap] [--allow-overlap]\n"
        "  liftmap to-paf     in.lmap [out.paf]        [--max-gap N]\n"
        "  liftmap to-chain   in.lmap [out.chain]      [--max-gap N]\n"
        "  liftmap lift       in.lmap in.bed [out.bed] [--from a|b]\n"
        "  liftmap dump       in.lmap [out.tsv]\n"
        "  liftmap info       in.lmap\n"
        "  liftmap verify     in.lmap\n"
        "\n"
        "Axis a is the sequence named first in each record: the PAF query, the chain\n"
        "target.  a -> b is each format's lift direction; --swap exchanges them on import.\n"
        "Import keeps aligned base pairs and strands only; export regroups runs into\n"
        "records joined across gaps of at most --max-gap bp (default 10000).\n");
    exit(2);
}

static lmap_file *open_or_die(const char *path) {
    lmap_file *f = lmap_open(path);
    if (!f) die("%s is not a readable liftmap file", path);
    return f;
}

int main(int argc, char **argv) {
    if (argc < 3) usage();
    const char *cmd = argv[1];
    const char *pos[3] = {0}; int npos = 0;
    int swap = 0, allow_overlap = 0, from_b = 0;
    int64_t max_gap = 10000;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--swap")) swap = 1;
        else if (!strcmp(argv[i], "--allow-overlap")) allow_overlap = 1;
        else if (!strcmp(argv[i], "--max-gap") && i + 1 < argc) max_gap = strtoll(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--from") && i + 1 < argc) {
            const char *v = argv[++i];
            if (!strcmp(v, "a")) from_b = 0; else if (!strcmp(v, "b")) from_b = 1; else usage();
        }
        else if (argv[i][0] == '-' && argv[i][1] == '-') usage();
        else if (npos < 3) pos[npos++] = argv[i];
        else usage();
    }
    if (max_gap < 0) usage();

    if (!strcmp(cmd, "from-paf") || !strcmp(cmd, "from-chain")) {
        if (npos != 2) usage();
        importer im; memset(&im, 0, sizeof im); im.swap = swap;
        int paf = !strcmp(cmd, "from-paf");
        if (paf) import_paf(&im, pos[0]); else import_chain(&im, pos[0]);
        write_import(&im, pos[1], paf ? "paf" : "chain", pos[0], allow_overlap);
        return 0;
    }
    if (!strcmp(cmd, "to-paf") || !strcmp(cmd, "to-chain")) {
        if (npos < 1 || npos > 2) usage();
        lmap_file *f = open_or_die(pos[0]);
        FILE *o = open_out(pos[1]);
        export_records(f, o, !strcmp(cmd, "to-chain"), max_gap);
        close_out(o, pos[1]); lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "lift")) {
        if (npos < 2 || npos > 3) usage();
        lmap_file *f = open_or_die(pos[0]);
        FILE *o = open_out(pos[2]);
        lift_bed(f, pos[1], o, from_b);
        close_out(o, pos[2]); lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "dump")) {
        if (npos < 1 || npos > 2) usage();
        lmap_file *f = open_or_die(pos[0]);
        FILE *o = open_out(pos[1]);
        for (uint32_t rank = 0; rank < lmap_n_members(f, 0); rank++) {
            uint32_t m = (uint32_t)lmap_member_by_rank(f, 0, rank);
            hitlist hl = member_runs(f, m);
            for (size_t k = 0; k < hl.n; k++)
                fprintf(o, "%s\t%" PRId64 "\t%s\t%" PRId64 "\t%" PRId64 "\t%c\n",
                        lmap_member_at(f, 0, m)->name, hl.h[k].a,
                        lmap_member_at(f, 1, hl.h[k].b_member)->name, hl.h[k].b, hl.h[k].len,
                        hl.h[k].strand ? '-' : '+');
            free(hl.h);
        }
        close_out(o, pos[1]); lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "info")) {
        if (npos != 1) usage();
        lmap_file *f = open_or_die(pos[0]);
        uint64_t runs = 0, bytes = 0;
        for (uint32_t i = 0; i < lmap_n_chunks(f); i++) { runs += lmap_chunk_at(f, i)->n_runs; bytes += lmap_chunk_at(f, i)->clen; }
        printf("profile\t%s\norder\t%c\naxis_a_overlaps\t%s\nsequences_a\t%u\nsequences_b\t%u\n"
               "chunks\t%u\nruns\t%" PRIu64 "\npayload_bytes\t%" PRIu64 "\nbytes_per_run\t%.3f\n",
               lmap_profile(f), lmap_order(f) == LMAP_ORDER_B ? 'b' : 'a', lmap_a_overlap(f) ? "yes" : "no",
               lmap_n_members(f, 0), lmap_n_members(f, 1), lmap_n_chunks(f), runs, bytes,
               runs ? (double)bytes / (double)runs : 0.0);
        uint8_t *meta; size_t n;
        if (lmap_read_section(f, "x.liftmap.meta", &meta, &n) == 0) { fwrite(meta, 1, n, stdout); free(meta); }
        lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "verify")) {
        if (npos != 1) usage();
        lmap_file *f = open_or_die(pos[0]);
        int bad = lmap_verify(f) != 0;
        for (uint32_t i = 0; !bad && i < lmap_n_chunks(f); i++) {
            lmap_run *r = xmalloc(lmap_chunk_at(f, i)->n_runs * sizeof *r);
            bad = lmap_read_chunk(f, i, r) != 0;
            free(r);
        }
        printf("%s\n", bad ? "CORRUPT" : "OK");
        lmap_close(f);
        return bad;
    }
    usage();
    return 2;
}
