/* liftmap -- command-line converters and lookups for liftmap files.
 *
 *   liftmap from-paf   in.paf[.gz]   out.lmap   [--swap] [--allow-overlap] [--group-sep C [--group-fields N]]
 *   liftmap from-chain in.chain[.gz] out.lmap   [--swap] [--allow-overlap] [--group-sep C [--group-fields N]]
 *   liftmap to-paf     in.lmap [out.paf]        [--max-gap N]
 *   liftmap to-chain   in.lmap [out.chain]      [--max-gap N]
 *   liftmap lift       in.lmap in.bed [out.bed] [--from a|b] [--max-gap N] [--min-match F]
 *   liftmap coarsen    in.lmap out.lmap --max-gap N [--key a|b] [--mem BYTES] [--tmp-dir DIR]
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
#include "../src/liftmap.h"

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
    lmap_builder *b; int axis;                      /* declare new names here, if set */
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
            if (s->length != length && r)
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
    if (t->b && lmap_builder_add_member(t->b, t->axis, name, (uint64_t)length) != (int32_t)t->n)
        die("cannot declare %s: %s", name, lmap_builder_error(t->b));
    return t->n++;
}

/* ------------------------------------------------------------------ import */

/* Sequences are interned per axis as records arrive and declared to the builder, which
 * takes the runs in any order and reduces them to canonical maximal runs at close. */
typedef struct {
    seqtab ax[2];
    lmap_builder *b;
    int swap;
} importer;

/* One aligned block: [p1, p1+len) on the first-named sequence pairs with [p2, p2+len) on
 * the second, reversed if strand.  Both are forward-strand coordinates. */
static void add_block(importer *im, uint32_t first, uint32_t second,
                      int64_t p1, int64_t p2, int64_t len, int strand) {
    if (len <= 0) return;
    int rc = im->swap ? lmap_builder_add_run(im->b, second, first, p2, p1, len, (uint8_t)(strand != 0))
                      : lmap_builder_add_run(im->b, first, second, p1, p2, len, (uint8_t)(strand != 0));
    if (rc) die("%s", lmap_builder_error(im->b));
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

/* Group each axis's sequences by the name up to the n-th `sep`: "hg38.chr1" -> "hg38"
 * ('.', 1), "GCA_000001635.9.chr1" -> "GCA_000001635.9" ('.', 2), "HG002#1#chr1" ->
 * "HG002#1" ('#', 2).  A name with fewer separators, or nothing after them, joins no
 * group.  No single rule fits every naming scheme, so the caller says which applies. */
static const char *nth_sep(const char *name, char sep, int n) {
    const char *p = name;
    for (int k = 0; k < n; k++) {
        p = strchr(k ? p + 1 : p, sep);
        if (!p) return NULL;
    }
    return (p == name || !p[1]) ? NULL : p;
}

static void add_groups(lmap_builder *w, importer *im, char sep, int fields) {
    for (int ax = 0; ax < 2; ax++) {
        seqtab g = {0};
        for (uint32_t i = 0; i < im->ax[ax].n; i++) {
            const char *name = im->ax[ax].seq[i].name, *cut = nth_sep(name, sep, fields);
            if (!cut) continue;
            char *prefix = xmalloc((size_t)(cut - name) + 1);
            memcpy(prefix, name, (size_t)(cut - name)); prefix[cut - name] = 0;
            uint32_t n0 = g.n;
            uint32_t gid = seqtab_id(&g, prefix, 0, NULL);
            if (g.n > n0 && lmap_builder_add_group(w, ax, prefix) != (int32_t)gid) die("%s", lmap_builder_error(w));
            if (lmap_builder_set_member_group(w, ax, i, gid)) die("%s", lmap_builder_error(w));
            free(prefix);
        }
        for (uint32_t i = 0; i < g.n; i++) free(g.seq[i].name);
        free(g.seq); free(g.slot);
    }
}

static void finish_import(importer *im, const char *out, const char *profile,
                          const char *source, char group_sep, int group_fields) {
    lmap_builder *w = im->b;
    if (group_sep) add_groups(w, im, group_sep, group_fields);
    int paf = !strcmp(profile, "paf");
    if (lmap_builder_set_meta(w, "source", source) ||
        lmap_builder_set_meta(w, "axis.a", paf ? (im->swap ? "target" : "query") : (im->swap ? "q" : "t")) ||
        lmap_builder_set_meta(w, "axis.b", paf ? (im->swap ? "query" : "target") : (im->swap ? "t" : "q")))
        die("%s", lmap_builder_error(w));
    lmap_build_stats st;
    if (lmap_builder_close(im->b, &st))
        die("%s%s", st.error, st.overlapping_runs ? " (--allow-overlap to accept it, stored in order b; "
                                                    "--swap to index the other side as axis a)" : "");
    fprintf(stderr, "liftmap: %" PRIu64 " blocks -> %" PRIu64 " runs, %u + %u sequences, order %c%s%s -> %s\n",
            st.input_runs, st.runs, im->ax[0].n, im->ax[1].n, st.order == LMAP_ORDER_B ? 'b' : 'a',
            st.overlapping_runs ? " (axis a overlaps)" : "",
            st.spill_segments ? " (spilled)" : "", out);
    for (int ax = 0; ax < 2; ax++) {
        for (uint32_t i = 0; i < im->ax[ax].n; i++) free(im->ax[ax].seq[i].name);
        free(im->ax[ax].seq); free(im->ax[ax].slot);
    }
}

/* ------------------------------------------------------------------ export */

/* Stream every run of axis-a member m in a order: a cursor that keeps no chunks, so
 * memory is the chunks overlapping the current position, not the member. */
static lmap_cursor *member_cursor(lmap_file *f, uint32_t m) {
    lmap_cursor *c = lmap_cursor_open(f, 0, m, NULL, 0, 0);
    if (!c || lmap_cursor_seek(c, 0, (int64_t)lmap_member_at(f, 0, m)->length))
        die("cannot read %s", lmap_member_at(f, 0, m)->name);
    return c;
}

static int next_run(lmap_cursor *c, lmap_hit *h) {
    int rc = lmap_cursor_next(c, h);
    if (rc < 0) die("corrupt chunk");
    return rc;
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

/* One record: runs that extend each other, written as a chain or a PAF line. */
static void write_record(lmap_file *f, FILE *o, int chain, const lmap_hit *h, size_t n, int64_t id) {
    const lmap_member *am = lmap_member_at(f, 0, h[0].a_member);
    const lmap_member *bm = lmap_member_at(f, 1, h[0].b_member);
    int s = h[0].strand;
    int64_t alen = 0;
    for (size_t k = 0; k < n; k++) alen += h[k].len;
    int64_t a0 = h[0].a, a1 = h[n-1].a + h[n-1].len;
    int64_t b0 = s ? h[n-1].b : h[0].b, b1 = s ? h[0].b + h[0].len : h[n-1].b + h[n-1].len;
    if (chain) {
        /* t = axis a (forward), q = axis b on the record's strand */
        int64_t qs = s ? (int64_t)bm->length - b1 : b0, qe = s ? (int64_t)bm->length - b0 : b1;
        fprintf(o, "chain %" PRId64 " %s %" PRIu64 " + %" PRId64 " %" PRId64 " %s %" PRIu64 " %c %"
                PRId64 " %" PRId64 " %" PRId64 "\n", alen, am->name, am->length, a0, a1,
                bm->name, bm->length, s ? '-' : '+', qs, qe, id);
        for (size_t k = 0; k < n; k++) {
            if (k + 1 < n) {
                int64_t dt = h[k+1].a - (h[k].a + h[k].len);
                int64_t dq = s ? h[k].b - (h[k+1].b + h[k+1].len) : h[k+1].b - (h[k].b + h[k].len);
                fprintf(o, "%" PRId64 "\t%" PRId64 "\t%" PRId64 "\n", h[k].len, dt, dq);
            } else fprintf(o, "%" PRId64 "\n\n", h[k].len);
        }
        return;
    }
    /* PAF: query = axis a, target = axis b.  The CIGAR walks the target forward, which on
     * '-' visits the runs last to first.  Block length = M + I + D. */
    fprintf(o, "%s\t%" PRIu64 "\t%" PRId64 "\t%" PRId64 "\t%c\t%s\t%" PRIu64 "\t%" PRId64 "\t%" PRId64
            "\t%" PRId64 "\t%" PRId64 "\t255\tcg:Z:", am->name, am->length, a0, a1, s ? '-' : '+',
            bm->name, bm->length, b0, b1, alen, (a1 - a0) + (b1 - b0) - alen);
    for (size_t k = 0; k < n; k++) {
        size_t c = s ? n - 1 - k : k, nx = s ? c - 1 : c + 1;
        fprintf(o, "%" PRId64 "M", h[c].len);
        if (k + 1 < n) {
            int64_t dq = s ? h[c].a - (h[nx].a + h[nx].len) : h[nx].a - (h[c].a + h[c].len);
            int64_t dt = h[nx].b - (h[c].b + h[c].len);
            if (dq) fprintf(o, "%" PRId64 "I", dq);
            if (dt) fprintf(o, "%" PRId64 "D", dt);
        }
    }
    fputc('\n', o);
}

/* Regroup each member's runs into records, streaming: only the open record is held. */
static void export_records(lmap_file *f, FILE *o, int chain, int64_t max_gap) {
    int64_t id = 0;
    lmap_hit *rec = NULL; size_t n = 0, cap = 0;
    for (uint32_t rank = 0; rank < lmap_n_members(f, 0); rank++) {
        lmap_cursor *c = member_cursor(f, (uint32_t)lmap_member_by_rank(f, 0, rank));
        lmap_hit h;
        while (next_run(c, &h)) {
            if (n > 0 && !extends(&rec[n-1], &h, max_gap)) { write_record(f, o, chain, rec, n, ++id); n = 0; }
            if (n == cap) { cap = cap ? cap * 2 : 256; rec = xrealloc(rec, cap * sizeof *rec); }
            rec[n++] = h;
        }
        if (n > 0) { write_record(f, o, chain, rec, n, ++id); n = 0; }
        lmap_cursor_close(c);
    }
    free(rec);
}

/* ------------------------------------------------------------------ lift */

static void lift_bed(lmap_file *f, const char *bed, FILE *o, int from_b, int64_t max_gap,
                     double min_match) {
    reader r; reader_open(&r, bed);
    int src = from_b ? 1 : 0, dst = 1 - src;
    char *line, *fl[64];
    int64_t nin = 0, nunmapped = 0, nout = 0, nlow = 0;
    while ((line = reader_line(&r))) {
        if (!*line || *line == '#' || !strncmp(line, "track", 5) || !strncmp(line, "browser", 7)) continue;
        int nf = split(line, fl, 64, 0);
        if (nf < 3) die("%s:%" PRId64 ": BED needs 3 columns", bed, r.lineno);
        int64_t s = parse_i64(fl[1], &r, "start"), e = parse_i64(fl[2], &r, "end");
        if (s < 0 || e < s) die("%s:%" PRId64 ": bad interval", bed, r.lineno);
        nin++;
        int32_t m = lmap_member_by_name(f, src, fl[0]);
        lmap_lifted *h = NULL; size_t n = 0;
        if (m >= 0 && lmap_lift(f, src, (uint32_t)m, s, e, max_gap, &h, &n))
            die("lift failed on %s:%" PRId64 "-%" PRId64, fl[0], s, e);
        int64_t kept = 0;
        for (size_t k = 0; k < n; k++) {
            if (min_match > 0 && e > s && (double)h[k].aligned_bp < min_match * (double)(e - s)) { nlow++; continue; }
            fprintf(o, "%s\t%" PRId64 "\t%" PRId64, lmap_member_at(f, dst, h[k].member)->name, h[k].start, h[k].end);
            for (int c = 3; c < nf; c++) {
                const char *v = fl[c];
                if (c == 5 && h[k].strand && (!strcmp(v, "+") || !strcmp(v, "-"))) v = v[0] == '+' ? "-" : "+";
                fprintf(o, "\t%s", v);
            }
            fputc('\n', o);
            nout++; kept++;
        }
        if (kept == 0) nunmapped++;
        free(h);
    }
    reader_close(&r);
    fprintf(stderr, "liftmap: %" PRId64 " intervals in, %" PRId64 " unmapped, %" PRId64 " lifted "
            "intervals out%s", nin, nunmapped, nout, nlow ? "" : "\n");
    if (nlow) fprintf(stderr, ", %" PRId64 " below --min-match\n", nlow);
}

/* ------------------------------------------------------------------ main */

static void usage(void) {
    fprintf(stderr,
        "usage:\n"
        "  liftmap from-paf   in.paf[.gz]   out.lmap   [--swap] [--allow-overlap]\n"
        "  liftmap from-chain in.chain[.gz] out.lmap   [--swap] [--allow-overlap]\n"
        "                     [--group-sep C [--group-fields N]] [--mem BYTES] [--tmp-dir DIR]\n"
        "  liftmap to-paf     in.lmap [out.paf]        [--max-gap N]\n"
        "  liftmap to-chain   in.lmap [out.chain]      [--max-gap N]\n"
        "  liftmap lift       in.lmap in.bed [out.bed] [--from a|b] [--max-gap N] [--min-match F]\n"
        "  liftmap coarsen    in.lmap out.lmap --max-gap N [--key a|b]\n"
        "  liftmap dump       in.lmap [out.tsv]\n"
        "  liftmap info       in.lmap\n"
        "  liftmap verify     in.lmap\n"
        "\n"
        "Axis a is the sequence named first in each record: the PAF query, the chain\n"
        "target.  a -> b is each format's lift direction; --swap exchanges them on import.\n"
        "Import keeps aligned base pairs and strands only; export regroups runs into\n"
        "records joined across gaps of at most --max-gap bp (default 10000).\n"
        "Import holds --mem bytes of runs (default 1 GiB), then spills sorted runs to a\n"
        "temporary file in --tmp-dir (default: the output's directory).\n"
        "lift writes one row per aligned piece; --max-gap N merges pieces on the same\n"
        "target and strand across gaps of at most N bp on both axes, and --min-match F\n"
        "drops a merged row whose aligned bases are under F of the input interval.\n"
        "coarsen chains runs across gaps of at most --max-gap bp on both axes into longer,\n"
        "approximate runs whose length is their span on the --key axis (default b).\n"
        "--group-sep C groups sequences into genomes by the name up to the N-th C\n"
        "(--group-fields, default 1): '.' 1 for hg38.chr1, '.' 2 for GCA_000001635.9.chr1,\n"
        "'#' 2 for PanSN HG002#1#chr1.\n");
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
    char group_sep = 0;
    int group_fields = 1;
    size_t mem_bytes = 0;
    const char *tmp_dir = NULL;
    int64_t max_gap = 10000;
    int key_b = 1, gap_given = 0;
    double min_match = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--swap")) swap = 1;
        else if (!strcmp(argv[i], "--allow-overlap")) allow_overlap = 1;
        else if (!strcmp(argv[i], "--group-sep") && i + 1 < argc) {
            const char *v = argv[++i];
            if (strlen(v) != 1) usage();
            group_sep = v[0];
        }
        else if (!strcmp(argv[i], "--mem") && i + 1 < argc) mem_bytes = (size_t)strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--tmp-dir") && i + 1 < argc) tmp_dir = argv[++i];
        else if (!strcmp(argv[i], "--group-fields") && i + 1 < argc) {
            group_fields = atoi(argv[++i]);
            if (group_fields < 1) usage();
        }
        else if (!strcmp(argv[i], "--max-gap") && i + 1 < argc) { max_gap = strtoll(argv[++i], NULL, 10); gap_given = 1; }
        else if (!strcmp(argv[i], "--min-match") && i + 1 < argc) min_match = atof(argv[++i]);
        else if (!strcmp(argv[i], "--key") && i + 1 < argc) {
            const char *v = argv[++i];
            if (!strcmp(v, "a")) key_b = 0; else if (!strcmp(v, "b")) key_b = 1; else usage();
        }
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
        int paf = !strcmp(cmd, "from-paf");
        importer im; memset(&im, 0, sizeof im); im.swap = swap;
        im.b = lmap_builder_open(pos[1], paf ? "paf" : "chain", mem_bytes, tmp_dir);
        if (!im.b) die("cannot write %s", pos[1]);
        lmap_builder_allow_overlap(im.b, allow_overlap);
        im.ax[0].b = im.ax[1].b = im.b; im.ax[1].axis = 1;
        if (paf) import_paf(&im, pos[0]); else import_chain(&im, pos[0]);
        finish_import(&im, pos[1], paf ? "paf" : "chain", pos[0], group_sep, group_fields);
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
    if (!strcmp(cmd, "coarsen")) {
        if (npos != 2 || !gap_given) usage();
        lmap_file *f = open_or_die(pos[0]);
        lmap_builder *b = lmap_builder_open(pos[1], lmap_profile(f), mem_bytes, tmp_dir);
        if (!b) die("cannot write %s", pos[1]);
        lmap_builder_set_order(b, lmap_order(f));
        lmap_coarsen_stats cs;
        if (lmap_builder_copy_layout(b, f) || lmap_coarsen(f, b, key_b, max_gap, &cs))
            die("%s", lmap_builder_error(b));
        lmap_build_stats st;
        if (lmap_builder_close(b, &st)) die("%s", st.error);
        fprintf(stderr, "liftmap: %" PRIu64 " runs -> %" PRIu64 " chains -> %" PRIu64 " runs "
                "(%.1fx), %" PRIu64 " bp clipped at sequence ends, order %c -> %s\n",
                cs.runs_in, cs.chains_out, st.runs,
                st.runs ? (double)cs.runs_in / (double)st.runs : 0.0, cs.clipped_bp,
                st.order == LMAP_ORDER_B ? 'b' : 'a', pos[1]);
        lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "lift")) {
        if (npos < 2 || npos > 3) usage();
        lmap_file *f = open_or_die(pos[0]);
        FILE *o = open_out(pos[2]);
        lift_bed(f, pos[1], o, from_b, gap_given ? max_gap : -1, min_match);
        close_out(o, pos[2]); lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "dump")) {
        if (npos < 1 || npos > 2) usage();
        lmap_file *f = open_or_die(pos[0]);
        FILE *o = open_out(pos[1]);
        for (uint32_t rank = 0; rank < lmap_n_members(f, 0); rank++) {
            uint32_t m = (uint32_t)lmap_member_by_rank(f, 0, rank);
            lmap_cursor *c = member_cursor(f, m);
            lmap_hit h;
            while (next_run(c, &h))
                fprintf(o, "%s\t%" PRId64 "\t%s\t%" PRId64 "\t%" PRId64 "\t%c\n",
                        lmap_member_at(f, 0, m)->name, h.a,
                        lmap_member_at(f, 1, h.b_member)->name, h.b, h.len, h.strand ? '-' : '+');
            lmap_cursor_close(c);
        }
        close_out(o, pos[1]); lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "info")) {
        if (npos != 1) usage();
        lmap_file *f = open_or_die(pos[0]);
        lmap_file_stats fs;
        lmap_get_stats(f, &fs);
        uint64_t runs = fs.runs, bytes = fs.payload_bytes;
        printf("profile\t%s\norder\t%c\naxis_a_overlaps\t%s\nsequences_a\t%u\nsequences_b\t%u\n"
               "chunks\t%u\nruns\t%" PRIu64 "\npayload_bytes\t%" PRIu64 "\nbytes_per_run\t%.3f\n",
               lmap_profile(f), lmap_order(f) == LMAP_ORDER_B ? 'b' : 'a', lmap_a_overlap(f) ? "yes" : "no",
               lmap_n_members(f, 0), lmap_n_members(f, 1), fs.chunks, runs, bytes,
               runs ? (double)bytes / (double)runs : 0.0);
        for (uint32_t i = 0; i < lmap_meta_count(f); i++) {
            const char *k, *v;
            lmap_meta_at(f, i, &k, &v);
            printf("meta.%s\t%s\n", k, v);
        }
        for (int ax = 0; ax < 2; ax++)
            for (uint32_t g = 0; g < lmap_n_groups(f, ax); g++) {
                uint32_t nm; uint64_t len;
                lmap_group_members(f, ax, g, NULL, &nm, &len);
                printf("group.%c\t%s\t%u sequences\t%" PRIu64 " bp\n", ax ? 'b' : 'a', lmap_group_name(f, ax, g), nm, len);
            }
        lmap_close(f);
        return 0;
    }
    if (!strcmp(cmd, "verify")) {
        if (npos != 1) usage();
        lmap_file *f = open_or_die(pos[0]);
        int bad = lmap_verify(f) != 0;
        printf("%s\n", bad ? "CORRUPT" : "OK");
        lmap_close(f);
        return bad;
    }
    usage();
    return 2;
}
