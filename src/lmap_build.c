/* lmap_build -- a writer front end that takes runs in any order.  See lmap_file.h.
 *
 * Runs are buffered and, over the memory budget, spilled as sorted segments to one
 * temporary file (unlinked as soon as it is created, so nothing outlives the process).
 * At close they are merged in (a_member, b_member, strand, diagonal, a) order, which puts
 * every run of one colinear stretch next to its neighbours: overlapping or abutting runs
 * on a diagonal are unioned into maximal runs -- the canonical form of the aligned
 * base-pair set.  The merged stream is grouped by a_member, so each member's runs are
 * then sorted on axis a in memory (one member at a time), checked for axis-a overlap,
 * and handed to the writer.  The storage order must be fixed before the writer's first
 * run, and whether axis a overlaps is only known once every member has been seen, so
 * when runs spilled the canonical runs go to a second temporary file first.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "lmap_file.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct { uint32_t am, bm; int64_t a, b, len; uint8_t strand; } brec;

struct lmap_builder {
    lmap_writer *w;
    char *path, *tmp_dir;
    size_t budget;                       /* bytes of runs held before spilling */
    int order;                           /* LMAP_ORDER_A / _B, or -1: choose */
    int allow_overlap;
    uint64_t *len[2]; uint32_t n_mem[2], cap_mem[2];   /* member lengths, for checks */
    brec *buf; size_t n, cap;
    /* pre-sorted mode: runs stream to the writer, merged with their predecessor */
    int presorted, started;
    brec pend; int have_pend;
    uint8_t *done_a;                     /* a-members already passed, to catch a revisit */
    uint32_t cur_a; int have_cur_a; int64_t cur_a_end;
    FILE *spill; uint64_t *seg_off, *seg_n; uint32_t n_seg, cap_seg;
    lmap_build_stats st;
    char err[512];
    int failed;
};

static int fail(lmap_builder *b, const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    vsnprintf(b->err, sizeof b->err, fmt, ap);
    va_end(ap);
    b->failed = 1;
    return -1;
}

const char *lmap_builder_error(const lmap_builder *b) { return b ? b->err : "no builder"; }

lmap_builder *lmap_builder_open(const char *path, const char *profile,
                                size_t mem_bytes, const char *tmp_dir) {
    if (!path) return NULL;
    lmap_builder *b = calloc(1, sizeof *b);
    if (!b) return NULL;
    b->w = lmap_writer_open(path, profile);
    b->path = strdup(path);
    if (tmp_dir && *tmp_dir) b->tmp_dir = strdup(tmp_dir);
    else {                                              /* the output's own directory */
        const char *sl = strrchr(path, '/');
        b->tmp_dir = sl ? strndup(path, (size_t)(sl - path) + (sl == path)) : strdup(".");
    }
    if (!b->w || !b->path || !b->tmp_dir) { lmap_builder_abort(b); return NULL; }
    b->budget = mem_bytes ? mem_bytes : (size_t)1 << 30;
    b->order = -1;
    return b;
}

/* ---- pass-throughs: everything but members, runs and the order ---- */

int lmap_builder_set_meta(lmap_builder *b, const char *key, const char *value) {
    if (!b || b->failed) return -1;
    return lmap_writer_set_meta(b->w, key, value) ? fail(b, "cannot set metadata %s", key ? key : "") : 0;
}
int32_t lmap_builder_add_group(lmap_builder *b, int axis, const char *name) {
    if (!b || b->failed) return -1;
    int32_t g = lmap_writer_add_group(b->w, axis, name);
    if (g < 0) return fail(b, "cannot add group %s (duplicate name?)", name ? name : "");
    return g;
}
int lmap_builder_set_member_group(lmap_builder *b, int axis, uint32_t member, uint32_t group) {
    if (!b || b->failed) return -1;
    return lmap_writer_set_member_group(b->w, axis, member, group)
         ? fail(b, "cannot put member %u in group %u", member, group) : 0;
}
int lmap_builder_add_section(lmap_builder *b, const char *id, const void *data, size_t n) {
    if (!b || b->failed) return -1;
    return lmap_writer_add_section(b->w, id, data, n)
         ? fail(b, "cannot add section %s (ids start x. and are unique)", id ? id : "") : 0;
}
int lmap_builder_set_params(lmap_builder *b, uint32_t count, uint64_t bspan, int codec) {
    if (!b || b->failed) return -1;
    return lmap_writer_set_params(b->w, count, bspan, codec) ? fail(b, "bad chunk parameters") : 0;
}

int lmap_builder_set_presorted(lmap_builder *b, int on) {
    if (!b || b->failed || b->n || b->n_seg || b->started) return -1;
    b->presorted = on != 0; return 0;
}

int32_t lmap_builder_add_member(lmap_builder *b, int axis, const char *name, uint64_t length) {
    if (!b || b->failed || (axis != 0 && axis != 1)) return -1;
    int32_t id = lmap_writer_add_member(b->w, axis, name, length);
    if (id < 0) return fail(b, "cannot add member %s", name ? name : "(null)"), -1;
    if (b->n_mem[axis] == b->cap_mem[axis]) {
        uint32_t cap = b->cap_mem[axis] ? b->cap_mem[axis] * 2 : 64;
        uint64_t *nl = realloc(b->len[axis], cap * sizeof *nl);
        if (!nl) return fail(b, "out of memory"), -1;
        b->len[axis] = nl; b->cap_mem[axis] = cap;
    }
    b->len[axis][b->n_mem[axis]++] = length;
    return id;
}

int lmap_builder_set_order(lmap_builder *b, int order) {
    if (!b || (order != LMAP_ORDER_A && order != LMAP_ORDER_B && order != -1)) return -1;
    b->order = order; return 0;
}

int lmap_builder_allow_overlap(lmap_builder *b, int allow) {
    if (!b) return -1;
    b->allow_overlap = allow != 0; return 0;
}

/* diagonal: constant along a colinear stretch (forward a+i~b+i, reverse a+i~b+len-1-i) */
static inline int64_t diag(const brec *r) { return r->strand ? r->a + r->b + r->len : r->b - r->a; }

static int cmp_diag(const brec *p, const brec *q) {
    if (p->am != q->am) return p->am < q->am ? -1 : 1;
    if (p->bm != q->bm) return p->bm < q->bm ? -1 : 1;
    if (p->strand != q->strand) return p->strand < q->strand ? -1 : 1;
    int64_t dp = diag(p), dq = diag(q);
    if (dp != dq) return dp < dq ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return p->len < q->len ? -1 : p->len > q->len;
}
static int qcmp_diag(const void *x, const void *y) { return cmp_diag(x, y); }

static int qcmp_a(const void *x, const void *y) {
    const brec *p = x, *q = y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->bm != q->bm) return p->bm < q->bm ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    if (p->len != q->len) return p->len < q->len ? -1 : 1;
    return (int)p->strand - (int)q->strand;
}

static FILE *temp_file(lmap_builder *b) {
    size_t n = strlen(b->tmp_dir) + 32;
    char *t = malloc(n);
    if (!t) return NULL;
    snprintf(t, n, "%s/lmap.spill.XXXXXX", b->tmp_dir);
    int fd = mkstemp(t);
    if (fd < 0) { free(t); return NULL; }
    unlink(t);                          /* gone from the directory now; freed on close/exit */
    free(t);
    FILE *f = fdopen(fd, "w+b");
    if (!f) close(fd);
    return f;
}

/* Sort the buffer by diagonal and append it to the spill file as one segment. */
static int spill(lmap_builder *b) {
    if (b->n == 0) return 0;
    if (!b->spill && !(b->spill = temp_file(b)))
        return fail(b, "cannot create a spill file in %s: %s", b->tmp_dir, strerror(errno));
    qsort(b->buf, b->n, sizeof *b->buf, qcmp_diag);
    if (fseeko(b->spill, 0, SEEK_END) != 0) return fail(b, "spill seek failed");
    off_t off = ftello(b->spill);
    if (fwrite(b->buf, sizeof *b->buf, b->n, b->spill) != b->n)
        return fail(b, "spill write failed (disk full?): %s", strerror(errno));
    if (b->n_seg == b->cap_seg) {
        uint32_t cap = b->cap_seg ? b->cap_seg * 2 : 16;
        uint64_t *no = realloc(b->seg_off, cap * sizeof *no), *nn;
        if (!no) return fail(b, "out of memory");
        b->seg_off = no;
        if (!(nn = realloc(b->seg_n, cap * sizeof *nn))) return fail(b, "out of memory");
        b->seg_n = nn; b->cap_seg = cap;
    }
    b->seg_off[b->n_seg] = (uint64_t)off; b->seg_n[b->n_seg] = b->n; b->n_seg++;
    b->st.spill_segments++;
    b->n = 0;
    return 0;
}

static int choose_order(lmap_builder *b);
static int presorted_run(lmap_builder *b, const brec *r);

int lmap_builder_add_run(lmap_builder *b, uint32_t a_member, uint32_t b_member,
                         int64_t a, int64_t bpos, int64_t len, uint8_t strand) {
    if (!b || b->failed) return -1;
    if (a_member >= b->n_mem[0] || b_member >= b->n_mem[1])
        return fail(b, "run names member %u/%u, which is not declared", a_member, b_member);
    if (a < 0 || bpos < 0 || len < 1 || len > INT64_MAX - a || len > INT64_MAX - bpos ||
        (uint64_t)(a + len) > b->len[0][a_member] || (uint64_t)(bpos + len) > b->len[1][b_member])
        return fail(b, "run a=%lld b=%lld len=%lld lies outside its members (lengths %llu, %llu)",
                    (long long)a, (long long)bpos, (long long)len,
                    (unsigned long long)b->len[0][a_member], (unsigned long long)b->len[1][b_member]);
    if (b->presorted) {
        brec r = { a_member, b_member, a, bpos, len, (uint8_t)(strand ? 1 : 0) };
        b->st.input_runs++;
        return presorted_run(b, &r);
    }
    if (b->n == b->cap) {
        size_t cap = b->cap ? b->cap * 2 : 4096;
        size_t lim = b->budget / sizeof *b->buf;
        if (lim < 1024) lim = 1024;
        if (b->cap >= lim) { if (spill(b)) return -1; }
        else {
            if (cap > lim) cap = lim;
            brec *nb = realloc(b->buf, cap * sizeof *nb);
            if (!nb) { if (spill(b)) return -1; }
            else { b->buf = nb; b->cap = cap; }
        }
    }
    brec *r = &b->buf[b->n++];
    r->am = a_member; r->bm = b_member; r->a = a; r->b = bpos; r->len = len; r->strand = strand ? 1 : 0;
    b->st.input_runs++;
    return 0;
}

/* ---- pre-sorted mode ----
 *
 * The caller promises what the sort would produce: each axis-a member's runs together,
 * in increasing a.  Nothing is buffered; a run that continues the previous one on its
 * diagonal is merged into it, and everything else goes straight to the writer.  The
 * order is fixed at the first run -- order b with A_OVERLAP declared when overlap is
 * allowed, since it cannot be known in advance -- and a broken promise fails the build
 * rather than writing something else. */

static int presorted_flush(lmap_builder *b) {
    if (!b->have_pend) return 0;
    b->have_pend = 0;
    brec *p = &b->pend;
    if (lmap_writer_add_run(b->w, p->am, p->bm, p->a, p->b, p->len, p->strand))
        return fail(b, "writer refused run a_member %u a=%lld len=%lld", p->am, (long long)p->a,
                    (long long)p->len);
    b->st.runs++;
    return 0;
}

static int presorted_start(lmap_builder *b) {
    if (b->started) return 0;
    b->started = 1;
    {
        if (b->allow_overlap) {
            if (lmap_writer_set_order(b->w, LMAP_ORDER_B) || lmap_writer_set_a_overlap(b->w))
                return fail(b, "cannot set the storage order");
            b->st.order = LMAP_ORDER_B;
        } else {
            int order = b->order < 0 ? LMAP_ORDER_A : b->order;
            if (lmap_writer_set_order(b->w, order)) return fail(b, "cannot set the storage order");
            b->st.order = order;
        }
        if (!(b->done_a = calloc(b->n_mem[0] ? b->n_mem[0] : 1, 1))) return fail(b, "out of memory");
    }
    return 0;
}

static int presorted_run(lmap_builder *b, const brec *r) {
    if (presorted_start(b)) return -1;
    if (!b->have_cur_a || r->am != b->cur_a) {
        if (r->am >= b->n_mem[0] || b->done_a[r->am])
            return fail(b, "pre-sorted runs of axis-a member %u are not contiguous", r->am);
        if (b->have_cur_a) b->done_a[b->cur_a] = 1;
        b->cur_a = r->am; b->have_cur_a = 1; b->cur_a_end = 0;
    } else if (r->a < b->cur_a_end) {
        if (!b->allow_overlap)
            return fail(b, "pre-sorted runs of axis-a member %u overlap or go backwards at a=%lld",
                        r->am, (long long)r->a);
        b->st.overlapping_runs++;
    }
    if (r->a + r->len > b->cur_a_end) b->cur_a_end = r->a + r->len;
    brec *p = &b->pend;
    if (b->have_pend && p->am == r->am && p->bm == r->bm && p->strand == r->strand &&
        p->a + p->len == r->a && diag(p) == diag(r)) {
        int64_t d = diag(p);
        p->len += r->len;
        p->b = p->strand ? d - p->a - p->len : d + p->a;
        return 0;
    }
    if (presorted_flush(b)) return -1;
    b->pend = *r; b->have_pend = 1;
    return 0;
}

/* ---- merged, canonical stream: the spill segments plus the (sorted) buffer ---- */

typedef struct { brec *buf; size_t n, i; uint64_t left, next_off; } segsrc;

typedef struct {
    lmap_builder *b;
    size_t seg_cap;                      /* records read ahead per segment */
    segsrc *src; uint32_t nsrc;
    uint32_t *heap; uint32_t nheap;
    brec cur; int have_cur;              /* the run being extended along its diagonal */
} merger;

enum { SEG_BUF = 4096 };

static int seg_fill(merger *m, segsrc *s) {
    size_t k = s->left < m->seg_cap ? (size_t)s->left : m->seg_cap;
    if (fseeko(m->b->spill, (off_t)s->next_off, SEEK_SET) != 0 ||
        fread(s->buf, sizeof *s->buf, k, m->b->spill) != k)
        return fail(m->b, "spill read failed");
    s->n = k; s->i = 0; s->left -= k; s->next_off += k * sizeof *s->buf;
    return 0;
}

static int heap_less(const merger *m, uint32_t x, uint32_t y) {
    const segsrc *p = &m->src[x], *q = &m->src[y];
    int c = cmp_diag(&p->buf[p->i], &q->buf[q->i]);
    return c ? c < 0 : x < y;
}
static void sift_down(merger *m, uint32_t i) {
    for (;;) {
        uint32_t l = 2 * i + 1, r = l + 1, s = i;
        if (l < m->nheap && heap_less(m, m->heap[l], m->heap[s])) s = l;
        if (r < m->nheap && heap_less(m, m->heap[r], m->heap[s])) s = r;
        if (s == i) return;
        uint32_t t = m->heap[i]; m->heap[i] = m->heap[s]; m->heap[s] = t; i = s;
    }
}

static int merger_open(merger *m, lmap_builder *b) {
    memset(m, 0, sizeof *m);
    m->b = b;
    qsort(b->buf, b->n, sizeof *b->buf, qcmp_diag);
    m->nsrc = b->n_seg + 1;
    /* the read-ahead buffers share the memory budget: SEG_BUF records each at most,
     * 64 at least (a budget that small means many segments, and tiny reads) */
    m->seg_cap = b->n_seg ? b->budget / ((size_t)b->n_seg * sizeof(brec)) : SEG_BUF;
    if (m->seg_cap > SEG_BUF) m->seg_cap = SEG_BUF;
    if (m->seg_cap < 64) m->seg_cap = 64;
    m->src = calloc(m->nsrc, sizeof *m->src);
    m->heap = malloc(m->nsrc * sizeof *m->heap);
    if (!m->src || !m->heap) return fail(b, "out of memory");
    for (uint32_t s = 0; s < b->n_seg; s++) {
        segsrc *x = &m->src[s];
        if (!(x->buf = malloc(m->seg_cap * sizeof *x->buf))) return fail(b, "out of memory");
        x->left = b->seg_n[s]; x->next_off = b->seg_off[s];
        if (seg_fill(m, x)) return -1;
    }
    segsrc *mem = &m->src[b->n_seg];     /* the in-memory remainder, already sorted */
    mem->buf = b->buf; mem->n = b->n; mem->i = 0; mem->left = 0;
    for (uint32_t s = 0; s < m->nsrc; s++) if (m->src[s].n > 0) m->heap[m->nheap++] = s;
    for (uint32_t i = m->nheap; i-- > 0; ) sift_down(m, i);
    return 0;
}

/* next record of the raw merge, in diagonal order; 0 at end, -1 on error */
static int merger_raw(merger *m, brec *out) {
    if (m->nheap == 0) return 0;
    uint32_t s = m->heap[0];
    segsrc *x = &m->src[s];
    *out = x->buf[x->i++];
    if (x->i == x->n) {
        if (x->left > 0 && s < m->b->n_seg) { if (seg_fill(m, x)) return -1; }
        else x->n = 0;
    }
    if (x->n == 0 || x->i == x->n) m->heap[0] = m->heap[--m->nheap];
    if (m->nheap) sift_down(m, 0);
    return 1;
}

/* next canonical run: raw runs on one diagonal that overlap or abut, unioned */
static int merger_next(merger *m, brec *out) {
    brec r;
    for (;;) {
        int rc = merger_raw(m, &r);
        if (rc < 0) return -1;
        if (rc == 0) {
            if (!m->have_cur) return 0;
            *out = m->cur; m->have_cur = 0; return 1;
        }
        if (m->have_cur) {
            brec *p = &m->cur;
            if (p->am == r.am && p->bm == r.bm && p->strand == r.strand && diag(p) == diag(&r) &&
                r.a <= p->a + p->len) {
                int64_t d = diag(p), end = r.a + r.len > p->a + p->len ? r.a + r.len : p->a + p->len;
                p->len = end - p->a;
                p->b = p->strand ? d - p->a - p->len : d + p->a;
                continue;
            }
            *out = m->cur; m->cur = r;
            return 1;
        }
        m->cur = r; m->have_cur = 1;
    }
}

static void merger_close(merger *m) {
    if (m->src) for (uint32_t s = 0; s + 1 < m->nsrc; s++) free(m->src[s].buf);
    free(m->src); free(m->heap);
}

/* ---- per a-member: sort on a, count overlap, hand on ---- */

typedef int (*member_sink)(lmap_builder *b, const brec *r, size_t n, void *ctx);

static int for_each_member(lmap_builder *b, member_sink sink, void *ctx) {
    merger m;
    int rc = merger_open(&m, b);
    brec *mem = NULL; size_t n = 0, cap = 0;
    brec r;
    while (rc == 0) {
        int got = merger_next(&m, &r);
        if (got < 0) { rc = -1; break; }
        if (n > 0 && (got == 0 || r.am != mem[0].am)) {
            qsort(mem, n, sizeof *mem, qcmp_a);
            for (size_t i = 1; i < n; i++)
                if (mem[i].a < mem[i-1].a + mem[i-1].len) b->st.overlapping_runs++;
            b->st.runs += n;
            if (sink(b, mem, n, ctx)) { rc = -1; break; }
            n = 0;
        }
        if (got == 0) break;
        if (n == cap) {
            cap = cap ? cap * 2 : 4096;
            brec *nm = realloc(mem, cap * sizeof *nm);
            if (!nm) { rc = fail(b, "out of memory"); break; }
            mem = nm;
        }
        mem[n++] = r;
    }
    free(mem);
    merger_close(&m);
    return rc;
}

static int sink_write(lmap_builder *b, const brec *r, size_t n, void *ctx) {
    (void)ctx;
    for (size_t i = 0; i < n; i++)
        if (lmap_writer_add_run(b->w, r[i].am, r[i].bm, r[i].a, r[i].b, r[i].len, r[i].strand))
            return fail(b, "writer refused run a_member %u a=%lld len=%lld", r[i].am,
                        (long long)r[i].a, (long long)r[i].len);
    return 0;
}

static int sink_file(lmap_builder *b, const brec *r, size_t n, void *ctx) {
    if (fwrite(r, sizeof *r, n, (FILE *)ctx) != n)
        return fail(b, "spill write failed (disk full?): %s", strerror(errno));
    return 0;
}

/* Decide the storage order once overlap is known, and configure the writer. */
static int choose_order(lmap_builder *b) {
    int order = b->order < 0 ? LMAP_ORDER_A : b->order;
    if (b->st.overlapping_runs) {
        if (!b->allow_overlap)
            return fail(b, "%llu runs overlap an earlier run of their axis-a member; axis a "
                        "maps each base at most once unless overlap is allowed",
                        (unsigned long long)b->st.overlapping_runs);
        order = LMAP_ORDER_B;              /* order a cannot encode an overlap */
    }
    if (lmap_writer_set_order(b->w, order) ||
        (b->st.overlapping_runs && lmap_writer_set_a_overlap(b->w)))
        return fail(b, "cannot set the storage order");
    b->st.order = order;
    return 0;
}

int lmap_builder_close(lmap_builder *b, lmap_build_stats *stats) {
    if (!b) return -1;
    int rc = b->failed ? -1 : 0;
    if (b->presorted) {
        if (rc == 0) rc = presorted_start(b);        /* no runs at all: still fix the order */
        if (rc == 0) rc = presorted_flush(b);
    } else if (rc == 0 && b->n_seg == 0) {
        /* everything in memory: canonicalize, sort each member, decide, write */
        merger m;
        size_t n = 0;
        if (merger_open(&m, b) != 0) rc = -1;
        brec r;
        int got;
        while (rc == 0 && (got = merger_next(&m, &r)) != 0) {
            if (got < 0) { rc = -1; break; }
            /* in place: each output consumes at least one input, so the write index
             * never reaches a record the merge has yet to read */
            b->buf[n++] = r;
        }
        merger_close(&m);
        if (rc == 0) {
            b->n = n;
            /* the merge is a_member-major; sort each member's runs on a */
            size_t i = 0;
            while (i < n) {
                size_t j = i;
                while (j < n && b->buf[j].am == b->buf[i].am) j++;
                qsort(b->buf + i, j - i, sizeof *b->buf, qcmp_a);
                for (size_t k = i + 1; k < j; k++)
                    if (b->buf[k].a < b->buf[k-1].a + b->buf[k-1].len) b->st.overlapping_runs++;
                i = j;
            }
            b->st.runs = n;
            if (choose_order(b) == 0) rc = sink_write(b, b->buf, n, NULL);
            else rc = -1;
        }
    } else if (rc == 0) {
        /* spilled: merge into a second temporary file of canonical, member-sorted runs,
         * then decide the order and replay it */
        if (spill(b)) rc = -1;
        FILE *canon = rc == 0 ? temp_file(b) : NULL;
        if (rc == 0 && !canon) rc = fail(b, "cannot create a spill file in %s", b->tmp_dir);
        if (rc == 0) rc = for_each_member(b, sink_file, canon);
        if (rc == 0 && choose_order(b) != 0) rc = -1;
        if (rc == 0) {
            rewind(canon);
            brec chunk[SEG_BUF];
            size_t k;
            while (rc == 0 && (k = fread(chunk, sizeof *chunk, SEG_BUF, canon)) > 0)
                rc = sink_write(b, chunk, k, NULL);
            if (rc == 0 && ferror(canon)) rc = fail(b, "spill read failed");
        }
        if (canon) fclose(canon);
    }
    if (rc == 0 && lmap_writer_close(b->w) != 0) rc = fail(b, "failed writing %s", b->path);
    else if (rc != 0) lmap_writer_abort(b->w);
    b->w = NULL;
    if (stats) {
        *stats = b->st;
        snprintf(stats->error, sizeof stats->error, "%s", rc ? b->err : "");
    }
    lmap_builder_abort(b);
    return rc;
}

int lmap_builder_copy_layout(lmap_builder *b, lmap_file *f) {
    if (!b || !f || b->failed) return -1;
    if (b->n_mem[0] || b->n_mem[1]) return fail(b, "copy_layout needs a fresh builder");
    for (int ax = 0; ax < 2; ax++)
        for (uint32_t i = 0; i < lmap_n_members(f, ax); i++) {
            const lmap_member *m = lmap_member_at(f, ax, i);
            if (lmap_builder_add_member(b, ax, m->name, m->length) != (int32_t)i) return -1;
        }
    for (int ax = 0; ax < 2; ax++)
        for (uint32_t g = 0; g < lmap_n_groups(f, ax); g++) {
            if (lmap_writer_add_group(b->w, ax, lmap_group_name(f, ax, g)) != (int32_t)g)
                return fail(b, "cannot copy group %s", lmap_group_name(f, ax, g));
            const uint32_t *mm; uint32_t n;
            lmap_group_members(f, ax, g, &mm, &n, NULL);
            for (uint32_t k = 0; k < n; k++)
                if (lmap_writer_set_member_group(b->w, ax, mm[k], g)) return fail(b, "cannot copy group");
        }
    for (uint32_t i = 0; i < lmap_meta_count(f); i++) {
        const char *k, *v;
        lmap_meta_at(f, i, &k, &v);
        if (lmap_writer_set_meta(b->w, k, v)) return fail(b, "cannot copy metadata %s", k);
    }
    for (uint32_t i = 0; i < lmap_n_app_sections(f); i++) {
        const char *id = lmap_app_section_id(f, i);
        uint8_t *p; size_t n;
        if (lmap_read_section(f, id, &p, &n) != 0) return fail(b, "cannot read section %s", id);
        int rc = lmap_writer_add_section(b->w, id, p, n);
        free(p);
        if (rc) return fail(b, "cannot copy section %s", id);
    }
    return 0;
}

/* ---- coarsening ---- */

typedef struct { int64_t k, o, len; uint8_t strand; } chain;

/* Clip a chain at the other member's end and hand it to the builder. */
static int emit_chain(lmap_builder *b, int K, uint32_t km, uint32_t om, chain c,
                      uint64_t olen, lmap_coarsen_stats *st) {
    if ((uint64_t)(c.o + c.len) > olen) {
        int64_t e = (int64_t)((uint64_t)(c.o + c.len) - olen);
        if (e >= c.len) { st->dropped++; st->clipped_bp += (uint64_t)c.len; return 0; }
        /* forward: the far end is the tail; reverse: it pairs with the key start */
        if (c.strand) c.k += e;
        c.len -= e; st->clipped_bp += (uint64_t)e;
    }
    st->chains_out++;
    return K ? lmap_builder_add_run(b, om, km, c.o, c.k, c.len, c.strand)
             : lmap_builder_add_run(b, km, om, c.k, c.o, c.len, c.strand);
}

int lmap_coarsen(lmap_file *f, lmap_builder *b, int K, int64_t max_gap,
                 lmap_coarsen_stats *stats) {
    lmap_coarsen_stats st = {0};
    if (!f || !b || b->failed || (K != 0 && K != 1) || max_gap < 0) return -1;
    int O = 1 - K;
    uint32_t no = lmap_n_members(f, O);
    chain *act = malloc((no ? no : 1) * sizeof *act);
    uint8_t *on = calloc(no ? no : 1, 1);
    uint32_t *list = malloc((no ? no : 1) * sizeof *list);   /* members with an open chain */
    if (!act || !on || !list) { free(act); free(on); free(list); return fail(b, "out of memory"); }
    lmap_builder_allow_overlap(b, 1);
    int rc = 0;
    for (uint32_t km = 0; rc == 0 && km < lmap_n_members(f, K); km++) {
        lmap_cursor *c = lmap_cursor_open(f, K, km, NULL, 0, 0);
        if (!c) { rc = fail(b, "cannot open a cursor"); break; }
        uint32_t nlist = 0;
        lmap_hit h;
        int got;
        rc = lmap_cursor_seek(c, 0, (int64_t)lmap_member_at(f, K, km)->length);
        while (rc == 0 && (got = lmap_cursor_next(c, &h)) != 0) {
            if (got < 0) { rc = fail(b, "corrupt chunk while coarsening"); break; }
            st.runs_in++;
            uint32_t om = K ? h.a_member : h.b_member;
            int64_t k = K ? h.b : h.a, o = K ? h.a : h.b;
            chain *ch = &act[om];
            if (on[om] && ch->strand == h.strand) {
                int64_t kg = k - (ch->k + ch->len);
                int64_t og = h.strand ? ch->o - (o + h.len) : o - (ch->o + ch->len);
                if (kg >= 0 && kg <= max_gap && og >= 0 && og <= max_gap) {
                    ch->len = k + h.len - ch->k;
                    if (h.strand) ch->o = o;
                    continue;
                }
            }
            if (on[om]) {
                if (emit_chain(b, K, km, om, *ch, lmap_member_at(f, O, om)->length, &st)) { rc = -1; break; }
            } else { on[om] = 1; list[nlist++] = om; }
            ch->k = k; ch->o = o; ch->len = h.len; ch->strand = h.strand;
        }
        lmap_cursor_close(c);
        for (uint32_t i = 0; i < nlist; i++) {           /* close this key member's chains */
            uint32_t om = list[i];
            if (rc == 0 && emit_chain(b, K, km, om, act[om], lmap_member_at(f, O, om)->length, &st)) rc = -1;
            on[om] = 0;
        }
    }
    free(act); free(on); free(list);
    if (rc == 0) {
        const char *old = lmap_meta_get(f, "max_gap");
        int64_t prev = old ? strtoll(old, NULL, 10) : 0;
        char v[32];
        snprintf(v, sizeof v, "%lld", (long long)(prev > max_gap ? prev : max_gap));
        if (lmap_writer_set_meta(b->w, "max_gap", v)) rc = fail(b, "cannot set max_gap");
    }
    if (stats) *stats = st;
    return rc == 0 ? 0 : -1;
}

void lmap_builder_abort(lmap_builder *b) {
    if (!b) return;
    if (b->w) lmap_writer_abort(b->w);
    if (b->spill) fclose(b->spill);
    free(b->len[0]); free(b->len[1]);
    free(b->buf); free(b->seg_off); free(b->seg_n); free(b->done_a);
    free(b->path); free(b->tmp_dir);
    free(b);
}
