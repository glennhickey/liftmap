/* lmap_chain -- chaining and the overlap-based paralogy filter.  See liftmap.h.
 *
 * A port of taffy's chain.c (itself paffy's chaining.c, generalised), reproducing its
 * results exactly -- partition order, tie-breaks, chain-id numbering -- but over plain
 * spans with member ids, without sonLib, and without reordering the caller's array.
 *
 * Per (q_member, strand) partition, spans are swept in q order (reverse-strand q
 * coordinates mirrored so the sweep runs forward).  Each span joins the best-scoring
 * active chain it can extend -- one ending before it on both axes, gaps within max_gap,
 * gap cost below the span's own score -- or starts a chain.  Chains are then claimed in
 * score order: a chain walks back through its links, and is cut where a link was already
 * claimed by a better chain.  The active set is ordered by (t_member, t_end, q_end, ...)
 * so a span's candidate predecessors are found by one search and a backward walk.
 */
#include "liftmap.h"

#include <stdlib.h>
#include <string.h>

static int icmp(int64_t a, int64_t b) { return a > b ? 1 : (a < b ? -1 : 0); }

/* A span as the sweep sees it: q mirrored on reverse, and its input index. */
typedef struct {
    int64_t qs, qe, ts, te, score;
    uint32_t qm, tm;
    int rev;
    size_t idx;
} cspan;

/* taffy orders partitions with the reverse strand first (its strand -1 < +1); kept so
 * results, including chain ids, match it. */
static int cmp_sweep(const void *x, const void *y) {
    const cspan *a = x, *b = y;
    int c = icmp(a->qm, b->qm); if (c) return c;
    if (a->rev != b->rev) return a->rev ? -1 : 1;
    if ((c = icmp(a->qs, b->qs))) return c;
    if ((c = icmp(a->tm, b->tm))) return c;
    if ((c = icmp(a->ts, b->ts))) return c;
    if ((c = icmp(a->te, b->te))) return c;
    if ((c = icmp(a->qe, b->qe))) return c;
    return icmp((int64_t)a->idx, (int64_t)b->idx);
}

/* ---- the active set: a treap over sweep positions, keyed (tm, te, qe, qs, ts, pos) ---- */

typedef struct { int32_t l, r; uint32_t pri; } tnode;
typedef struct { const cspan *s; tnode *t; int32_t root; } treap;

static int key_cmp(const treap *T, int32_t x, int32_t y) {
    const cspan *a = &T->s[x], *b = &T->s[y];
    int c;
    if ((c = icmp(a->tm, b->tm))) return c;
    if ((c = icmp(a->te, b->te))) return c;
    if ((c = icmp(a->qe, b->qe))) return c;
    if ((c = icmp(a->qs, b->qs))) return c;
    if ((c = icmp(a->ts, b->ts))) return c;
    return icmp(x, y);
}
static void split(treap *T, int32_t n, int32_t x, int32_t *l, int32_t *r) {   /* l: keys < x */
    if (n < 0) { *l = *r = -1; return; }
    if (key_cmp(T, n, x) < 0) { split(T, T->t[n].r, x, &T->t[n].r, r); *l = n; }
    else { split(T, T->t[n].l, x, l, &T->t[n].l); *r = n; }
}
static int32_t merge(treap *T, int32_t a, int32_t b) {
    if (a < 0) return b;
    if (b < 0) return a;
    if (T->t[a].pri > T->t[b].pri) { T->t[a].r = merge(T, T->t[a].r, b); return a; }
    T->t[b].l = merge(T, a, T->t[b].l); return b;
}
static void t_insert(treap *T, int32_t x) {
    int32_t l, r;
    T->t[x].l = T->t[x].r = -1;
    split(T, T->root, x, &l, &r);
    T->root = merge(T, merge(T, l, x), r);
}
static int32_t t_remove(treap *T, int32_t n, int32_t x) {
    if (n < 0) return -1;
    if (n == x) return merge(T, T->t[n].l, T->t[n].r);
    if (key_cmp(T, x, n) < 0) T->t[n].l = t_remove(T, T->t[n].l, x);
    else T->t[n].r = t_remove(T, T->t[n].r, x);
    return n;
}
/* Walk the nodes with (tm, te) <= (tm, ts) in descending key order: every predecessor
 * candidate, nearest first.  A stack of the unvisited nodes at or below the bound makes
 * each step O(1) amortised (the tree is not modified during a walk). */
typedef struct { int32_t *st; int32_t n; } rwalk;
static void rwalk_push_le(const treap *T, rwalk *w, int32_t n, uint32_t tm, int64_t ts) {
    while (n >= 0) {
        const cspan *a = &T->s[n];
        if (a->tm < tm || (a->tm == tm && a->te <= ts)) { w->st[w->n++] = n; n = T->t[n].r; }
        else n = T->t[n].l;
    }
}
static int32_t rwalk_next(const treap *T, rwalk *w) {
    if (w->n == 0) return -1;
    int32_t x = w->st[--w->n];
    for (int32_t m = T->t[x].l; m >= 0; m = T->t[m].r) w->st[w->n++] = m;   /* its left subtree */
    return x;
}

typedef struct { int64_t score; int32_t prev; } dnode;

static int64_t gap_cost(const lmap_chain_params *p, int64_t qg, int64_t tg) {
    return qg + tg == 0 ? 0 : p->chain_open + p->chain_extend * (qg + tg);
}

/* claim order: score, then (qs, ts, qe, te), then position -- taffy's pop order */
typedef struct { const cspan *s; const dnode *d; } rank_ctx;
static int cmp_rank(const void *x, const void *y, const rank_ctx *rc) {
    int32_t a = *(const int32_t *)x, b = *(const int32_t *)y;
    int c;
    if ((c = icmp(rc->d[a].score, rc->d[b].score))) return c;
    const cspan *p = &rc->s[a], *q = &rc->s[b];
    if ((c = icmp(p->qs, q->qs))) return c;
    if ((c = icmp(p->ts, q->ts))) return c;
    if ((c = icmp(p->qe, q->qe))) return c;
    if ((c = icmp(p->te, q->te))) return c;
    return icmp(a, b);
}
/* insertion-free merge sort with a context (qsort has none, and a global would race) */
static void msort_rank(int32_t *v, int32_t *tmp, size_t n, const rank_ctx *rc) {
    if (n < 2) return;
    size_t h = n / 2;
    msort_rank(v, tmp, h, rc); msort_rank(v + h, tmp, n - h, rc);
    size_t i = 0, j = h, k = 0;
    while (i < h && j < n) tmp[k++] = cmp_rank(&v[i], &v[j], rc) <= 0 ? v[i++] : v[j++];
    while (i < h) tmp[k++] = v[i++];
    while (j < n) tmp[k++] = v[j++];
    memcpy(v, tmp, n * sizeof *v);
}

static int cmp_info(const void *x, const void *y) {
    const lmap_chain_info *a = x, *b = y;
    if (a->score != b->score) return a->score > b->score ? -1 : 1;
    return icmp(a->id, b->id);
}

int lmap_chain(const lmap_span *in, size_t n, const lmap_chain_params *params,
               int64_t *chain_id, lmap_chain_info **chains, size_t *n_chains) {
    if (!chains || !n_chains || (n && (!in || !chain_id))) return -1;
    *chains = NULL; *n_chains = 0;
    if (n == 0) return 0;
    if (n > INT32_MAX) return -1;
    lmap_chain_params p = params ? *params : (lmap_chain_params)LMAP_CHAIN_PARAMS_DEFAULT;
    cspan *s = malloc(n * sizeof *s);
    tnode *tn = malloc(n * sizeof *tn);
    dnode *d = malloc(n * sizeof *d);
    uint8_t *taken = calloc(n, 1);
    int32_t *ord = malloc(n * sizeof *ord), *tmp = malloc(n * sizeof *tmp);
    int32_t *drop = malloc(n * sizeof *drop), *walk = malloc(n * sizeof *walk);
    lmap_chain_info *info = malloc(n * sizeof *info);
    if (!s || !tn || !d || !taken || !ord || !tmp || !drop || !walk || !info) {
        free(s); free(tn); free(d); free(taken); free(ord); free(tmp); free(drop); free(walk); free(info);
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        const lmap_span *x = &in[i];
        cspan *c = &s[i];
        c->rev = x->strand != 0;
        c->qs = c->rev ? -x->q_end : x->q_start;
        c->qe = c->rev ? -x->q_start : x->q_end;
        c->ts = x->t_start; c->te = x->t_end;
        c->qm = x->q_member; c->tm = x->t_member;
        c->score = x->score; c->idx = i;
    }
    qsort(s, n, sizeof *s, cmp_sweep);
    for (size_t i = 0; i < n; i++) {                   /* deterministic treap priorities: */
        uint64_t z = (uint64_t)i + 0x9E3779B97F4A7C15ull;  /* splitmix64, so they behave as */
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;       /* random -- a plain multiplicative */
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;       /* hash of 0,1,2.. is correlated */
        tn[i].pri = (uint32_t)((z ^ (z >> 31)) >> 32);     /* and left the tree deep */
    }
    int64_t next_id = 1;
    size_t ni = 0;
    for (size_t lo = 0; lo < n; ) {                    /* one partition: (qm, strand) */
        size_t hi = lo + 1;
        while (hi < n && s[hi].qm == s[lo].qm && s[hi].rev == s[lo].rev) hi++;
        treap T = { s, tn, -1 };
        for (size_t i = lo; i < hi; i++) {
            const cspan *a = &s[i];
            d[i].score = a->score; d[i].prev = -1;
            size_t nd = 0;
            rwalk rw = { walk, 0 };
            rwalk_push_le(&T, &rw, T.root, a->tm, a->ts);
            for (int32_t pc = rwalk_next(&T, &rw); pc >= 0; pc = rwalk_next(&T, &rw)) {
                const cspan *b = &s[pc];
                if (b->tm != a->tm) break;              /* the key leads with tm */
                if (a->qs < b->qe) continue;            /* q overlap: skip, keep scanning */
                int64_t qg = a->qs - b->qe;
                if (qg > p.max_gap) { drop[nd++] = pc; continue; }   /* never extendable again */
                if (a->ts < b->te) continue;
                int64_t tg = a->ts - b->te;
                if (tg > p.max_gap) break;              /* earlier te only widens the gap */
                int64_t g = gap_cost(&p, qg, tg), cand = a->score + d[pc].score - g;
                if (g < a->score && cand > d[i].score) { d[i].score = cand; d[i].prev = pc; }
            }
            t_insert(&T, (int32_t)i);
            for (size_t k = 0; k < nd; k++) T.root = t_remove(&T, T.root, drop[k]);
        }
        /* claim chains best first; a link already claimed cuts the chain there */
        size_t m = hi - lo;
        for (size_t k = 0; k < m; k++) ord[k] = (int32_t)(lo + k);
        rank_ctx rc = { s, d };
        msort_rank(ord, tmp, m, &rc);
        for (size_t k = m; k-- > 0; ) {
            int32_t root = ord[k];
            if (taken[root]) continue;
            int64_t cid = next_id++, bp = 0, cnt = 0, total = d[root].score;
            for (int32_t w = root; ; ) {
                taken[w] = 1;
                chain_id[s[w].idx] = cid;
                bp += s[w].qe - s[w].qs; cnt++;
                int32_t pv = d[w].prev;
                if (pv < 0) break;
                if (taken[pv]) {
                    /* the chain's score included everything below the cut; take it out:
                     * w's contribution falls from score + prev - gap to its own score */
                    total -= d[pv].score - gap_cost(&p, s[w].qs - s[pv].qe, s[w].ts - s[pv].te);
                    break;
                }
                w = pv;
            }
            info[ni].id = cid; info[ni].score = total; info[ni].bp = bp; info[ni].n_spans = cnt;
            ni++;
        }
        lo = hi;
    }
    qsort(info, ni, sizeof *info, cmp_info);
    free(s); free(tn); free(d); free(taken); free(ord); free(tmp); free(drop); free(walk);
    *chains = info; *n_chains = ni;
    return 0;
}

/* ---- paralogy filter ---- */

typedef struct { int64_t lo, hi; } qiv;
static int cmp_qiv(const void *x, const void *y) {
    int64_t a = ((const qiv *)x)->lo, b = ((const qiv *)y)->lo;
    return (a > b) - (a < b);
}

int lmap_chain_select(const lmap_span *s, size_t n, const int64_t *chain_id,
                      const lmap_chain_info *chains, size_t n_chains, double overlap_frac,
                      size_t cap, uint8_t *keep, size_t keep_len) {
    if (!keep || (n && (!s || !chain_id)) || (n_chains && !chains)) return -1;
    if (n == 0 || n_chains == 0) return 0;
    int64_t max_id = 0;
    for (size_t i = 0; i < n; i++) {
        if (chain_id[i] < 1) return -1;
        if (chain_id[i] > max_id) max_id = chain_id[i];
    }
    if ((size_t)max_id >= keep_len) return -1;
    /* spans bucketed by chain */
    size_t *off = calloc((size_t)max_id + 2, sizeof *off), *bucket = malloc(n * sizeof *bucket);
    size_t *cur = malloc(((size_t)max_id + 2) * sizeof *cur);
    qiv *kept = malloc(n * sizeof *kept), *next = malloc(n * sizeof *next), *cand = malloc(n * sizeof *cand);
    if (!off || !bucket || !cur || !kept || !next || !cand) {
        free(off); free(bucket); free(cur); free(kept); free(next); free(cand); return -1;
    }
    for (size_t i = 0; i < n; i++) off[chain_id[i] + 1]++;
    for (int64_t k = 1; k <= max_id + 1; k++) off[k] += off[k - 1];
    memcpy(cur, off, ((size_t)max_id + 2) * sizeof *cur);
    for (size_t i = 0; i < n; i++) bucket[cur[chain_id[i]]++] = i;
    size_t nkept = 0, accepted = 0, limit = cap ? cap : n_chains;
    for (size_t k = 0; k < n_chains && accepted < limit; k++) {
        int64_t cid = chains[k].id;
        if (cid < 1 || cid > max_id) continue;
        size_t nc = 0;
        for (size_t z = off[cid]; z < off[cid + 1]; z++) {
            cand[nc].lo = s[bucket[z]].q_start; cand[nc].hi = s[bucket[z]].q_end; nc++;
        }
        if (nc == 0) continue;
        qsort(cand, nc, sizeof *cand, cmp_qiv);        /* reverse chains arrive q-descending */
        size_t w = 0;
        for (size_t r = 0; r < nc; r++) {
            if (w && cand[r].lo <= cand[w-1].hi) { if (cand[r].hi > cand[w-1].hi) cand[w-1].hi = cand[r].hi; }
            else cand[w++] = cand[r];
        }
        nc = w;
        int64_t cbp = 0, obp = 0;
        for (size_t i = 0; i < nc; i++) cbp += cand[i].hi - cand[i].lo;
        if (cbp == 0) continue;
        for (size_t i = 0, j = 0; i < nkept && j < nc; ) {
            int64_t lo = kept[i].lo > cand[j].lo ? kept[i].lo : cand[j].lo;
            int64_t hi = kept[i].hi < cand[j].hi ? kept[i].hi : cand[j].hi;
            if (hi > lo) obp += hi - lo;
            if (kept[i].hi < cand[j].hi) i++; else j++;
        }
        if ((double)obp > overlap_frac * (double)cbp) continue;    /* a paralog of kept: drop */
        size_t i = 0, j = 0, o = 0;                                 /* kept |= cand */
        while (i < nkept || j < nc) {
            qiv pk = (j >= nc || (i < nkept && kept[i].lo <= cand[j].lo)) ? kept[i++] : cand[j++];
            if (o && pk.lo <= next[o-1].hi) { if (pk.hi > next[o-1].hi) next[o-1].hi = pk.hi; }
            else next[o++] = pk;
        }
        qiv *t = kept; kept = next; next = t; nkept = o;
        keep[cid] = 1;
        accepted++;
    }
    free(off); free(bucket); free(cur); free(kept); free(next); free(cand);
    return 0;
}
