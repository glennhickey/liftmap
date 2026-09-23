/* query -- check lmap_query_a / lmap_query_b and the cursor against a brute-force oracle.
 *
 * in: binary {int64 a; int64 b; int64 len; int64 strand; int32 amem; int32 bmem} x N,
 *     grouped by amem with a increasing within a member.
 * Axis-a member ids are remapped (reversed) before writing, so the input no longer
 * arrives in id order and the writer's dir.a sort is exercised.
 * The oracle maps a window to the other axis via the explicit base pairing
 * pair(i) = strand ? other+len-1-i : other+i, not via the library's clip formula. */
#include "../src/lmap_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int64_t a,b,len,strand; int32_t am,bm; } inrec;

static int cmp_a(const void *x, const void *y) {
    const lmap_hit *p=x,*q=y; return (p->a>q->a)-(p->a<q->a);
}
static int cmp_b(const void *x, const void *y) {
    const lmap_hit *p=x,*q=y;
    if (p->b!=q->b) return (p->b>q->b)-(p->b<q->b);
    if (p->a_member!=q->a_member) return (p->a_member>q->a_member)-(p->a_member<q->a_member);
    return (p->a>q->a)-(p->a<q->a);
}

/* cursor order: key position, other member, other position, len, strand */
static int cur_axis;
static int cmp_cur(const void *x, const void *y) {
    const lmap_hit *p=x,*q=y; int K=cur_axis;
    int64_t kp=K?p->b:p->a, kq=K?q->b:q->a;
    if (kp!=kq) return (kp>kq)-(kp<kq);
    uint32_t mp=K?p->a_member:p->b_member, mq=K?q->a_member:q->b_member;
    if (mp!=mq) return (mp>mq)-(mp<mq);
    int64_t op=K?p->a:p->b, oq=K?q->a:q->b;
    if (op!=oq) return (op>oq)-(op<oq);
    if (p->len!=q->len) return (p->len>q->len)-(p->len<q->len);
    return (int)p->strand-(int)q->strand;
}
static int same_hit(const lmap_hit *x, const lmap_hit *y) {
    return x->a==y->a && x->b==y->b && x->len==y->len && x->strand==y->strand &&
           x->a_member==y->a_member && x->b_member==y->b_member;
}
static void as_hit(const inrec *r, lmap_hit *h) {
    h->a=r->a; h->b=r->b; h->len=r->len; h->strand=(uint8_t)r->strand;
    h->a_member=(uint32_t)r->am; h->b_member=(uint32_t)r->bm;
}

/* oracle clip, derived independently of the library */
static int oracle(const inrec *r, int axis, int64_t lo, int64_t hi, lmap_hit *h) {
    int64_t s0 = axis ? r->b : r->a, o0 = axis ? r->a : r->b;
    int64_t s = s0 > lo ? s0 : lo, e = s0 + r->len < hi ? s0 + r->len : hi;
    if (s >= e) return 0;
    int64_t i1 = s - s0, i2 = e - 1 - s0;               /* first and last offsets kept */
    int64_t p1 = r->strand ? o0 + r->len - 1 - i1 : o0 + i1;
    int64_t p2 = r->strand ? o0 + r->len - 1 - i2 : o0 + i2;
    int64_t other = p1 < p2 ? p1 : p2;
    if (axis) { h->b = s; h->a = other; } else { h->a = s; h->b = other; }
    h->len = e - s; h->strand = (uint8_t)r->strand;
    h->a_member = (uint32_t)r->am; h->b_member = (uint32_t)r->bm;
    return 1;
}

static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

int main(int argc, char **argv) {
    if (argc != 4 && argc != 5) { fprintf(stderr,"usage: query runs.bin out.lmap nqueries [a|b]\n"); return 2; }
    int order = (argc == 5 && argv[4][0] == 'b') ? LMAP_ORDER_B : LMAP_ORDER_A;
    FILE *fp = fopen(argv[1],"rb"); if (!fp) { perror("runs"); return 2; }
    fseek(fp,0,SEEK_END); long sz=ftell(fp); fseek(fp,0,SEEK_SET);
    long n = sz/(long)sizeof(inrec);
    inrec *r = malloc((size_t)n*sizeof *r);
    if (fread(r,sizeof *r,(size_t)n,fp)!=(size_t)n) return 2;
    fclose(fp);
    int32_t na=0, nb=0;
    for (long i=0;i<n;i++){ if(r[i].am+1>na) na=r[i].am+1; if(r[i].bm+1>nb) nb=r[i].bm+1; }
    for (long i=0;i<n;i++) r[i].am = na-1-r[i].am;          /* remap: ids no longer in input order */

    lmap_writer *w = lmap_writer_open(argv[2],"hal2.edge","liftmap 1\n");
    if (lmap_writer_set_order(w, order) != 0) { fprintf(stderr,"set_order failed\n"); return 1; }
    char nm[64];
    for (int32_t i=0;i<na;i++){ snprintf(nm,sizeof nm,"a%d",i); lmap_writer_add_member(w,0,nm,1ULL<<40); }
    for (int32_t i=0;i<nb;i++){ snprintf(nm,sizeof nm,"b%d",i); lmap_writer_add_member(w,1,nm,1ULL<<40); }
    for (long i=0;i<n;i++)
        if (lmap_writer_add_run(w,(uint32_t)r[i].am,(uint32_t)r[i].bm,r[i].a,r[i].b,r[i].len,(uint8_t)r[i].strand)) {
            fprintf(stderr,"add_run %ld failed\n",i); return 1; }
    if (lmap_writer_close(w)) { fprintf(stderr,"close failed\n"); return 1; }
    lmap_file *f = lmap_open(argv[2]);
    if (!f) { fprintf(stderr,"open failed (invariant check?)\n"); return 1; }
    if (lmap_order(f) != order) { fprintf(stderr,"order not round-tripped\n"); return 1; }

    /* index input runs by member on each axis, for the oracle */
    long *ca = calloc((size_t)na+1,sizeof *ca), *cb = calloc((size_t)nb+1,sizeof *cb);
    for (long i=0;i<n;i++){ ca[r[i].am+1]++; cb[r[i].bm+1]++; }
    for (int32_t i=0;i<na;i++) ca[i+1]+=ca[i];
    for (int32_t i=0;i<nb;i++) cb[i+1]+=cb[i];
    long *ia = malloc((size_t)n*sizeof *ia), *ib = malloc((size_t)n*sizeof *ib);
    long *fa = calloc((size_t)na,sizeof *fa), *fb = calloc((size_t)nb,sizeof *fb);
    for (long i=0;i<n;i++){ ia[ca[r[i].am]+fa[r[i].am]++]=i; ib[cb[r[i].bm]+fb[r[i].bm]++]=i; }

    const int64_t widths[] = { 1, 37, 1000, 100000, 5000000 };
    const int NW = 5;
    long nq = atol(argv[3]), fails = 0, total_hits = 0;
    double dec_sum[2][5] = {{0}}; long dec_n[2][5] = {{0}};
    long expcap = 1L << 20;
    lmap_hit *exp = malloc((size_t)expcap * sizeof *exp);

    for (long q=0; q<nq; q++) {
        long i = (long)(xr() % (uint64_t)n);
        int axis = (int)(xr() & 1);
        int wi = (int)(xr() % NW);
        int64_t c = (axis ? r[i].b : r[i].a) + (int64_t)(xr() % (uint64_t)r[i].len);
        int64_t lo = c - widths[wi], hi = c + widths[wi];
        if (lo < 0) lo = 0;
        if (q % 7 == 0) { lo = axis ? r[i].b : r[i].a; hi = lo + r[i].len; }   /* exact run */
        if (q % 11 == 0) { hi = lo; }                                        /* empty window */
        uint32_t m = (uint32_t)(axis ? r[i].bm : r[i].am);

        long ne = 0;
        long *idx = axis ? ib : ia; long *cc = axis ? cb : ca;
        for (long k = cc[m]; k < cc[m+1]; k++) {
            if (ne == expcap) { expcap *= 2; exp = realloc(exp, (size_t)expcap*sizeof *exp); }
            if (oracle(&r[idx[k]], axis, lo, hi, &exp[ne])) ne++;
        }
        qsort(exp, (size_t)ne, sizeof *exp, axis ? cmp_b : cmp_a);

        lmap_hit *got; size_t ng; lmap_query_stats st;
        int rc = axis ? lmap_query_b(f,m,lo,hi,&got,&ng,&st) : lmap_query_a(f,m,lo,hi,&got,&ng,&st);
        int ok = rc == 0 && (long)ng == ne;
        for (long k = 0; ok && k < ne; k++)
            ok = got[k].a==exp[k].a && got[k].b==exp[k].b && got[k].len==exp[k].len &&
                 got[k].strand==exp[k].strand && got[k].a_member==exp[k].a_member &&
                 got[k].b_member==exp[k].b_member;
        if (!ok) {
            if (fails < 5) fprintf(stderr,"MISMATCH q=%ld axis=%c m=%u [%lld,%lld): rc=%d got %zu want %ld\n",
                                   q, axis?'b':'a', m, (long long)lo, (long long)hi, rc, ng, ne);
            fails++;
        }
        if (hi > lo && q % 7 && q % 11) { dec_sum[axis][wi] += st.chunks_decoded; dec_n[axis][wi]++; }
        total_hits += ne;
        free(got);
    }
    printf("order %c: %ld queries, %ld hits compared, %ld mismatches\n", order ? 'b' : 'a', nq, total_hits, fails);
    printf("mean chunks decoded per query, by half-width:\n  %-9s", "");
    for (int k=0;k<NW;k++) printf("%10lld", (long long)widths[k]);
    for (int ax=0; ax<2; ax++) {
        printf("\n  axis %c   ", ax?'b':'a');
        for (int k=0;k<NW;k++) printf("%10.2f", dec_n[ax][k] ? dec_sum[ax][k]/dec_n[ax][k] : 0);
    }
    printf("\n");

    /* ---- cursor: whole runs in key order, other-axis filters, point lookups, reuse */
    long cfails = 0, chits = 0, cq = nq / 2 + 1;
    uint8_t *keep = malloc((size_t)(na > nb ? na : nb));
    uint32_t *flt = malloc((size_t)(na > nb ? na : nb) * sizeof *flt);
    lmap_hit *got = malloc((size_t)expcap * sizeof *got);
    long gotcap = expcap;
    for (long q=0; q<cq; q++) {
        long i = (long)(xr() % (uint64_t)n);
        int axis = (int)(xr() & 1), O = 1 - axis;
        uint32_t m = (uint32_t)(axis ? r[i].bm : r[i].am);
        int32_t no = O ? nb : na;
        int filtered = (int)(xr() % 3) == 0;
        uint32_t nf = 0;
        if (filtered) {                             /* keep ~half the other members, and i's */
            for (int32_t k=0;k<no;k++) { keep[k] = (uint8_t)(xr() & 1); }
            keep[O ? r[i].bm : r[i].am] = 1;
            for (int32_t k=0;k<no;k++) if (keep[k]) flt[nf++] = (uint32_t)k;
        }
        size_t budget = (size_t[]){ 0, 1u << 20, LMAP_CACHE_ALL }[xr() % 3];
        lmap_cursor *c = lmap_cursor_open(f, axis, m, filtered ? flt : NULL, nf, budget);
        if (!c) { fprintf(stderr,"cursor_open failed\n"); return 1; }
        long *idx = axis ? ib : ia; long *cc = axis ? cb : ca;
        for (int pass = 0; pass < 2; pass++) {       /* two seeks on one cursor */
            int wi = (int)(xr() % NW);
            int64_t ctr = (axis ? r[i].b : r[i].a) + (int64_t)(xr() % (uint64_t)r[i].len);
            int64_t lo = ctr - widths[wi], hi = ctr + widths[wi];
            if (lo < 0) lo = 0;
            long ne = 0;
            for (long k = cc[m]; k < cc[m+1]; k++) {
                const inrec *x = &r[idx[k]];
                int64_t s0 = axis ? x->b : x->a;
                if (s0 >= hi || s0 + x->len <= lo) continue;
                if (filtered && !keep[O ? x->bm : x->am]) continue;
                if (ne == expcap) { expcap *= 2; exp = realloc(exp, (size_t)expcap*sizeof *exp); }
                as_hit(x, &exp[ne++]);
            }
            cur_axis = axis;
            qsort(exp, (size_t)ne, sizeof *exp, cmp_cur);
            long ng = 0; int rc = lmap_cursor_seek(c, lo, hi);
            lmap_hit h;
            while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
                rc = 0;
                if (ng == gotcap) { gotcap *= 2; got = realloc(got, (size_t)gotcap*sizeof *got); }
                got[ng++] = h;
            }
            int ok = rc == 0 && ng == ne;
            for (long k = 0; ok && k < ne; k++) ok = same_hit(&got[k], &exp[k]);
            if (!ok) {
                if (cfails < 5) fprintf(stderr,"CURSOR MISMATCH q=%ld axis=%c m=%u filt=%d [%lld,%lld): got %ld want %ld\n",
                                        q, axis?'b':'a', m, filtered, (long long)lo, (long long)hi, ng, ne);
                cfails++;
            }
            chits += ne;
        }
        /* point lookups */
        for (int t = 0; t < 4; t++) {
            int64_t pos = (axis ? r[i].b : r[i].a) + (int64_t)(xr() % (uint64_t)r[i].len);
            long ne = 0;
            for (long k = cc[m]; k < cc[m+1]; k++) {
                const inrec *x = &r[idx[k]];
                int64_t s0 = axis ? x->b : x->a;
                if (pos < s0 || pos >= s0 + x->len) continue;
                if (filtered && !keep[O ? x->bm : x->am]) continue;
                as_hit(x, &exp[ne++]);
            }
            cur_axis = axis;
            qsort(exp, (size_t)ne, sizeof *exp, cmp_cur);
            int np = lmap_cursor_point(c, pos, got, (int)gotcap);
            int ok = np == ne;
            for (long k = 0; ok && k < ne; k++) ok = same_hit(&got[k], &exp[k]);
            if (!ok) {
                if (cfails < 5) fprintf(stderr,"POINT MISMATCH q=%ld axis=%c m=%u pos=%lld: got %d want %ld\n",
                                        q, axis?'b':'a', m, (long long)pos, np, ne);
                cfails++;
            }
            if (np > 0 && np > 1) {              /* a short buffer still reports the total */
                int n1 = lmap_cursor_point(c, pos, got, 1);
                if (n1 != np) { if (cfails < 5) fprintf(stderr,"POINT cap=1 returned %d, want %d\n", n1, np); cfails++; }
            }
        }
        lmap_cursor_close(c);
    }

    /* streaming full scans: every run of the member, in order, holding only what overlaps */
    long sfails = 0; size_t worst_peak = 0, worst_total = 0;
    for (int axis = 0; axis < 2; axis++) {
        uint32_t nm = lmap_n_members(f, axis);
        for (uint32_t t = 0; t < 4 && t < nm; t++) {
            uint32_t m = (uint32_t)(xr() % nm);
            long want = (axis ? cb : ca)[m+1] - (axis ? cb : ca)[m];
            lmap_cursor *c = lmap_cursor_open(f, axis, m, NULL, 0, 0);
            lmap_hit h, prev = {0}; long ng = 0; int rc = lmap_cursor_seek(c, 0, (int64_t)1 << 40), ordered = 1;
            while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
                rc = 0;
                cur_axis = axis;
                if (ng && cmp_cur(&prev, &h) > 0) ordered = 0;
                prev = h; ng++;
            }
            lmap_cursor_stats cs; lmap_cursor_get_stats(c, &cs);
            if (rc != 0 || ng != want || !ordered || cs.cached_bytes != 0) {
                fprintf(stderr,"SCAN axis=%c m=%u: rc=%d got %ld want %ld ordered=%d cached_after=%zu\n",
                        axis?'b':'a', m, rc, ng, want, ordered, cs.cached_bytes);
                sfails++;
            }
            if (cs.peak_cached_bytes > worst_peak) { worst_peak = cs.peak_cached_bytes; worst_total = (size_t)want * 32; }
            lmap_cursor_close(c);
        }
    }
    printf("cursor: %ld windows+points checked (%ld hits), %ld mismatches; 8 streaming scans, %ld failures; "
           "largest scan peak cache %zu bytes for %zu bytes of runs\n",
           cq, chits, cfails, sfails, worst_peak, worst_total);
    fails += cfails + sfails;
    free(got); free(keep); free(flt);
    lmap_close(f);
    free(exp); free(ia); free(ib); free(fa); free(fb); free(ca); free(cb); free(r);
    return fails ? 1 : 0;
}
