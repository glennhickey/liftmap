/* query -- check imap_query_a / imap_query_b against a brute-force oracle.
 *
 * in: binary {int64 a; int64 b; int64 len; int64 strand; int32 amem; int32 bmem} x N,
 *     grouped by amem with a increasing within a member.
 * Axis-a member ids are remapped (reversed) before writing, so the input no longer
 * arrives in id order and the writer's dir.a sort is exercised.
 * The oracle maps a window to the other axis via the explicit base pairing
 * pair(i) = strand ? other+len-1-i : other+i, not via the library's clip formula. */
#include "../src/imap_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int64_t a,b,len,strand; int32_t am,bm; } inrec;

static int cmp_a(const void *x, const void *y) {
    const imap_hit *p=x,*q=y; return (p->a>q->a)-(p->a<q->a);
}
static int cmp_b(const void *x, const void *y) {
    const imap_hit *p=x,*q=y;
    if (p->b!=q->b) return (p->b>q->b)-(p->b<q->b);
    if (p->a_member!=q->a_member) return (p->a_member>q->a_member)-(p->a_member<q->a_member);
    return (p->a>q->a)-(p->a<q->a);
}

/* oracle clip, derived independently of the library */
static int oracle(const inrec *r, int axis, int64_t lo, int64_t hi, imap_hit *h) {
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
    if (argc != 4) { fprintf(stderr,"usage: query runs.bin out.imap nqueries\n"); return 2; }
    FILE *fp = fopen(argv[1],"rb"); if (!fp) { perror("runs"); return 2; }
    fseek(fp,0,SEEK_END); long sz=ftell(fp); fseek(fp,0,SEEK_SET);
    long n = sz/(long)sizeof(inrec);
    inrec *r = malloc((size_t)n*sizeof *r);
    if (fread(r,sizeof *r,(size_t)n,fp)!=(size_t)n) return 2;
    fclose(fp);
    int32_t na=0, nb=0;
    for (long i=0;i<n;i++){ if(r[i].am+1>na) na=r[i].am+1; if(r[i].bm+1>nb) nb=r[i].bm+1; }
    for (long i=0;i<n;i++) r[i].am = na-1-r[i].am;          /* remap: ids no longer in input order */

    imap_writer *w = imap_writer_open(argv[2],"hal2.edge","imap 1\n");
    char nm[64];
    for (int32_t i=0;i<na;i++){ snprintf(nm,sizeof nm,"a%d",i); imap_writer_add_member(w,0,nm,1ULL<<40); }
    for (int32_t i=0;i<nb;i++){ snprintf(nm,sizeof nm,"b%d",i); imap_writer_add_member(w,1,nm,1ULL<<40); }
    for (long i=0;i<n;i++)
        if (imap_writer_add_run(w,(uint32_t)r[i].am,(uint32_t)r[i].bm,r[i].a,r[i].b,r[i].len,(uint8_t)r[i].strand)) {
            fprintf(stderr,"add_run %ld failed\n",i); return 1; }
    if (imap_writer_close(w)) { fprintf(stderr,"close failed\n"); return 1; }
    imap_file *f = imap_open(argv[2]);
    if (!f) { fprintf(stderr,"open failed (invariant check?)\n"); return 1; }

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
    imap_hit *exp = malloc((size_t)expcap * sizeof *exp);

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

        imap_hit *got; size_t ng; imap_query_stats st;
        int rc = axis ? imap_query_b(f,m,lo,hi,&got,&ng,&st) : imap_query_a(f,m,lo,hi,&got,&ng,&st);
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
    printf("%ld queries, %ld hits compared, %ld mismatches\n", nq, total_hits, fails);
    printf("mean chunks decoded per query, by half-width:\n  %-9s", "");
    for (int k=0;k<NW;k++) printf("%10lld", (long long)widths[k]);
    for (int ax=0; ax<2; ax++) {
        printf("\n  axis %c   ", ax?'b':'a');
        for (int k=0;k<NW;k++) printf("%10.2f", dec_n[ax][k] ? dec_sum[ax][k]/dec_n[ax][k] : 0);
    }
    printf("\n");
    imap_close(f);
    free(exp); free(ia); free(ib); free(fa); free(fb); free(ca); free(cb); free(r);
    return fails ? 1 : 0;
}
