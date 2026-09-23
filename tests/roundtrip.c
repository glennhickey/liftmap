/* roundtrip -- write runs through the container, read them back, compare.
 *
 * in: binary {int64 a; int64 b; int64 len; int64 strand; int32 amem; int32 bmem} x N
 * Exits non-zero on any mismatch. */
#include "../src/lmap_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int64_t a,b,len,strand; int32_t am,bm; } inrec;

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr,"usage: roundtrip runs.bin out.lmap count\n"); return 2; }
    FILE *f = fopen(argv[1],"rb");
    if (!f) { perror("runs"); return 2; }
    fseek(f,0,SEEK_END); long sz=ftell(f); fseek(f,0,SEEK_SET);
    long n = sz / (long)sizeof(inrec);
    inrec *r = malloc((size_t)n * sizeof *r);
    if (fread(r,sizeof *r,(size_t)n,f)!=(size_t)n) { fprintf(stderr,"short read\n"); return 2; }
    fclose(f);

    int32_t namem=0, nbmem=0;
    for (long i=0;i<n;i++){ if(r[i].am+1>namem) namem=r[i].am+1; if(r[i].bm+1>nbmem) nbmem=r[i].bm+1; }

    lmap_writer *w = lmap_writer_open(argv[2], "hal2.edge");
    if (!w) { fprintf(stderr,"writer open failed\n"); return 1; }
    if (lmap_writer_set_params(w,(uint32_t)atoi(argv[3]),1000000,LMAP_CODEC_DEFLATE)!=0) return 1;
    char nm[64];
    for (int32_t i=0;i<namem;i++){ snprintf(nm,sizeof nm,"achr%d",i); if(lmap_writer_add_member(w,0,nm,1ULL<<40)<0) return 1; }
    for (int32_t i=0;i<nbmem;i++){ snprintf(nm,sizeof nm,"bscaf%d",i); if(lmap_writer_add_member(w,1,nm,1ULL<<40)<0) return 1; }
    for (long i=0;i<n;i++)
        if (lmap_writer_add_run(w,(uint32_t)r[i].am,(uint32_t)r[i].bm,r[i].a,r[i].b,r[i].len,(uint8_t)r[i].strand)!=0) {
            fprintf(stderr,"add_run failed at %ld\n",i); return 1; }
    if (lmap_writer_close(w)!=0) { fprintf(stderr,"writer close failed\n"); return 1; }

    lmap_file *m = lmap_open(argv[2]);
    if (!m) { fprintf(stderr,"open failed\n"); return 1; }
    if (strcmp(lmap_profile(m),"hal2.edge")) { fprintf(stderr,"profile mismatch\n"); return 1; }

    long total=0; uint32_t nc=lmap_n_chunks(m);
    lmap_run *buf=NULL; uint32_t bufn=0;
    for (uint32_t ci=0; ci<nc; ci++) {
        const lmap_chunk *c = lmap_chunk_at(m,ci);
        if (c->n_runs>bufn) { buf=realloc(buf,c->n_runs*sizeof *buf); bufn=c->n_runs; }
        if (lmap_read_chunk(m,ci,buf)!=0) { fprintf(stderr,"read_chunk %u failed\n",ci); return 1; }
        for (uint32_t k=0;k<c->n_runs;k++,total++) {
            inrec *o=&r[total];
            if (buf[k].a!=o->a||buf[k].b!=o->b||buf[k].len!=o->len||buf[k].strand!=(uint8_t)o->strand
                || c->a_member!=(uint32_t)o->am || c->b_member!=(uint32_t)o->bm) {
                fprintf(stderr,"MISMATCH chunk %u run %u (global %ld)\n",ci,k,total); return 1;
            }
        }
    }
    if (total!=n) { fprintf(stderr,"run count %ld != %ld\n",total,n); return 1; }

    /* member lookup by name must round-trip */
    for (int32_t i=0;i<nbmem;i++){ snprintf(nm,sizeof nm,"bscaf%d",i);
        if (lmap_member_by_name(m,1,nm)!=i){ fprintf(stderr,"name lookup failed %s\n",nm); return 1; } }
    if (lmap_member_by_name(m,1,"nope")!=-1){ fprintf(stderr,"absent name found\n"); return 1; }

    fprintf(stderr,"OK %ld runs, %u chunks, %u a-members, %u b-members\n",
            n,nc,lmap_n_members(m,0),lmap_n_members(m,1));
    free(buf); lmap_close(m); free(r);
    return 0;
}
