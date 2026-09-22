/* robust -- position independence and corruption detection. */
#include "../src/imap_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

static uint8_t *slurp(const char *p, long *n) {
    FILE *f=fopen(p,"rb"); if(!f) return NULL;
    fseek(f,0,SEEK_END); *n=ftell(f); fseek(f,0,SEEK_SET);
    uint8_t *b=malloc((size_t)*n); if(fread(b,1,(size_t)*n,f)!=(size_t)*n){free(b);fclose(f);return NULL;}
    fclose(f); return b;
}
/* Recompute the header CRC in the trailer, so an edited header reaches the check it is
 * meant to exercise instead of being stopped by the checksum. */
static void reseal(uint8_t *b, long n) {
    uLong c = crc32(crc32(0L, Z_NULL, 0), b, 64);
    b[n-12] = (uint8_t)c; b[n-11] = (uint8_t)(c >> 8); b[n-10] = (uint8_t)(c >> 16); b[n-9] = (uint8_t)(c >> 24);
}
static int opens(const uint8_t *buf, long n) {
    imap_io *io = imap_io_open_mem(buf,n); if(!io) return 0;
    imap_file *f = imap_open_io(io,0);
    int ok = f!=NULL; if(f) imap_close(f); imap_io_close(io);
    return ok;
}

int main(int argc,char**argv){
    if(argc!=2){fprintf(stderr,"usage: robust file.imap\n");return 2;}
    long n; uint8_t *orig=slurp(argv[1],&n);
    if(!orig){perror("slurp");return 2;}
    int fails=0;
    #define CHECK(what,cond) do{ if(!(cond)){ printf("  FAIL %s\n",what); fails++; } \
                                 else printf("  ok   %s\n",what); }while(0)

    printf("baseline:\n");
    CHECK("opens from memory", opens(orig,n));

    printf("position independence (the container-embedding claim):\n");
    for (long pad=1; pad<=100000; pad*=10) {
        uint8_t *big=calloc((size_t)(n+pad+7),1);
        memcpy(big+pad,orig,(size_t)n);
        imap_io *outer=imap_io_open_mem(big,n+pad+7);
        imap_io *sl=imap_io_slice(outer,pad,n);
        imap_file *f= sl? imap_open_io(sl,0):NULL;
        uint32_t nc = f? imap_n_chunks(f):0;
        char lbl[64]; snprintf(lbl,sizeof lbl,"opens at offset %ld inside a larger buffer",pad);
        CHECK(lbl, f!=NULL && nc>0);
        if(f) imap_close(f);
        if(sl) imap_io_close(sl);
        imap_io_close(outer); free(big);
    }

    printf("corruption detection:\n");
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n); c[0]^=0xFF;
      CHECK("rejects bad header magic", !opens(c,n)); free(c); }
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n); c[n-1]^=0xFF;
      CHECK("rejects bad trailer magic", !opens(c,n)); free(c); }
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n); c[n-40]^=0xFF;   /* inside footer */
      CHECK("rejects footer bit flip (crc)", !opens(c,n)); free(c); }
    CHECK("rejects truncation (half)", !opens(orig,n/2));
    CHECK("rejects truncation (trailer only)", !opens(orig,32));
    CHECK("rejects empty", !opens(orig,0));
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n);
      c[16]^=0x01;                                   /* flip the order bit, no reseal */
      CHECK("rejects a flipped header feature bit (header crc)", !opens(c,n)); free(c); }
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n);
      c[24]^=0x20;                                   /* profile byte, no reseal */
      CHECK("rejects a header profile bit flip (header crc)", !opens(c,n)); free(c); }
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n);
      c[16]|=0x20; reseal(c,n);                      /* an unknown feature bit, sealed */
      CHECK("refuses unknown feature flag", !opens(c,n)); free(c); }
    { uint8_t *c=malloc((size_t)n); memcpy(c,orig,(size_t)n);
      c[8]=0x09; reseal(c,n);                        /* wrong major version, sealed */
      CHECK("refuses wrong major version", !opens(c,n)); free(c); }

    free(orig);
    printf("%s (%d failures)\n", fails?"FAILURES":"ALL PASS", fails);
    return fails?1:0;
}
