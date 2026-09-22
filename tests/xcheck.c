/* xcheck -- encode runs with the C codec, verify self round-trip, dump streams.
 *
 * in:  binary file of records {int64 a; int64 b; int64 len; int64 strand}
 * out: binary file  {uint64 n; then per stream: uint64 len, bytes}
 * Exits non-zero if the C codec cannot round-trip its own output. */
#include "../src/imap_codec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: xcheck runs.bin streams.bin\n"); return 2; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("open runs"); return 2; }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    long n = sz / 32;
    imap_run *runs = malloc((size_t)n * sizeof *runs);
    if (!runs) return 2;
    for (long i = 0; i < n; i++) {
        int64_t rec[4];
        if (fread(rec, sizeof rec, 1, f) != 1) { fprintf(stderr, "short read\n"); return 2; }
        runs[i].a = rec[0]; runs[i].b = rec[1];
        runs[i].len = rec[2]; runs[i].strand = (uint8_t)rec[3];
    }
    fclose(f);

    imap_buf st[IMAP_N_STREAMS];
    if (imap_encode_chunk(runs, n, st) != 0) { fprintf(stderr, "encode failed\n"); return 1; }

    imap_run *back = malloc((size_t)n * sizeof *back);
    int64_t b_enter0 = imap_b_enter(runs[0].b, runs[0].len, runs[0].strand);
    if (imap_decode_chunk(st, n, runs[0].a, b_enter0, back) != 0) {
        fprintf(stderr, "decode failed\n"); return 1;
    }
    for (long i = 0; i < n; i++) {
        if (back[i].a != runs[i].a || back[i].b != runs[i].b ||
            back[i].len != runs[i].len || back[i].strand != runs[i].strand) {
            fprintf(stderr, "C self round-trip MISMATCH at run %ld\n", i); return 1;
        }
    }

    FILE *o = fopen(argv[2], "wb");
    if (!o) { perror("open out"); return 2; }
    uint64_t nn = (uint64_t)n;
    fwrite(&nn, sizeof nn, 1, o);
    for (int i = 0; i < IMAP_N_STREAMS; i++) {
        uint64_t L = st[i].n;
        fwrite(&L, sizeof L, 1, o);
        if (L) fwrite(st[i].p, 1, L, o);
    }
    fclose(o);
    fprintf(stderr, "C: %ld runs, streams a=%zu b=%zu len=%zu strand=%zu, self round-trip OK\n",
            n, st[0].n, st[1].n, st[2].n, st[3].n);
    for (int i = 0; i < IMAP_N_STREAMS; i++) imap_buf_free(&st[i]);
    free(runs); free(back);
    return 0;
}
