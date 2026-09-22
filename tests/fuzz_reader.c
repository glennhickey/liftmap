/* fuzz_reader -- open one .lmap and exercise every read path: all chunks (up to a
 * cap), and axis-a and axis-b queries of several widths on several members.  Build
 * with -fsanitize=address,undefined and drive it with tests/crcfuzz.py.
 * Exit 1 if any decoded run or query hit violates the format's invariants. */
#include "../src/lmap_file.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: fuzz_reader file.lmap\n"); return 2; }
    lmap_file *f = lmap_open(argv[1]);
    if (!f) { printf("  REJECTED at open\n"); return 0; }
    uint32_t nc = lmap_n_chunks(f);
    long bad = 0, decoded = 0, hits = 0;
    for (uint32_t i = 0; i < nc && i < 200; i++) {
        const lmap_chunk *c = lmap_chunk_at(f, i);
        if (!c || c->n_runs == 0 || c->n_runs > 5000000u) continue;
        lmap_run *r = malloc((size_t)c->n_runs * sizeof *r);
        if (r && lmap_read_chunk(f, i, r) == 0) {
            decoded++;
            for (uint32_t k = 0; k < c->n_runs; k++)
                if (r[k].len < 1 || r[k].a < 0 || r[k].b < 0) bad++;
        }
        free(r);
    }
    const int64_t widths[] = { 1, 1000, 1000000, 100000000 };
    for (int ax = 0; ax < 2; ax++) {
        uint32_t nm = lmap_n_members(f, ax);
        for (uint32_t m = 0; m < nm && m < 8; m++) {
            const lmap_member *mem = lmap_member_at(f, ax, m);
            int64_t mid = (int64_t)((mem->axis_min + mem->axis_max) / 2);
            for (int w = 0; w < 4; w++) {
                lmap_hit *h; size_t n;
                int64_t lo = mid - widths[w] < 0 ? 0 : mid - widths[w];
                int rc = ax ? lmap_query_b(f, m, lo, mid + widths[w], &h, &n, NULL)
                            : lmap_query_a(f, m, lo, mid + widths[w], &h, &n, NULL);
                if (rc != 0) continue;
                for (size_t k = 0; k < n; k++) {
                    if (h[k].len < 1 || h[k].a < 0 || h[k].b < 0) bad++;
                    int64_t s = ax ? h[k].b : h[k].a;
                    if (s < lo || s + h[k].len > mid + widths[w]) bad++;   /* clipped? */
                }
                hits += (long)n;
                free(h);
            }
        }
    }
    printf("  opened: %u chunks, %u/%u members, %ld chunks decoded, %ld hits, %ld invalid\n",
           nc, lmap_n_members(f, 0), lmap_n_members(f, 1), decoded, hits, bad);
    lmap_close(f);
    return bad ? 1 : 0;
}
