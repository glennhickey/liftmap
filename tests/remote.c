/* remote -- reading through a custom backend and the block cache, as a remote file is read.
 *
 * Opens one liftmap file twice: directly, and through a counting backend (built with
 * lmap_io_open_backend, the way an application plugs in HTTP) behind lmap_io_cache.  Runs
 * the same queries, cursor scans and point sweeps on both and requires identical results,
 * then reports how many backend reads the cache needed -- the number of HTTP range
 * requests a remote reader would make.  Optional per-read latency (LMAP_TEST_LATENCY_US)
 * makes the effect visible in wall time.
 *
 * usage: remote file.lmap [nqueries]
 */
#define _POSIX_C_SOURCE 200809L
#include "../src/liftmap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { lmap_io *inner; long reads; long long bytes; long latency_us; } counting;

static int64_t counting_read(void *ctx, void *buf, int64_t off, int64_t len) {
    counting *c = ctx;
    __atomic_add_fetch(&c->reads, 1, __ATOMIC_RELAXED);
    __atomic_add_fetch(&c->bytes, len, __ATOMIC_RELAXED);
    if (c->latency_us) { struct timespec t = { 0, c->latency_us * 1000 }; nanosleep(&t, NULL); }
    return lmap_pread(c->inner, buf, off, len) == 0 ? len : -1;
}
static void counting_close(void *ctx) { counting *c = ctx; lmap_io_close(c->inner); free(c); }
static const lmap_io_backend COUNTING = { "counting", counting_read, counting_close };

static uint64_t rng = 0x243F6A8885A308D3ull;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }

/* One workload, as a checksum of every hit, so two readers can be compared. */
static uint64_t workload(lmap_file *f, long nq, uint64_t seed) {
    uint64_t save = rng, h = 1469598103934665603ull;
    rng = seed;
#define MIX(v) (h = (h ^ (uint64_t)(v)) * 1099511628211ull)
    for (long q = 0; q < nq; q++) {
        int ax = (int)(xr() & 1);
        uint32_t nm = lmap_n_members(f, ax);
        if (!nm) continue;
        uint32_t m = (uint32_t)(xr() % nm);
        const lmap_member *mm = lmap_member_at(f, ax, m);
        int64_t len = (int64_t)mm->length ? (int64_t)mm->length : 1;
        int64_t lo = (int64_t)(xr() % (uint64_t)len), w = (int64_t[]){ 10, 1000, 100000 }[xr() % 3];
        lmap_hit *hit; size_t n;
        int rc = ax ? lmap_query_b(f, m, lo, lo + w, &hit, &n, NULL) : lmap_query_a(f, m, lo, lo + w, &hit, &n, NULL);
        MIX(rc);
        if (rc == 0) { for (size_t k = 0; k < n; k++) { MIX(hit[k].a); MIX(hit[k].b); MIX(hit[k].len); MIX(hit[k].strand); } free(hit); }
        if (q % 10 == 0) {                       /* a point sweep through a cached cursor */
            lmap_cursor *c = lmap_cursor_open(f, ax, m, NULL, 0, (size_t)8 << 20);
            lmap_hit ph[16];
            for (int64_t p = lo; c && p < lo + 2000 && p < len; p += 7) {
                int k = lmap_cursor_point(c, p, ph, 16);
                MIX(k);
                for (int j = 0; j < k && j < 16; j++) { MIX(ph[j].a); MIX(ph[j].b); }
            }
            lmap_cursor_close(c);
        }
    }
    rng = save;
    return h;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: remote file.lmap [nqueries]\n"); return 2; }
    long nq = argc > 2 ? atol(argv[2]) : 500;
    const char *lat = getenv("LMAP_TEST_LATENCY_US");
    int fails = 0;
#define CHECK(what, cond) do { if (!(cond)) { printf("  FAIL %s\n", what); fails++; } else printf("  ok   %s\n", what); } while (0)

    double t0 = now();
    lmap_file *direct = lmap_open(argv[1]);
    if (!direct) { printf("cannot open %s\n", argv[1]); return 1; }
    uint64_t want = workload(direct, nq, 42);
    double t_direct = now() - t0;

    /* uncached: every read is a request */
    counting *u = calloc(1, sizeof *u);
    u->inner = lmap_io_open_file(argv[1]); u->latency_us = lat ? atol(lat) : 0;
    lmap_io *uio = lmap_io_open_backend(&COUNTING, u, lmap_io_size(u->inner));
    t0 = now();
    lmap_file *fu = lmap_open_io(uio, 1);
    CHECK("opens through a custom backend", fu != NULL);
    long open_reads_u = u->reads;
    uint64_t got_u = fu ? workload(fu, nq, 42) : 0;
    double t_u = now() - t0;
    long reads_u = u->reads; long long bytes_u = u->bytes;
    CHECK("same results through the backend", got_u == want);

    /* cached */
    counting *c = calloc(1, sizeof *c);
    c->inner = lmap_io_open_file(argv[1]); c->latency_us = u->latency_us;
    const char *eb = getenv("LMAP_TEST_BLOCK"), *er = getenv("LMAP_TEST_READAHEAD");
    lmap_io *cio = lmap_io_cache(lmap_io_open_backend(&COUNTING, c, lmap_io_size(c->inner)),
                                 eb ? (size_t)atol(eb) : 0, 0, er ? (size_t)atol(er) : 0, 1);
    t0 = now();
    lmap_file *fc = lmap_open_io(cio, 0);
    CHECK("opens through the cache", fc != NULL);
    long open_reads_c = c->reads;
    uint64_t got_c = fc ? workload(fc, nq, 42) : 0;
    double t_c = now() - t0;
    CHECK("same results through the cache", got_c == want);
    lmap_io_cache_stats cs;
    CHECK("cache stats are available", lmap_io_get_cache_stats(cio, &cs) == 0 && (long)cs.requests == c->reads);
    CHECK("the cache needs fewer backend reads", c->reads < reads_u);
    {   /* a small workload, then the same again on a fresh cache: the repeat needs no reads */
        counting *r = calloc(1, sizeof *r);
        r->inner = lmap_io_open_file(argv[1]);
        lmap_io *rio = lmap_io_cache(lmap_io_open_backend(&COUNTING, r, lmap_io_size(r->inner)), 0, 0, 0, 1);
        lmap_file *fr = lmap_open_io(rio, 0);
        uint64_t h1 = fr ? workload(fr, 20, 7) : 0;
        long before = r->reads;
        uint64_t h2 = fr ? workload(fr, 20, 7) : 1;
        CHECK("repeating a small workload is served from the cache (0 new reads)", h1 == h2 && r->reads == before);
        if (fr) lmap_close(fr);
        lmap_io_close(rio);
    }

    printf("open:     %ld backend reads uncached, %ld cached\n", open_reads_u, open_reads_c);
    printf("workload: %ld queries; backend reads %ld uncached (%.1f MB) vs %ld cached (%.1f MB, %llu hits / %llu misses)\n",
           nq, reads_u, bytes_u / 1e6, (long)cs.requests, cs.bytes_fetched / 1e6,
           (unsigned long long)cs.hits, (unsigned long long)cs.misses);
    printf("time:     direct %.2f s, uncached %.2f s, cached %.2f s%s\n", t_direct, t_u, t_c,
           lat ? " (with simulated latency)" : "");
    if (fu) lmap_close(fu);
    if (fc) lmap_close(fc);
    lmap_io_close(cio);
    lmap_close(direct);
    printf("%s (%d failures)\n", fails ? "FAILURES" : "ALL PASS", fails);
    return fails ? 1 : 0;
}
