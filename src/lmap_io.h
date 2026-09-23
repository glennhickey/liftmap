/* lmap_io -- the single I/O seam.  See doc/SPEC.md "must not preclude" #2.
 *
 * Every read in the library goes through lmap_pread() on a handle that carries a
 * base offset and a length.  That is what lets one file live standalone, inside a
 * container at some offset, behind an HTTP range fetcher, or in memory, without the
 * rest of the library knowing which.  Offsets passed in are always relative to the
 * mapped region's own byte 0.
 */
#ifndef LMAP_IO_H
#define LMAP_IO_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct lmap_io lmap_io;

/* Backend vtable.  read() returns bytes read, or -1; it must fill the whole
 * request or return short only at end of region. */
typedef struct {
    const char *name;
    int64_t (*read)(void *ctx, void *buf, int64_t off, int64_t len);
    void    (*close)(void *ctx);
} lmap_io_backend;

/* Open a plain local file read-only. Returns NULL on failure. */
lmap_io *lmap_io_open_file(const char *path);

/* Wrap a caller-owned memory buffer (not copied; must outlive the handle). */
lmap_io *lmap_io_open_mem(const void *data, int64_t len);

/* Restrict an existing handle to [base, base+len) -- container slice.  The
 * returned handle borrows the parent, which must outlive it. */
lmap_io *lmap_io_slice(lmap_io *parent, int64_t base, int64_t len);

/* Any source: a backend (read may be called from several threads at once unless it is
 * wrapped in lmap_io_cache, which serialises its own reads), its context and the size
 * of what it reads.  The handle owns ctx: lmap_io_close calls be->close(ctx).  This is
 * how an application plugs in its own remote access -- htslib's hFILE, UCSC's udc, a
 * cloud SDK -- without liftmap depending on any of them. */
lmap_io *lmap_io_open_backend(const lmap_io_backend *be, void *ctx, int64_t size);

/* A block cache in front of any handle, for sources where each read is a round trip
 * (HTTP range requests).  Reads are served from aligned blocks of block_bytes (0: 64 KiB),
 * at most capacity_bytes of them (0: 64 MiB), least recently used evicted first.  Missing
 * blocks adjacent to each other are fetched in one inner read, and a miss that follows
 * the previous one sequentially reads ahead, doubling up to readahead_max_bytes (0: 4 MiB)
 * -- so a scan becomes a few large requests.  A read larger than a quarter of the cache
 * goes straight through.  Thread-safe.  With own_inner the cache closes `inner`. */
lmap_io *lmap_io_cache(lmap_io *inner, size_t block_bytes, size_t capacity_bytes,
                       size_t readahead_max_bytes, int own_inner);

typedef struct {
    uint64_t requests;         /* reads issued to the inner handle */
    uint64_t bytes_fetched;    /* bytes they returned */
    uint64_t hits, misses;     /* blocks served from the cache / fetched */
} lmap_io_cache_stats;
/* 0 and *st filled if io is a cache (or a slice of one); -1 otherwise. */
int lmap_io_get_cache_stats(const lmap_io *io, lmap_io_cache_stats *st);

/* HTTP(S) through libcurl range requests -- only in a build with HTTP=1 (-DLMAP_HTTP,
 * -lcurl); otherwise NULL.  lmap_http_available says which.  Wrap it in lmap_io_cache;
 * lmap_open does both for a path containing "://". */
lmap_io *lmap_io_open_url(const char *url);
int      lmap_http_available(void);

/* Read exactly len bytes at off (relative to this handle's region).
 * Returns 0 on success, -1 on short read or error. */
int lmap_pread(lmap_io *io, void *buf, int64_t off, int64_t len);

int64_t lmap_io_size(const lmap_io *io);
void    lmap_io_close(lmap_io *io);

#ifdef __cplusplus
}
#endif

#endif /* LMAP_IO_H */
