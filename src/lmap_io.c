/* POSIX.1-2008 for pread/strdup/fstat, whatever -std= the embedding build uses. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "lmap_io.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <pthread.h>

struct lmap_io {
    const lmap_io_backend *be;
    void   *ctx;
    int64_t base;      /* offset of our region within the backend */
    int64_t len;       /* size of our region */
    lmap_io *parent;   /* non-NULL for slices; borrowed */
};

/* ------------------------------------------------------------------ file */

typedef struct { int fd; } file_ctx;

static int64_t file_read(void *vctx, void *buf, int64_t off, int64_t len) {
    file_ctx *c = (file_ctx *)vctx;
    int64_t got = 0;
    while (got < len) {
        ssize_t n = pread(c->fd, (char *)buf + got, (size_t)(len - got), (off_t)(off + got));
        if (n < 0) { if (errno == EINTR) continue; return -1; }
        if (n == 0) break;                    /* end of file */
        got += n;
    }
    return got;
}

static void file_close(void *vctx) {
    file_ctx *c = (file_ctx *)vctx;
    if (c->fd >= 0) close(c->fd);
    free(c);
}

static const lmap_io_backend FILE_BE = { "file", file_read, file_close };

lmap_io *lmap_io_open_file(const char *path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) return NULL;
    struct stat st;
    if (fstat(fd, &st) != 0) { close(fd); return NULL; }
    file_ctx *c = (file_ctx *)malloc(sizeof *c);
    lmap_io *io = (lmap_io *)calloc(1, sizeof *io);
    if (!c || !io) { free(c); free(io); close(fd); return NULL; }
    c->fd = fd;
    io->be = &FILE_BE; io->ctx = c; io->base = 0; io->len = (int64_t)st.st_size;
    return io;
}

/* ---------------------------------------------------------------- memory */

typedef struct { const uint8_t *p; int64_t n; } mem_ctx;

static int64_t mem_read(void *vctx, void *buf, int64_t off, int64_t len) {
    mem_ctx *c = (mem_ctx *)vctx;
    if (off < 0 || off > c->n) return -1;
    int64_t avail = c->n - off;
    int64_t n = len < avail ? len : avail;
    memcpy(buf, c->p + off, (size_t)n);
    return n;
}

static void mem_close(void *vctx) { free(vctx); }

static const lmap_io_backend MEM_BE = { "mem", mem_read, mem_close };

lmap_io *lmap_io_open_mem(const void *data, int64_t len) {
    mem_ctx *c = (mem_ctx *)malloc(sizeof *c);
    lmap_io *io = (lmap_io *)calloc(1, sizeof *io);
    if (!c || !io) { free(c); free(io); return NULL; }
    c->p = (const uint8_t *)data; c->n = len;
    io->be = &MEM_BE; io->ctx = c; io->base = 0; io->len = len;
    return io;
}

/* ----------------------------------------------------------------- slice */

lmap_io *lmap_io_slice(lmap_io *parent, int64_t base, int64_t len) {
    /* base+len would be signed overflow (UB) for large inputs; compare by subtraction. */
    if (!parent || base < 0 || len < 0) return NULL;
    if (base > parent->len || len > parent->len - base) return NULL;
    lmap_io *io = (lmap_io *)calloc(1, sizeof *io);
    if (!io) return NULL;
    io->be = parent->be; io->ctx = parent->ctx;
    io->base = parent->base + base; io->len = len;
    io->parent = parent;
    return io;
}

/* --------------------------------------------------------------- backend */

lmap_io *lmap_io_open_backend(const lmap_io_backend *be, void *ctx, int64_t size) {
    if (!be || !be->read || size < 0) return NULL;
    lmap_io *io = (lmap_io *)calloc(1, sizeof *io);
    if (!io) return NULL;
    io->be = be; io->ctx = ctx; io->base = 0; io->len = size;
    return io;
}

/* ----------------------------------------------------------------- cache */

typedef struct cblock {
    int64_t idx;                         /* block number; -1 if the slot is free */
    uint8_t *data; int64_t n;            /* n < block only for the last block */
    struct cblock *prev, *next;          /* LRU, head = most recent */
    struct cblock *hnext;
} cblock;

typedef struct {
    lmap_io *inner; int own;
    int64_t block, size;                 /* block bytes; inner size */
    size_t nslots;
    cblock *slots, *lru_head, *lru_tail, *free_list;
    cblock **htab; size_t hcap;
    int64_t last_miss; int64_t ra;       /* sequential detection, current read-ahead (blocks) */
    int64_t ra_max;
    uint8_t *scratch; int64_t scratch_n;
    pthread_mutex_t mu;
    lmap_io_cache_stats st;
} cache_ctx;

static size_t hslot(const cache_ctx *c, int64_t idx) { return (size_t)((uint64_t)idx * 0x9E3779B97F4A7C15ull >> 20) & (c->hcap - 1); }

static cblock *c_find(cache_ctx *c, int64_t idx) {
    for (cblock *b = c->htab[hslot(c, idx)]; b; b = b->hnext) if (b->idx == idx) return b;
    return NULL;
}
static void c_unlink(cache_ctx *c, cblock *b) {
    if (b->prev) b->prev->next = b->next; else c->lru_head = b->next;
    if (b->next) b->next->prev = b->prev; else c->lru_tail = b->prev;
    b->prev = b->next = NULL;
}
static void c_front(cache_ctx *c, cblock *b) {
    b->next = c->lru_head; b->prev = NULL;
    if (c->lru_head) c->lru_head->prev = b; else c->lru_tail = b;
    c->lru_head = b;
}
static void c_unhash(cache_ctx *c, cblock *b) {
    cblock **pp = &c->htab[hslot(c, b->idx)];
    while (*pp && *pp != b) pp = &(*pp)->hnext;
    if (*pp) *pp = b->hnext;
}
/* a slot for a new block: a free one, else the least recently used (never `keep`) */
static cblock *c_slot(cache_ctx *c) {
    cblock *b = c->free_list;
    if (b) { c->free_list = b->hnext; return b; }
    b = c->lru_tail;
    c_unlink(c, b); c_unhash(c, b);
    return b;
}
static void c_install(cache_ctx *c, int64_t idx, const uint8_t *src, int64_t n) {
    cblock *b = c_slot(c);
    b->idx = idx; b->n = n;
    memcpy(b->data, src, (size_t)n);
    b->hnext = c->htab[hslot(c, idx)]; c->htab[hslot(c, idx)] = b;
    c_front(c, b);
}

/* Fetch blocks [first, first+count) from the inner handle in one read and cache them. */
static int c_fetch(cache_ctx *c, int64_t first, int64_t count) {
    int64_t off = first * c->block;
    int64_t end = (first + count) * c->block;
    if (end > c->size) end = c->size;
    int64_t n = end - off;
    if (n > c->scratch_n) {
        uint8_t *ns = realloc(c->scratch, (size_t)n);
        if (!ns) return -1;
        c->scratch = ns; c->scratch_n = n;
    }
    if (lmap_pread(c->inner, c->scratch, off, n) != 0) return -1;
    c->st.requests++; c->st.bytes_fetched += (uint64_t)n;
    for (int64_t k = 0; k < count && k * c->block < n; k++) {
        int64_t bn = n - k * c->block < c->block ? n - k * c->block : c->block;
        if (!c_find(c, first + k)) c_install(c, first + k, c->scratch + k * c->block, bn);
    }
    return 0;
}

static int64_t cache_read(void *vctx, void *buf, int64_t off, int64_t len) {
    cache_ctx *c = (cache_ctx *)vctx;
    if (off >= c->size) return 0;
    if (len > c->size - off) len = c->size - off;
    /* big reads bypass the cache: they would only evict what is worth keeping */
    if ((uint64_t)len > (uint64_t)c->nslots * (uint64_t)c->block / 4) {
        if (lmap_pread(c->inner, buf, off, len) != 0) return -1;
        pthread_mutex_lock(&c->mu);
        c->st.requests++; c->st.bytes_fetched += (uint64_t)len;
        pthread_mutex_unlock(&c->mu);
        return len;
    }
    pthread_mutex_lock(&c->mu);
    int64_t first = off / c->block, last = (off + len - 1) / c->block;
    int64_t got = 0;
    int rc = 0;
    for (int64_t i = first; i <= last && rc == 0; ) {
        cblock *b = c_find(c, i);
        if (!b) {
            /* a run of missing blocks, fetched together, plus read-ahead if this miss
             * continues the last one */
            int64_t j = i;
            while (j <= last && !c_find(c, j)) j++;
            int64_t count = j - i;
            if (i == c->last_miss + 1) {
                c->ra = c->ra ? c->ra * 2 : 1;
                if (c->ra * c->block > c->ra_max) c->ra = c->ra_max / c->block;
            } else c->ra = 0;
            int64_t nblocks = (c->size + c->block - 1) / c->block;
            int64_t want = count + c->ra;
            if (want > (int64_t)c->nslots / 2) want = (int64_t)c->nslots / 2 > count ? (int64_t)c->nslots / 2 : count;
            if (i + want > nblocks) want = nblocks - i;
            c->st.misses += (uint64_t)count;
            if (c_fetch(c, i, want) != 0) { rc = -1; break; }
            c->last_miss = i + want - 1;
            b = c_find(c, i);
            if (!b) { rc = -1; break; }
        } else {
            c->st.hits++;
            c_unlink(c, b); c_front(c, b);
        }
        int64_t bstart = i * c->block;
        int64_t from = off + got - bstart;            /* offset within the block */
        int64_t take = b->n - from;
        if (take > len - got) take = len - got;
        if (take <= 0) { rc = -1; break; }
        memcpy((uint8_t *)buf + got, b->data + from, (size_t)take);
        got += take;
        i++;
    }
    pthread_mutex_unlock(&c->mu);
    return rc == 0 ? got : -1;
}

static void cache_close(void *vctx) {
    cache_ctx *c = (cache_ctx *)vctx;
    if (c->own) lmap_io_close(c->inner);
    for (size_t k = 0; k < c->nslots; k++) free(c->slots[k].data);
    free(c->slots); free(c->htab); free(c->scratch);
    pthread_mutex_destroy(&c->mu);
    free(c);
}

static const lmap_io_backend CACHE_BE = { "cache", cache_read, cache_close };

lmap_io *lmap_io_cache(lmap_io *inner, size_t block_bytes, size_t capacity_bytes,
                       size_t readahead_max_bytes, int own_inner) {
    if (!inner) return NULL;
    if (!block_bytes) block_bytes = (size_t)64 << 10;
    if (!capacity_bytes) capacity_bytes = (size_t)64 << 20;
    if (!readahead_max_bytes) readahead_max_bytes = (size_t)4 << 20;
    size_t nslots = capacity_bytes / block_bytes;
    if (nslots < 4) nslots = 4;
    cache_ctx *c = calloc(1, sizeof *c);
    if (!c) return NULL;
    c->inner = inner; c->own = own_inner;
    c->block = (int64_t)block_bytes; c->size = lmap_io_size(inner);
    c->nslots = nslots; c->last_miss = -2;
    c->ra_max = (int64_t)(readahead_max_bytes > block_bytes ? readahead_max_bytes : block_bytes);
    c->hcap = 1; while (c->hcap < nslots * 2) c->hcap <<= 1;
    c->slots = calloc(nslots, sizeof *c->slots);
    c->htab = calloc(c->hcap, sizeof *c->htab);
    if (!c->slots || !c->htab || pthread_mutex_init(&c->mu, NULL) != 0) {
        free(c->slots); free(c->htab); free(c); return NULL;
    }
    for (size_t k = 0; k < nslots; k++) {
        if (!(c->slots[k].data = malloc(block_bytes))) {
            for (size_t j = 0; j < k; j++) free(c->slots[j].data);
            free(c->slots); free(c->htab); pthread_mutex_destroy(&c->mu); free(c); return NULL;
        }
        c->slots[k].idx = -1;
        c->slots[k].hnext = c->free_list; c->free_list = &c->slots[k];
    }
    lmap_io *io = lmap_io_open_backend(&CACHE_BE, c, c->size);
    if (!io) { cache_close(c); return NULL; }
    return io;
}

int lmap_io_get_cache_stats(const lmap_io *io, lmap_io_cache_stats *st) {
    if (!io || !st || io->be != &CACHE_BE) return -1;
    cache_ctx *c = (cache_ctx *)io->ctx;
    pthread_mutex_lock(&c->mu);
    *st = c->st;
    pthread_mutex_unlock(&c->mu);
    return 0;
}

/* ------------------------------------------------------------------- api */

int lmap_pread(lmap_io *io, void *buf, int64_t off, int64_t len) {
    if (!io || off < 0 || len < 0) return -1;
    if (off > io->len || len > io->len - off) return -1;   /* no overflow: both >= 0 */
    if (len == 0) return 0;
    int64_t n = io->be->read(io->ctx, buf, io->base + off, len);
    return (n == len) ? 0 : -1;
}

int64_t lmap_io_size(const lmap_io *io) { return io ? io->len : -1; }

void lmap_io_close(lmap_io *io) {
    if (!io) return;
    if (!io->parent && io->be->close) io->be->close(io->ctx);   /* slices borrow */
    free(io);
}
