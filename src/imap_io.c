#include "imap_io.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

struct imap_io {
    const imap_io_backend *be;
    void   *ctx;
    int64_t base;      /* offset of our region within the backend */
    int64_t len;       /* size of our region */
    imap_io *parent;   /* non-NULL for slices; borrowed */
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

static const imap_io_backend FILE_BE = { "file", file_read, file_close };

imap_io *imap_io_open_file(const char *path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) return NULL;
    struct stat st;
    if (fstat(fd, &st) != 0) { close(fd); return NULL; }
    file_ctx *c = (file_ctx *)malloc(sizeof *c);
    imap_io *io = (imap_io *)calloc(1, sizeof *io);
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

static const imap_io_backend MEM_BE = { "mem", mem_read, mem_close };

imap_io *imap_io_open_mem(const void *data, int64_t len) {
    mem_ctx *c = (mem_ctx *)malloc(sizeof *c);
    imap_io *io = (imap_io *)calloc(1, sizeof *io);
    if (!c || !io) { free(c); free(io); return NULL; }
    c->p = (const uint8_t *)data; c->n = len;
    io->be = &MEM_BE; io->ctx = c; io->base = 0; io->len = len;
    return io;
}

/* ----------------------------------------------------------------- slice */

imap_io *imap_io_slice(imap_io *parent, int64_t base, int64_t len) {
    if (!parent || base < 0 || len < 0 || base + len > parent->len) return NULL;
    imap_io *io = (imap_io *)calloc(1, sizeof *io);
    if (!io) return NULL;
    io->be = parent->be; io->ctx = parent->ctx;
    io->base = parent->base + base; io->len = len;
    io->parent = parent;
    return io;
}

/* ------------------------------------------------------------------- api */

int imap_pread(imap_io *io, void *buf, int64_t off, int64_t len) {
    if (!io || off < 0 || len < 0) return -1;
    if (off > io->len || len > io->len - off) return -1;   /* no overflow: both >= 0 */
    if (len == 0) return 0;
    int64_t n = io->be->read(io->ctx, buf, io->base + off, len);
    return (n == len) ? 0 : -1;
}

int64_t imap_io_size(const imap_io *io) { return io ? io->len : -1; }

void imap_io_close(imap_io *io) {
    if (!io) return;
    if (!io->parent && io->be->close) io->be->close(io->ctx);   /* slices borrow */
    free(io);
}
