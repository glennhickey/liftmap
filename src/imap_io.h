/* imap_io -- the single I/O seam.  See doc/SPEC.md "must not preclude" #2.
 *
 * Every read in the library goes through imap_pread() on a handle that carries a
 * base offset and a length.  That is what lets one file live standalone, inside a
 * container at some offset, behind an HTTP range fetcher, or in memory, without the
 * rest of the library knowing which.  Offsets passed in are always relative to the
 * mapped region's own byte 0.
 */
#ifndef IMAP_IO_H
#define IMAP_IO_H

#include <stddef.h>
#include <stdint.h>

typedef struct imap_io imap_io;

/* Backend vtable.  read() returns bytes read, or -1; it must fill the whole
 * request or return short only at end of region. */
typedef struct {
    const char *name;
    int64_t (*read)(void *ctx, void *buf, int64_t off, int64_t len);
    void    (*close)(void *ctx);
} imap_io_backend;

/* Open a plain local file read-only. Returns NULL on failure. */
imap_io *imap_io_open_file(const char *path);

/* Wrap a caller-owned memory buffer (not copied; must outlive the handle). */
imap_io *imap_io_open_mem(const void *data, int64_t len);

/* Restrict an existing handle to [base, base+len) -- container slice.  The
 * returned handle borrows the parent, which must outlive it. */
imap_io *imap_io_slice(imap_io *parent, int64_t base, int64_t len);

/* Read exactly len bytes at off (relative to this handle's region).
 * Returns 0 on success, -1 on short read or error. */
int imap_pread(imap_io *io, void *buf, int64_t off, int64_t len);

int64_t imap_io_size(const imap_io *io);
void    imap_io_close(imap_io *io);

#endif /* IMAP_IO_H */
