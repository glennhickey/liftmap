/* imap_file -- container framing for the interval map format. See doc/SPEC.md S2.
 *
 * Layout:   [64 B header] [sections...] [footer] [32 B trailer]
 * Open is two reads: the fixed trailer, then the footer it points at.
 * Every offset stored in the file is relative to the file's own byte 0.
 */
#ifndef IMAP_FILE_H
#define IMAP_FILE_H

#include "imap_codec.h"
#include "imap_io.h"

#define IMAP_MAGIC        "IMAP\x1A\x0D\x0A"   /* 7 chars + implicit NUL = 8 bytes */
#define IMAP_HEADER_SIZE  64
#define IMAP_TRAILER_SIZE 32
#define IMAP_DIRA_ENTRY   64
#define IMAP_DIRB_ENTRY   12
#define IMAP_MEM_ENTRY    40

/* Format limits (SPEC 2.2).  A chunk's decode buffers are sized from its directory
 * entry, so these, together with the byte-count checks at open, are what stop a few
 * hundred bytes of hostile file from demanding gigabytes. */
#define IMAP_MAX_CHUNK_RUNS  (1u << 24)
#define IMAP_MAX_CHUNK_RAW   (1u << 30)
#define IMAP_ZLIB_MAX_RATIO  1032u      /* deflate's maximum expansion */

#define IMAP_CODEC_NONE 0
#define IMAP_CODEC_ZLIB 1

/* One chunk's directory entry, in memory (native types; serialized explicitly). */
typedef struct {
    uint32_t a_member, b_member;
    uint64_t a_min;  uint32_t a_span;
    uint64_t b_min;  uint32_t b_span;
    uint64_t b_enter0;          /* first run's axis-b entry point -- a decode base */
    uint64_t off;               /* within the runs section */
    uint32_t clen, rawlen, n_runs;
    uint32_t codec;            /* IMAP_CODEC_* for this chunk's blob */
} imap_chunk;

/* One member (sequence) of an axis. */
typedef struct {
    char    *name;
    uint64_t length;
    uint64_t axis_min, axis_max;
    uint32_t first_chunk, n_chunks;
} imap_member;

/* ------------------------------------------------------------------ writer */

typedef struct imap_writer imap_writer;

/* `profile` is e.g. "hal2.edge"; `schema` is the schema text, stored verbatim. */
imap_writer *imap_writer_open(const char *path, const char *profile, const char *schema);

/* Members must be declared before any run referencing them. Returns member id. */
int32_t imap_writer_add_member(imap_writer *w, int axis, const char *name, uint64_t length);

/* Runs must arrive in axis-a order within a member. The writer closes chunks per
 * SPEC S1.3 (count / bspan / member change) and never lets a chunk span members. */
int imap_writer_add_run(imap_writer *w, uint32_t a_member, uint32_t b_member,
                        int64_t a, int64_t b, int64_t len, uint8_t strand);

int imap_writer_set_params(imap_writer *w, uint32_t count, uint64_t bspan, int codec);
int imap_writer_close(imap_writer *w);          /* writes footer+trailer, frees */
void imap_writer_abort(imap_writer *w);

/* ------------------------------------------------------------------ reader */

typedef struct imap_file imap_file;

imap_file *imap_open(const char *path);
imap_file *imap_open_io(imap_io *io, int own_io);
void       imap_close(imap_file *f);

const char *imap_profile(const imap_file *f);
const char *imap_schema_text(const imap_file *f);
uint32_t    imap_n_chunks(const imap_file *f);
uint32_t    imap_n_members(const imap_file *f, int axis);
const imap_member *imap_member_at(const imap_file *f, int axis, uint32_t i);
int32_t     imap_member_by_name(const imap_file *f, int axis, const char *name);
const imap_chunk  *imap_chunk_at(const imap_file *f, uint32_t i);

/* Decode chunk i into out[], which must hold imap_chunk_at(f,i)->n_runs runs. */
int imap_read_chunk(imap_file *f, uint32_t i, imap_run *out);

/* ------------------------------------------------------------------ queries */

/* One run clipped to a query window. [a, a+len) on axis a of a_member pairs with
 * [b, b+len) on axis b of b_member; strand 1 pairs a+i with b+len-1-i. */
typedef struct {
    uint32_t a_member, b_member;
    int64_t  a, b, len;
    uint8_t  strand;
} imap_hit;

typedef struct {
    uint32_t chunks_examined;   /* directory entries considered */
    uint32_t chunks_decoded;    /* chunks actually read and inflated */
} imap_query_stats;

/* Every run of axis-a member m overlapping [lo,hi) on axis a, clipped to the window.
 * Results are in axis-a order.  *out is malloc'd; caller frees.  stats may be NULL.
 * Axis a is a partial function, so results never overlap on axis a. */
int imap_query_a(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *stats);

/* Every run of axis-b member m overlapping [lo,hi) on axis b, clipped to the window.
 * Results are sorted by (b, a_member, a).  Axis b is a multimap, so results may
 * overlap on axis b -- that is paralogy, and it is reported, not resolved. */
int imap_query_b(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *stats);

#endif /* IMAP_FILE_H */
