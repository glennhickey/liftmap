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
#define IMAP_DIRA_ENTRY   72
#define IMAP_DIRB_ENTRY   12
#define IMAP_MEM_ENTRY    40

/* Format limits (SPEC 2.2).  A chunk's decode buffers are sized from its directory
 * entry, so these, together with the byte-count checks at open, are what stop a few
 * hundred bytes of hostile file from demanding gigabytes. */
#define IMAP_MAX_CHUNK_RUNS  (1u << 24)
#define IMAP_MAX_CHUNK_RAW   (1u << 30)
#define IMAP_ZLIB_MAX_RATIO  1032u      /* deflate's maximum expansion */

/* Header feature flags (SPEC 2).  A reader refuses any bit it does not know, so a
 * file using a newer feature is rejected by an older reader, never misread. */
#define IMAP_FEAT_ORDER_B   (1ull << 0)     /* runs within a chunk are in order b */
#define IMAP_FEAT_A_OVERLAP (1ull << 1)     /* axis a may overlap (order b only; SPEC 1.2) */
#define IMAP_FEAT_KNOWN     (IMAP_FEAT_ORDER_B | IMAP_FEAT_A_OVERLAP)

/* Chunk codecs (SPEC 2.2).  A chunk's four streams are compressed independently. */
#define IMAP_CODEC_NONE     0     /* every stream stored verbatim */
#define IMAP_CODEC_DEFLATE  1     /* raw deflate per stream, verbatim where it would not shrink */

/* On-disk format version.  There is no reader for earlier drafts: a file written
 * before a bump must be regenerated (the .tui policy). */
#define IMAP_FORMAT_MAJOR 4

/* One chunk's directory entry, in memory (native types; serialized explicitly). */
typedef struct {
    uint32_t a_member, b_member;
    uint64_t a_min;  uint32_t a_span;
    uint64_t b_min;  uint32_t b_span;
    uint64_t base;              /* decode base of the NON-ordering axis (SPEC 2.2):
                                 * order a: the first run's axis-b entry point;
                                 * order b: the first run's a.  The ordering axis's
                                 * base is its min, since the first run is the min. */
    uint64_t off;               /* within the runs section */
    uint32_t clen, rawlen, n_runs;
    uint32_t codec;            /* IMAP_CODEC_* for this chunk's blob */
    uint32_t crc;              /* CRC-32 of the chunk's clen stored bytes, checked on read */
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

/* Runs must arrive in axis-a order within a member, without overlap (unless
 * imap_writer_set_a_overlap), and lie inside both members' lengths; a violating run returns -1 and fails the writer.  The writer
 * closes chunks per SPEC S1.3 (count / bspan / member change) and never lets a chunk
 * span members. */
int imap_writer_add_run(imap_writer *w, uint32_t a_member, uint32_t b_member,
                        int64_t a, int64_t b, int64_t len, uint8_t strand);

int imap_writer_set_params(imap_writer *w, uint32_t count, uint64_t bspan, int codec);

/* Application section (SPEC 2.5): id must begin "x." and be unique.  Stored verbatim and
 * checksummed; the library never interprets it.  Data is copied. */
int imap_writer_add_section(imap_writer *w, const char *id, const void *data, size_t n);

/* IMAP_ORDER_A (default) or IMAP_ORDER_B.  Only before the first run.  Under order b the
 * writer buffers each axis-a member's runs, sorts them by (b_member, b, a) and cuts
 * chunks along axis b (SPEC S1.3); runs still arrive in axis-a order per member. */
int imap_writer_set_order(imap_writer *w, int order);

/* Declare that runs of one axis-a member may overlap on axis a (SPEC 1.2), for maps that
 * are approximate by construction -- e.g. runs merged across gaps, whose lengths no longer
 * match on both axes.  Recorded in the header so readers know axis a is not a partial
 * function.  Requires order b (set it first), and only before the first run. */
int imap_writer_set_a_overlap(imap_writer *w);
int imap_writer_close(imap_writer *w);          /* writes footer+trailer, frees */
void imap_writer_abort(imap_writer *w);

/* ------------------------------------------------------------------ reader */

typedef struct imap_file imap_file;

imap_file *imap_open(const char *path);
imap_file *imap_open_io(imap_io *io, int own_io);
void       imap_close(imap_file *f);

const char *imap_profile(const imap_file *f);
const char *imap_schema_text(const imap_file *f);
int         imap_order(const imap_file *f);          /* IMAP_ORDER_A or IMAP_ORDER_B */
int         imap_a_overlap(const imap_file *f);      /* 1 if axis a may overlap */
uint32_t    imap_n_chunks(const imap_file *f);
uint32_t    imap_n_members(const imap_file *f, int axis);
const imap_member *imap_member_at(const imap_file *f, int axis, uint32_t i);
int32_t     imap_member_by_name(const imap_file *f, int axis, const char *name);
const imap_chunk  *imap_chunk_at(const imap_file *f, uint32_t i);

/* Read an application section written with imap_writer_add_section, checksum-verified.
 * Returns 0 with *out malloc'd (caller frees), 1 if the file has no such section, -1 on
 * error. */
int imap_read_section(imap_file *f, const char *id, uint8_t **out, size_t *n);

/* Members whose name begins with `prefix` form a contiguous range in name order:
 * ranks [*first_rank, *first_rank + *count).  imap_member_by_rank maps a rank to a
 * member id. */
int     imap_member_prefix(const imap_file *f, int axis, const char *prefix,
                           uint32_t *first_rank, uint32_t *count);
int32_t imap_member_by_rank(const imap_file *f, int axis, uint32_t rank);

/* Decode chunk i into out[], which must hold imap_chunk_at(f,i)->n_runs runs.  The
 * chunk's bytes are checked against its directory CRC first; -1 on any mismatch. */
int imap_read_chunk(imap_file *f, uint32_t i, imap_run *out);

/* Full integrity check: the runs section's CRC, streamed in bounded reads, and every
 * chunk's own CRC.  Open does not read the payload -- it is the one section whose size
 * scales with the data -- so this is the way to check a whole file.  0 if intact. */
int imap_verify(imap_file *f);

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
 * Results are sorted by (a, b, len).  *out is malloc'd; caller frees.  stats may be NULL.
 * Axis a is a partial function, so results never overlap on axis a -- unless the file
 * declares IMAP_FEAT_A_OVERLAP (imap_a_overlap). */
int imap_query_a(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *stats);

/* Every run of axis-b member m overlapping [lo,hi) on axis b, clipped to the window.
 * Results are sorted by (b, a_member, a).  Axis b is a multimap, so results may
 * overlap on axis b -- that is paralogy, and it is reported, not resolved. */
int imap_query_b(imap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 imap_hit **out, size_t *n, imap_query_stats *stats);

#endif /* IMAP_FILE_H */
