/* lmap_file -- container framing for the interval map format. See doc/SPEC.md S2.
 *
 * Layout:   [64 B header] [sections...] [footer] [32 B trailer]
 * Open is two reads: the fixed trailer, then the footer it points at.
 * Every offset stored in the file is relative to the file's own byte 0.
 */
#ifndef LMAP_FILE_H
#define LMAP_FILE_H

#include "lmap_codec.h"
#include "lmap_io.h"

#define LMAP_MAGIC        "LMAP\x1A\x0D\x0A"   /* 7 chars + implicit NUL = 8 bytes */
#define LMAP_HEADER_SIZE  64
#define LMAP_TRAILER_SIZE 32
#define LMAP_DIRA_ENTRY   72
#define LMAP_DIRB_ENTRY   12
#define LMAP_MEM_ENTRY    40

/* Format limits (SPEC 2.2).  A chunk's decode buffers are sized from its directory
 * entry, so these, together with the byte-count checks at open, are what stop a few
 * hundred bytes of hostile file from demanding gigabytes. */
#define LMAP_MAX_CHUNK_RUNS  (1u << 24)
#define LMAP_MAX_CHUNK_RAW   (1u << 30)
#define LMAP_ZLIB_MAX_RATIO  1032u      /* deflate's maximum expansion */

/* Header feature flags (SPEC 2).  A reader refuses any bit it does not know, so a
 * file using a newer feature is rejected by an older reader, never misread. */
#define LMAP_FEAT_ORDER_B   (1ull << 0)     /* runs within a chunk are in order b */
#define LMAP_FEAT_A_OVERLAP (1ull << 1)     /* axis a may overlap (order b only; SPEC 1.2) */
#define LMAP_FEAT_KNOWN     (LMAP_FEAT_ORDER_B | LMAP_FEAT_A_OVERLAP)

/* Chunk codecs (SPEC 2.2).  A chunk's four streams are compressed independently. */
#define LMAP_CODEC_NONE     0     /* every stream stored verbatim */
#define LMAP_CODEC_DEFLATE  1     /* raw deflate per stream, verbatim where it would not shrink */

/* On-disk format version.  There is no reader for earlier drafts: a file written
 * before a bump must be regenerated.  Major 1 was libintervalmap's format 4 renamed;
 * major 2 replaced the free-text schema section with metadata and added groups. */
#define LMAP_FORMAT_MAJOR 2

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
    uint32_t codec;            /* LMAP_CODEC_* for this chunk's blob */
    uint32_t crc;              /* CRC-32 of the chunk's clen stored bytes, checked on read */
} lmap_chunk;

/* One member (sequence) of an axis. */
typedef struct {
    char    *name;
    uint64_t length;
    uint64_t axis_min, axis_max;
    uint32_t first_chunk, n_chunks;
} lmap_member;

/* ------------------------------------------------------------------ writer */

typedef struct lmap_writer lmap_writer;

#define LMAP_NO_GROUP 0xFFFFFFFFu

/* `profile` names what the file is, e.g. "hal2.edge" (at most 31 bytes kept).  The
 * file is written as `path`.partial and renamed to `path` only by a successful close;
 * a failed close or an abort removes it. */
lmap_writer *lmap_writer_open(const char *path, const char *profile);

/* Metadata (SPEC 2.5): key -> value strings, loaded at open.  Setting a key again
 * replaces it.  Keys are 1..255 bytes; see the SPEC for the standard keys (axis names,
 * max_gap, tree).  Applications prefix their own keys, e.g. "tui.format". */
int lmap_writer_set_meta(lmap_writer *w, const char *key, const char *value);

/* Groups (SPEC 2.4a): a named set of members of one axis, e.g. one genome's sequences.
 * add_group returns the group id (declaration order; names unique per axis);
 * set_member_group puts a member in it.  A member is in at most one group. */
int32_t lmap_writer_add_group(lmap_writer *w, int axis, const char *name);
int     lmap_writer_set_member_group(lmap_writer *w, int axis, uint32_t member, uint32_t group);

/* Members must be declared before any run referencing them. Returns member id. */
int32_t lmap_writer_add_member(lmap_writer *w, int axis, const char *name, uint64_t length);

/* Runs must arrive in axis-a order within a member, without overlap (unless
 * lmap_writer_set_a_overlap), and lie inside both members' lengths; a violating run returns -1 and fails the writer.  The writer
 * closes chunks per SPEC S1.3 (count / bspan / member change) and never lets a chunk
 * span members. */
int lmap_writer_add_run(lmap_writer *w, uint32_t a_member, uint32_t b_member,
                        int64_t a, int64_t b, int64_t len, uint8_t strand);

int lmap_writer_set_params(lmap_writer *w, uint32_t count, uint64_t bspan, int codec);

/* Application section (SPEC 2.5): id must begin "x." and be unique.  Stored verbatim and
 * checksummed; the library never interprets it.  Data is copied. */
int lmap_writer_add_section(lmap_writer *w, const char *id, const void *data, size_t n);

/* LMAP_ORDER_A (default) or LMAP_ORDER_B.  Only before the first run.  Under order b the
 * writer buffers each axis-a member's runs, sorts them by (b_member, b, a) and cuts
 * chunks along axis b (SPEC S1.3); runs still arrive in axis-a order per member. */
int lmap_writer_set_order(lmap_writer *w, int order);

/* Declare that runs of one axis-a member may overlap on axis a (SPEC 1.2), for maps that
 * are approximate by construction -- e.g. runs merged across gaps, whose lengths no longer
 * match on both axes.  Recorded in the header so readers know axis a is not a partial
 * function.  Requires order b (set it first), and only before the first run. */
int lmap_writer_set_a_overlap(lmap_writer *w);
int lmap_writer_close(lmap_writer *w);          /* writes footer+trailer, frees */
void lmap_writer_abort(lmap_writer *w);

/* ------------------------------------------------------------------ reader */

typedef struct lmap_file lmap_file;

lmap_file *lmap_open(const char *path);
lmap_file *lmap_open_io(lmap_io *io, int own_io);
void       lmap_close(lmap_file *f);

const char *lmap_profile(const lmap_file *f);
const char *lmap_meta_get(const lmap_file *f, const char *key);    /* NULL if absent */
uint32_t    lmap_meta_count(const lmap_file *f);                    /* in key order */
int         lmap_meta_at(const lmap_file *f, uint32_t i, const char **key, const char **value);

/* Groups of one axis: ids 0..n-1 in declaration order.  lmap_group_members gives the
 * group's member ids in name order (borrowed, valid until close), their count and
 * summed length -- pass the ids straight to lmap_cursor_open to restrict a cursor to
 * the group.  lmap_member_group is -1 for a member in no group. */
uint32_t    lmap_n_groups(const lmap_file *f, int axis);
const char *lmap_group_name(const lmap_file *f, int axis, uint32_t g);
int32_t     lmap_group_by_name(const lmap_file *f, int axis, const char *name);
int32_t     lmap_member_group(const lmap_file *f, int axis, uint32_t m);
int         lmap_group_members(const lmap_file *f, int axis, uint32_t g,
                               const uint32_t **members, uint32_t *n, uint64_t *total_length);
int         lmap_order(const lmap_file *f);          /* LMAP_ORDER_A or LMAP_ORDER_B */
int         lmap_a_overlap(const lmap_file *f);      /* 1 if axis a may overlap */
uint32_t    lmap_n_chunks(const lmap_file *f);
uint32_t    lmap_n_members(const lmap_file *f, int axis);
const lmap_member *lmap_member_at(const lmap_file *f, int axis, uint32_t i);
int32_t     lmap_member_by_name(const lmap_file *f, int axis, const char *name);
const lmap_chunk  *lmap_chunk_at(const lmap_file *f, uint32_t i);

/* Read an application section written with lmap_writer_add_section, checksum-verified.
 * Returns 0 with *out malloc'd (caller frees), 1 if the file has no such section, -1 on
 * error. */
int lmap_read_section(lmap_file *f, const char *id, uint8_t **out, size_t *n);

/* Members whose name begins with `prefix` form a contiguous range in name order:
 * ranks [*first_rank, *first_rank + *count).  lmap_member_by_rank maps a rank to a
 * member id. */
int     lmap_member_prefix(const lmap_file *f, int axis, const char *prefix,
                           uint32_t *first_rank, uint32_t *count);
int32_t lmap_member_by_rank(const lmap_file *f, int axis, uint32_t rank);

/* Decode chunk i into out[], which must hold lmap_chunk_at(f,i)->n_runs runs.  The
 * chunk's bytes are checked against its directory CRC first; -1 on any mismatch. */
int lmap_read_chunk(lmap_file *f, uint32_t i, lmap_run *out);

/* Full integrity check: the runs section's CRC, streamed in bounded reads, and every
 * chunk's own CRC.  Open does not read the payload -- it is the one section whose size
 * scales with the data -- so this is the way to check a whole file.  0 if intact. */
int lmap_verify(lmap_file *f);

/* ------------------------------------------------------------------ queries */

/* One run clipped to a query window. [a, a+len) on axis a of a_member pairs with
 * [b, b+len) on axis b of b_member; strand 1 pairs a+i with b+len-1-i. */
typedef struct {
    uint32_t a_member, b_member;
    int64_t  a, b, len;
    uint8_t  strand;
} lmap_hit;

typedef struct {
    uint32_t chunks_examined;   /* directory entries considered */
    uint32_t chunks_decoded;    /* chunks actually read and inflated */
} lmap_query_stats;

/* Every run of axis-a member m overlapping [lo,hi) on axis a, clipped to the window.
 * Results are sorted by (a, b, len).  *out is malloc'd; caller frees.  stats may be NULL.
 * Axis a is a partial function, so results never overlap on axis a -- unless the file
 * declares LMAP_FEAT_A_OVERLAP (lmap_a_overlap). */
int lmap_query_a(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *stats);

/* Every run of axis-b member m overlapping [lo,hi) on axis b, clipped to the window.
 * Results are sorted by (b, a_member, a).  Axis b is a multimap, so results may
 * overlap on axis b -- that is paralogy, and it is reported, not resolved. */
int lmap_query_b(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *stats);

/* ------------------------------------------------------------------ builder
 *
 * A front end to the writer for runs that arrive in any order.  Runs may repeat, overlap,
 * or split one colinear stretch into pieces: at close they are reduced to the canonical
 * form of the aligned base-pair set -- maximal runs, unioned along each diagonal -- so the
 * same base pairs always give the same file.  Each axis-a member's runs are then sorted
 * on a and the storage order chosen: as set (default order a), except that if axis a
 * overlaps the file is written in order b with LMAP_FEAT_A_OVERLAP, or close fails if
 * overlap was not allowed.
 *
 * Memory: runs are held up to mem_bytes (0 = 1 GiB), then spilled as sorted segments to
 * a temporary file in tmp_dir (NULL = the output's directory).  The file is unlinked on
 * creation, so nothing is left behind even on a crash.  Close holds one axis-a member's
 * runs at a time.
 *
 * Members are declared through the builder (it checks runs against their lengths).
 * Metadata, groups, sections and chunk parameters go straight to the underlying writer,
 * lmap_builder_writer(b); do not add members, runs, or set the order on it directly. */

#define LMAP_ORDER_AUTO (-1)

typedef struct lmap_builder lmap_builder;

typedef struct {
    uint64_t input_runs;          /* as added */
    uint64_t runs;                /* canonical, as written */
    uint64_t overlapping_runs;    /* runs overlapping an earlier one of their member on a */
    uint32_t spill_segments;
    int      order;               /* the order written */
    char     error[512];          /* why close failed, "" on success */
} lmap_build_stats;

lmap_builder *lmap_builder_open(const char *path, const char *profile,
                                size_t mem_bytes, const char *tmp_dir);
lmap_writer  *lmap_builder_writer(lmap_builder *b);
int32_t lmap_builder_add_member(lmap_builder *b, int axis, const char *name, uint64_t length);
int     lmap_builder_add_run(lmap_builder *b, uint32_t a_member, uint32_t b_member,
                             int64_t a, int64_t bpos, int64_t len, uint8_t strand);
int     lmap_builder_set_order(lmap_builder *b, int order);   /* LMAP_ORDER_A/B/AUTO */
int     lmap_builder_allow_overlap(lmap_builder *b, int allow);
/* Writes the file and frees b: 0 on success, -1 on failure with the reason in
 * stats->error (pass stats to get it).  Any earlier call that returned -1 failed the
 * builder; lmap_builder_error says why.  Abort discards everything. */
int     lmap_builder_close(lmap_builder *b, lmap_build_stats *stats);
void    lmap_builder_abort(lmap_builder *b);
const char *lmap_builder_error(const lmap_builder *b);   /* last failure, "" if none */

/* ------------------------------------------------------------------ cursor
 *
 * A cursor walks one member of a KEY axis -- in key order, across chunks that overlap --
 * optionally restricted to a set of members of the other axis (e.g. one genome's
 * sequences).  It caches decoded chunks, so repeated point lookups and windows in one
 * region decode each chunk once.
 *
 *   cache_bytes   budget for decoded chunks kept between calls, evicted least recently
 *                 used first.  0 streams: a chunk is dropped as soon as the walk passes
 *                 it, so a full scan holds only the chunks overlapping its position.
 *                 LMAP_CACHE_ALL keeps everything decoded.
 *
 * Hits are whole runs (not clipped to the window; see lmap_hit_clip), ordered by
 * (key position, other member, other position, len, strand).  A cursor is not
 * thread-safe; the file it reads is, so use one cursor per thread. */

#define LMAP_CACHE_ALL ((size_t)-1)

typedef struct lmap_cursor lmap_cursor;

typedef struct {
    uint64_t chunks_examined;     /* directory entries considered */
    uint64_t chunks_decoded;      /* chunks read and inflated */
    uint64_t chunk_reuses;        /* lookups served from the cache */
    size_t   cached_bytes, peak_cached_bytes;
} lmap_cursor_stats;

/* other = NULL for every member of the other axis, else n_other member ids to keep.
 * NULL on bad arguments or no memory. */
lmap_cursor *lmap_cursor_open(lmap_file *f, int key_axis, uint32_t key_member,
                              const uint32_t *other, uint32_t n_other, size_t cache_bytes);

/* Iterate the runs overlapping [lo, hi) on the key axis: seek, then next until it
 * returns 0 (1 = *out filled, -1 = error).  Seek again at any time. */
int  lmap_cursor_seek(lmap_cursor *c, int64_t lo, int64_t hi);
int  lmap_cursor_next(lmap_cursor *c, lmap_hit *out);

/* Every run covering key position pos.  Fills up to cap hits and returns how many
 * there are in total (may exceed cap; then the ones filled are an arbitrary subset,
 * so retry with a larger buffer), or -1 on error. */
int  lmap_cursor_point(lmap_cursor *c, int64_t pos, lmap_hit *out, int cap);

/* Key range the cursor's chunks cover, [*lo, *hi); 0 if it has none (then *lo = *hi = 0).
 * From the directory: nothing is decoded. */
int  lmap_cursor_extent(const lmap_cursor *c, int64_t *lo, int64_t *hi);
uint32_t lmap_cursor_n_chunks(const lmap_cursor *c);

void lmap_cursor_get_stats(const lmap_cursor *c, lmap_cursor_stats *stats);
void lmap_cursor_close(lmap_cursor *c);

/* Clip a hit to [lo, hi) on `axis`, mapping the other axis through the strand.  Returns
 * 1 with *out filled, 0 if the hit does not overlap the window. */
int  lmap_hit_clip(const lmap_hit *h, int axis, int64_t lo, int64_t hi, lmap_hit *out);

#endif /* LMAP_FILE_H */
