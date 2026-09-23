/* liftmap -- indexed pairwise coordinate maps between genomes.  The public API.
 *
 * A liftmap is a set of colinear ungapped runs between two axes, each made of named
 * members (sequences):
 *
 *     (a_member, a, b_member, b, len, strand)    [a, a+len) <-> [b, b+len)
 *
 * with strand 1 pairing a+i with b+len-1-i.  Axis a is a function (each base maps at
 * most once) unless the file declares otherwise; axis b is a multimap (paralogy is
 * reported, not resolved).  Files are written with a builder, read with a cursor or
 * the one-shot queries, and described by metadata and groups.  See doc/SPEC.md.
 *
 * Errors: calls return -1 (or NULL) on failure.  The builder keeps a message
 * (lmap_builder_error).  Reading: a file is refused at open if any structure it
 * depends on is inconsistent, and a chunk whose bytes fail their checksum is refused
 * when read.
 *
 * Threads: an open lmap_file is read-only and safe to share; cursors and builders
 * are not -- one per thread. */
#ifndef LIFTMAP_H
#define LIFTMAP_H

#include <stddef.h>
#include <stdint.h>

#include "lmap_io.h"

#define LMAP_FORMAT_MAJOR 3

#define LMAP_ORDER_A     0          /* chunks cut and runs stored along axis a */
#define LMAP_ORDER_B     1          /* along axis b */
#define LMAP_ORDER_AUTO  (-1)       /* builder: a, or b when axis a overlaps */

#define LMAP_CODEC_NONE     0       /* chunk streams stored verbatim */
#define LMAP_CODEC_DEFLATE  1       /* raw deflate per stream where it shrinks (default) */

#define LMAP_NO_GROUP  0xFFFFFFFFu
#define LMAP_CACHE_ALL ((size_t)-1)

/* ------------------------------------------------------------------ reading */

typedef struct lmap_file lmap_file;

/* One member (sequence) of an axis.  axis_min/axis_max: the extent its runs cover. */
typedef struct {
    char    *name;
    uint64_t length;
    uint64_t axis_min, axis_max;
} lmap_member;

lmap_file *lmap_open(const char *path);
lmap_file *lmap_open_io(lmap_io *io, int own_io);     /* any lmap_io: memory, slice, ... */
void       lmap_close(lmap_file *f);

const char *lmap_profile(const lmap_file *f);        /* e.g. "taffy.tui", "hal2.edge" */
int         lmap_order(const lmap_file *f);           /* LMAP_ORDER_A or LMAP_ORDER_B */
int         lmap_a_overlap(const lmap_file *f);       /* 1 if axis a may overlap */

/* Members of one axis (0 = a, 1 = b): ids 0..n-1.  Names are unique per axis; by_rank
 * walks them in name order, and a prefix names a contiguous range of ranks. */
uint32_t           lmap_n_members(const lmap_file *f, int axis);
const lmap_member *lmap_member_at(const lmap_file *f, int axis, uint32_t i);
int32_t            lmap_member_by_name(const lmap_file *f, int axis, const char *name);
int32_t            lmap_member_by_rank(const lmap_file *f, int axis, uint32_t rank);
int                lmap_member_prefix(const lmap_file *f, int axis, const char *prefix,
                                      uint32_t *first_rank, uint32_t *count);

/* Metadata: string key -> value, loaded at open.  Standard keys: axis.a, axis.b (what
 * each axis is), max_gap (non-zero: coordinates approximate), tree, source, program,
 * description.  Applications prefix their own ("tui.format"). */
const char *lmap_meta_get(const lmap_file *f, const char *key);    /* NULL if absent */
uint32_t    lmap_meta_count(const lmap_file *f);                    /* in key order */
int         lmap_meta_at(const lmap_file *f, uint32_t i, const char **key, const char **value);

/* Groups: named sets of members of one axis, e.g. one genome's sequences; ids 0..n-1 in
 * declaration order.  lmap_group_members gives the member ids in name order (borrowed,
 * valid until close), their count and summed length -- pass the ids to
 * lmap_cursor_open to restrict a cursor to the group.  lmap_member_group is -1 for a
 * member in no group. */
uint32_t    lmap_n_groups(const lmap_file *f, int axis);
const char *lmap_group_name(const lmap_file *f, int axis, uint32_t g);
int32_t     lmap_group_by_name(const lmap_file *f, int axis, const char *name);
int32_t     lmap_member_group(const lmap_file *f, int axis, uint32_t m);
int         lmap_group_members(const lmap_file *f, int axis, uint32_t g,
                               const uint32_t **members, uint32_t *n, uint64_t *total_length);

/* Application sections ("x." ids): bytes the library stores and checksums but never
 * interprets.  lmap_read_section returns 0 with *out malloc'd (caller frees), 1 if
 * absent, -1 on error. */
int         lmap_read_section(lmap_file *f, const char *id, uint8_t **out, size_t *n);
uint32_t    lmap_n_app_sections(const lmap_file *f);
const char *lmap_app_section_id(const lmap_file *f, uint32_t i);

typedef struct {
    uint64_t runs;               /* runs stored */
    uint64_t payload_bytes;      /* compressed run data */
    uint32_t chunks;
} lmap_file_stats;
void lmap_get_stats(const lmap_file *f, lmap_file_stats *stats);    /* from the directory */

/* Full integrity check: every checksum, and every chunk decoded and checked against its
 * directory entry.  Open never reads the run payload, so this is how to check it all.
 * 0 if intact. */
int lmap_verify(lmap_file *f);

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
 * (key position, other member, other position, len, strand). */

/* A run: [a, a+len) on axis a of a_member pairs with [b, b+len) on axis b of b_member;
 * strand 1 pairs a+i with b+len-1-i. */
typedef struct {
    uint32_t a_member, b_member;
    int64_t  a, b, len;
    uint8_t  strand;
} lmap_hit;

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
 * so retry with a larger buffer), or -1 on error.  A sweep of increasing pos costs
 * O(1) per call. */
int  lmap_cursor_point(lmap_cursor *c, int64_t pos, lmap_hit *out, int cap);

/* Key range the cursor's chunks cover, [*lo, *hi); 0 if it has none.  From the
 * directory: nothing is decoded. */
int      lmap_cursor_extent(const lmap_cursor *c, int64_t *lo, int64_t *hi);
uint32_t lmap_cursor_n_chunks(const lmap_cursor *c);

void lmap_cursor_get_stats(const lmap_cursor *c, lmap_cursor_stats *stats);
void lmap_cursor_close(lmap_cursor *c);

/* Clip a hit to [lo, hi) on `axis`, mapping the other axis through the strand.  Returns
 * 1 with *out filled, 0 if the hit does not overlap the window. */
int  lmap_hit_clip(const lmap_hit *h, int axis, int64_t lo, int64_t hi, lmap_hit *out);

/* ------------------------------------------------------------------ one-shot queries */

typedef struct {
    uint32_t chunks_examined;
    uint32_t chunks_decoded;
} lmap_query_stats;

/* Every run of axis-a member m overlapping [lo,hi) on axis a, clipped to the window,
 * sorted by (a, b, len, strand).  *out is malloc'd; caller frees.  stats may be NULL.
 * A streaming cursor underneath; use a cursor directly for repeated lookups. */
int lmap_query_a(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *stats);

/* The same from axis b, sorted by (b, a_member, a).  Results may overlap on axis b --
 * that is paralogy, and it is reported, not resolved. */
int lmap_query_b(lmap_file *f, uint32_t m, int64_t lo, int64_t hi,
                 lmap_hit **out, size_t *n, lmap_query_stats *stats);

/* ------------------------------------------------------------------ chaining and paralogy
 *
 * Chain aligned spans into collinear chains and keep the best copy where copies overlap
 * -- how liftOver picks one target among paralogs, and taffy's lift/view/blockViz dupe
 * filter.  Spans are general: the two extents may differ in length (a gapped block),
 * so a lifted run, a MAF block or a PAF record can all be chained.
 *
 * lmap_chain partitions spans by (q_member, strand) and, sweeping each in q order, lets
 * a span extend the best chain ending before it on both axes (gaps within max_gap, gap
 * cost below the span's own score): cost = chain_open + chain_extend * (q_gap + t_gap)
 * for a non-zero gap.  Chains are then claimed best first, and a chain is cut where a
 * link was already claimed.  chain_id[i] (1-based) is span i's chain; *chains is sorted
 * by score, best first.  The input is not reordered.
 *
 * lmap_chain_select walks the chains best first and keeps one unless its q coverage
 * overlaps the kept chains' by more than overlap_frac of its own: a paralog (the same
 * q bases landing again elsewhere) overlaps ~100% and drops; an inversion or a disjoint
 * piece overlaps 0% and stays.  keep[id] is set for kept chains (keep_len > max id,
 * zeroed by the caller); cap > 0 stops after that many. */

typedef struct {
    uint32_t q_member, t_member;
    int64_t  q_start, q_end;      /* forward coordinates on both axes */
    int64_t  t_start, t_end;
    uint8_t  strand;              /* 1: q reversed relative to t */
    int64_t  score;               /* typically the aligned length */
} lmap_span;

typedef struct {
    int64_t chain_open, chain_extend;
    int64_t max_gap;              /* a gap beyond this on either axis breaks a chain */
} lmap_chain_params;
#define LMAP_CHAIN_PARAMS_DEFAULT { 0, 1, 10000000 }

typedef struct {
    int64_t id, score, bp, n_spans;   /* bp: summed q extent */
} lmap_chain_info;

int lmap_chain(const lmap_span *spans, size_t n, const lmap_chain_params *params,
               int64_t *chain_id, lmap_chain_info **chains, size_t *n_chains);
int lmap_chain_select(const lmap_span *spans, size_t n, const int64_t *chain_id,
                      const lmap_chain_info *chains, size_t n_chains, double overlap_frac,
                      size_t cap, uint8_t *keep, size_t keep_len);

/* ------------------------------------------------------------------ lift
 *
 * Lift [lo, hi) of member m on axis `from` to the other axis.  The aligned pieces are
 * grouped into target intervals: pieces on the same target member and strand merge when
 * the gaps between them are within [0, max_gap] on both axes (walking the source
 * forward, the target forward on strand 0 and backward on strand 1) -- the same rule
 * lmap_coarsen chains by.  max_gap < 0 merges nothing: one interval per piece.  Results
 * are in source order (of their first piece); *out is malloc'd.
 * aligned_bp counts the bases actually aligned inside an interval, so aligned_bp /
 * (hi - lo) is liftOver's -minMatch fraction.
 *
 * With overlap_frac >= 0 the pieces are first chained (on the source axis, with
 * `chain`) and only the chains lmap_chain_select keeps are lifted: the best copy where
 * the source lands more than once -- liftOver's single-target behaviour at 0. */
typedef struct {
    uint32_t member;              /* on the target axis */
    int64_t  start, end;          /* target interval */
    uint8_t  strand;              /* 1: the source interval lands reversed */
    int64_t  src_start, src_end;  /* the part of [lo, hi) it came from */
    int64_t  aligned_bp;
} lmap_lifted;

typedef struct {
    int64_t max_gap;              /* merge across gaps up to this; < 0: one interval per piece */
    double  overlap_frac;         /* paralogy filter; < 0: off, 0: strict */
    lmap_chain_params chain;      /* for the filter's chaining */
} lmap_lift_opts;
#define LMAP_LIFT_OPTS_DEFAULT { -1, -1.0, LMAP_CHAIN_PARAMS_DEFAULT }

int lmap_lift(lmap_file *f, int from, uint32_t m, int64_t lo, int64_t hi,
              const lmap_lift_opts *opts, lmap_lifted **out, size_t *n);

/* ------------------------------------------------------------------ writing: the builder
 *
 * The one way to write a liftmap file.  Runs arrive in any order.  Runs may repeat,
 * overlap, or split one colinear stretch into pieces: at close they are reduced to the
 * canonical form of the aligned base-pair set -- maximal runs, unioned along each
 * diagonal -- so the same base pairs always give the same file.  Each axis-a member's
 * runs are then sorted on a and the storage order chosen: as set (default order a),
 * except that if axis a overlaps the file is written in order b with the overlap flag,
 * or close fails if overlap was not allowed.
 *
 * Memory: runs are held up to mem_bytes (0 = 1 GiB), then spilled as sorted segments to
 * a temporary file in tmp_dir (NULL = the output's directory), unlinked on creation so
 * nothing is left behind even on a crash.  Close holds the budget plus one axis-a
 * member's runs.  The output is written as `path`.partial and renamed into place only
 * by a successful close.
 *
 * Pre-sorted mode (lmap_builder_set_presorted, before the first run) is for callers
 * whose runs already arrive as the sort would leave them -- each axis-a member's runs
 * together, in increasing a.  Nothing is buffered or spilled: runs are merged with a
 * predecessor they continue and streamed to the file.  The order is fixed at the first
 * run (order b with the overlap flag if overlap is allowed, since overlap cannot be
 * known in advance), and a run that breaks the promise fails the build.
 *
 * Members must be declared before runs that name them; metadata, groups and sections
 * may be set any time before close.  Every call returns -1 on failure, fails the
 * builder, and leaves the reason in lmap_builder_error. */

typedef struct lmap_builder lmap_builder;

typedef struct {
    uint64_t input_runs;          /* as added */
    uint64_t runs;                /* canonical, as written */
    uint64_t overlapping_runs;    /* runs overlapping an earlier one of their member on a */
    uint32_t spill_segments;
    int      order;               /* the order written */
    char     error[512];          /* why close failed, "" on success */
} lmap_build_stats;

/* `profile` names what the file is, e.g. "hal2.edge" (at most 31 bytes kept). */
lmap_builder *lmap_builder_open(const char *path, const char *profile,
                                size_t mem_bytes, const char *tmp_dir);
int32_t lmap_builder_add_member(lmap_builder *b, int axis, const char *name, uint64_t length);
int     lmap_builder_add_run(lmap_builder *b, uint32_t a_member, uint32_t b_member,
                             int64_t a, int64_t bpos, int64_t len, uint8_t strand);
int     lmap_builder_set_meta(lmap_builder *b, const char *key, const char *value);
int32_t lmap_builder_add_group(lmap_builder *b, int axis, const char *name);
int     lmap_builder_set_member_group(lmap_builder *b, int axis, uint32_t member, uint32_t group);
int     lmap_builder_add_section(lmap_builder *b, const char *id, const void *data, size_t n);
/* chunk caps: at most `count` runs and `bspan` of axis b per chunk (defaults 8192, 10^6) */
int     lmap_builder_set_params(lmap_builder *b, uint32_t count, uint64_t bspan, int codec);
int     lmap_builder_set_order(lmap_builder *b, int order);   /* LMAP_ORDER_A/B/AUTO */
int     lmap_builder_allow_overlap(lmap_builder *b, int allow);
int     lmap_builder_set_presorted(lmap_builder *b, int on);
/* Writes the file and frees b: 0 on success, -1 on failure with the reason in
 * stats->error (pass stats to get it).  Abort discards everything. */
int     lmap_builder_close(lmap_builder *b, lmap_build_stats *stats);
void    lmap_builder_abort(lmap_builder *b);
const char *lmap_builder_error(const lmap_builder *b);   /* last failure, "" if none */

/* Declare f's members (same ids, both axes), groups, metadata and application sections
 * on a fresh builder: the start of any file derived from f. */
int lmap_builder_copy_layout(lmap_builder *b, lmap_file *f);

/* ------------------------------------------------------------------ coarsening
 *
 * Chain f's runs across small gaps into fewer, longer, APPROXIMATE runs -- a zoomed-out
 * level of detail (taffy's tui-chain, a hal2 LOD sidecar).  Along each member of
 * key_axis, in key order, each member of the other axis keeps one open chain; a run
 * extends it when it has the same strand and the gap to it is between 0 and max_gap on
 * both axes.  A chain's length is its span on the key axis, so its extent on the other
 * axis is approximate (it drifts wherever the two gaps it bridged differ), and chains of
 * one member may overlap on the other axis.  A chain is clipped where it would pass the
 * other member's end.  max_gap 0 merges only exact continuations.
 *
 * The chains go to b, which should hold f's layout (lmap_builder_copy_layout); coarsen
 * allows axis-a overlap on it and raises its max_gap metadata to max_gap. */
typedef struct {
    uint64_t runs_in, chains_out;
    uint64_t clipped_bp, dropped;   /* bases cut at the other member's end; chains lost whole */
} lmap_coarsen_stats;

int lmap_coarsen(lmap_file *f, lmap_builder *b, int key_axis, int64_t max_gap,
                 lmap_coarsen_stats *stats);

#endif /* LIFTMAP_H */
