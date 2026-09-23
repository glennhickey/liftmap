/* lmap_file -- liftmap INTERNALS: the container framing (doc/SPEC.md S2), the writer the
 * builder drives, and the chunk directory.  For the library and its tests; the public
 * API is liftmap.h.
 *
 * Layout:   [64 B header] [sections...] [footer] [32 B trailer]
 * Every offset stored in the file is relative to the file's own byte 0.
 */
#ifndef LMAP_FILE_H
#define LMAP_FILE_H

#include "liftmap.h"
#include "lmap_codec.h"

#define LMAP_MAGIC        "LMAP\x1A\x0D\x0A"   /* 7 chars + implicit NUL = 8 bytes */
#define LMAP_HEADER_SIZE  64
#define LMAP_TRAILER_SIZE 32
#define LMAP_DIRA_ENTRY   72
#define LMAP_DIRB_ENTRY   12
#define LMAP_DIR_PAGE     4096u   /* directory entries per page, each with its own CRC */
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


/* A member's range in its axis's directory: dir.a positions for axis a, dir.b for b. */
typedef struct { uint32_t first, n; } lmap_mrange;

/* ------------------------------------------------------------------ writer
 * Driven by the builder (lmap_build.c); files are written through the builder. */

typedef struct lmap_writer lmap_writer;


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

/* ------------------------------------------------------------------ chunk directory */

uint32_t          lmap_n_chunks(const lmap_file *f);
const lmap_chunk *lmap_chunk_at(const lmap_file *f, uint32_t i);

/* Decode chunk i into out[], which must hold lmap_chunk_at(f,i)->n_runs runs.  The
 * chunk's bytes are checked against its directory CRC first; -1 on any mismatch. */
int lmap_read_chunk(lmap_file *f, uint32_t i, lmap_run *out);

#endif /* LMAP_FILE_H */
