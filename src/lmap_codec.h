/* lmap_codec -- run codec for the interval map format.  See doc/SPEC.md S1.2.
 *
 * A chunk's runs are encoded as four independent streams, in declaration order:
 *
 *   a       uvarint   a - prev_a_exit          (prev_a_exit = prev a + prev len)
 *   b       zigzag    b_enter - prev_b_exit    (strand-aware, see lmap_b_enter/exit)
 *   len     uvarint   len - 1
 *   strand  bitmap    1 bit per run, LSB first
 *
 * Every delta restarts at zero at a chunk boundary, which is what makes chunks
 * independently decodable; the bases come from the chunk directory.
 */
#ifndef LMAP_CODEC_H
#define LMAP_CODEC_H

#include <stddef.h>
#include <stdint.h>

#define LMAP_STREAM_A       0
#define LMAP_STREAM_B       1
#define LMAP_STREAM_LEN     2
#define LMAP_STREAM_STRAND  3
#define LMAP_N_STREAMS      4

/* Run order within a chunk (SPEC 1.3).  Order a: axis a is a uvarint delta from the
 * previous run's end and axis b follows the traversal (strand-aware).  Order b: runs
 * sorted by (b, a); both axes are zigzag deltas from the previous run's end. */
#include "liftmap.h"      /* LMAP_ORDER_A / LMAP_ORDER_B */

/* One run, in the coordinates of its chunk's members. */
typedef struct {
    int64_t a;
    int64_t b;
    int64_t len;
    uint8_t strand;   /* 0 forward, 1 reverse */
} lmap_run;

/* A buffer the encoder appends to; caller frees .p with free(). */
typedef struct {
    uint8_t *p;
    size_t   n;     /* bytes used */
    size_t   cap;
} lmap_buf;

/* --- varint primitives (exposed for tests) --- */
size_t  lmap_put_uvarint(uint8_t *dst, uint64_t v);   /* dst needs >= 10 bytes */
int     lmap_get_uvarint(const uint8_t **pp, const uint8_t *end, uint64_t *out);

static inline uint64_t lmap_zigzag(int64_t v)   { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
static inline int64_t  lmap_unzigzag(uint64_t z){ return (int64_t)(z >> 1) ^ -(int64_t)(z & 1); }

/* --- strand-aware axis-b traversal points (SPEC S1.2) --- */
static inline int64_t lmap_b_enter(int64_t b, int64_t len, uint8_t s) { return s ? b + len : b; }
static inline int64_t lmap_b_exit (int64_t b, int64_t len, uint8_t s) { return s ? b : b + len; }

/* --- chunk encode / decode ---
 * encode: runs must already be in `order` (the writer sorts); fills streams[], each
 *         buffer malloc'd, caller frees.  Returns 0, or -1 on bad input or no memory.
 * decode: writes n runs into out[].  base_a is the first run's a.  base_b is the first
 *         run's axis-b entry point: lmap_b_enter(b,len,strand) under order a, plain b
 *         under order b.  Returns 0, or -1 on malformed input. */
int lmap_encode_chunk(const lmap_run *runs, int64_t n, int order,
                      lmap_buf streams[LMAP_N_STREAMS]);
int lmap_decode_chunk(const lmap_buf streams[LMAP_N_STREAMS], int64_t n, int order,
                      int64_t base_a, int64_t base_b, lmap_run *out);

void lmap_buf_free(lmap_buf *b);

#endif /* LMAP_CODEC_H */
