/* imap_codec -- run codec for the interval map format.  See doc/SPEC.md S1.2.
 *
 * A chunk's runs are encoded as four independent streams, in declaration order:
 *
 *   a       uvarint   a - prev_a_exit          (prev_a_exit = prev a + prev len)
 *   b       zigzag    b_enter - prev_b_exit    (strand-aware, see imap_b_enter/exit)
 *   len     uvarint   len - 1
 *   strand  bitmap    1 bit per run, LSB first
 *
 * Every delta restarts at zero at a chunk boundary, which is what makes chunks
 * independently decodable; the bases come from the chunk directory.
 */
#ifndef IMAP_CODEC_H
#define IMAP_CODEC_H

#include <stddef.h>
#include <stdint.h>

#define IMAP_STREAM_A       0
#define IMAP_STREAM_B       1
#define IMAP_STREAM_LEN     2
#define IMAP_STREAM_STRAND  3
#define IMAP_N_STREAMS      4

/* One run, in the coordinates of its chunk's members. */
typedef struct {
    int64_t a;
    int64_t b;
    int64_t len;
    uint8_t strand;   /* 0 forward, 1 reverse */
} imap_run;

/* A buffer the encoder appends to; caller frees .p with free(). */
typedef struct {
    uint8_t *p;
    size_t   n;     /* bytes used */
    size_t   cap;
} imap_buf;

/* --- varint primitives (exposed for tests) --- */
size_t  imap_put_uvarint(uint8_t *dst, uint64_t v);   /* dst needs >= 10 bytes */
int     imap_get_uvarint(const uint8_t **pp, const uint8_t *end, uint64_t *out);

static inline uint64_t imap_zigzag(int64_t v)   { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
static inline int64_t  imap_unzigzag(uint64_t z){ return (int64_t)(z >> 1) ^ -(int64_t)(z & 1); }

/* --- strand-aware axis-b traversal points (SPEC S1.2) --- */
static inline int64_t imap_b_enter(int64_t b, int64_t len, uint8_t s) { return s ? b + len : b; }
static inline int64_t imap_b_exit (int64_t b, int64_t len, uint8_t s) { return s ? b : b + len; }

/* --- chunk encode / decode ---
 * encode: fills streams[IMAP_N_STREAMS]; each buffer is malloc'd, caller frees.
 *         returns 0 on success, -1 on allocation failure.
 * decode: writes n runs into out[]; a0/b_enter0 are the chunk bases.
 *         returns 0 on success, -1 on malformed input. */
int imap_encode_chunk(const imap_run *runs, int64_t n, imap_buf streams[IMAP_N_STREAMS]);
int imap_decode_chunk(const imap_buf streams[IMAP_N_STREAMS], int64_t n,
                      int64_t a0, int64_t b_enter0, imap_run *out);

void imap_buf_free(imap_buf *b);

#endif /* IMAP_CODEC_H */
