"""Reference implementation of the imap run codec and chunking (draft 0.1).

Deliberately literal and slow: it follows doc/SPEC.md line by line so that any divergence
between the spec text and working code is visible. It is the thing the C library must agree
with, not a performance target.

Covers SPEC.md sections 1.2 (field semantics) and 1.3 (chunking). The container framing of
section 2 is not implemented here.
"""
import numpy as np
import zlib

__all__ = ["u_enc", "u_dec", "zz_enc", "zz_dec", "bitmap_enc", "bitmap_dec",
           "b_enter_exit", "chunk_bounds", "encode_chunk", "decode_chunk",
           "compress_streams", "decompress_streams"]


# --- encodings (SPEC 1.2) ---------------------------------------------------

def u_enc(vals):
    """LEB128 of non-negative values."""
    out = bytearray()
    for v in vals:
        v = int(v)
        if v < 0:
            raise ValueError("uvarint of negative value %d" % v)
        while v >= 0x80:
            out.append((v & 0x7F) | 0x80)
            v >>= 7
        out.append(v)
    return bytes(out)


def u_dec(buf, n):
    out = np.empty(n, dtype=np.int64)
    p = 0
    for i in range(n):
        v = 0
        s = 0
        while True:
            b = buf[p]
            p += 1
            v |= (b & 0x7F) << s
            if not b & 0x80:
                break
            s += 7
        out[i] = v
    return out, p


def zz_enc(vals):
    """LEB128 of (v << 1) ^ (v >> 63)."""
    return u_enc([((int(v) << 1) ^ (int(v) >> 63)) & ((1 << 64) - 1) for v in vals])


def zz_dec(buf, n):
    v, p = u_dec(buf, n)
    return (v >> 1) ^ -(v & 1), p


def bitmap_enc(bits):
    """One bit per run, LSB-first, zero-padded to a byte."""
    return np.packbits(np.asarray(bits, dtype=np.uint8), bitorder="little").tobytes()


def bitmap_dec(buf, n):
    return np.unpackbits(np.frombuffer(buf, dtype=np.uint8), count=n, bitorder="little")


# --- strand-aware axis-b delta (SPEC 1.2, `delta prev_exit`) ----------------

def b_enter_exit(b, ln, rev):
    """Where a run enters and leaves axis b, following its traversal direction.

    Forward: enter = b,       exit = b + len
    Reverse: enter = b + len, exit = b
    """
    enter = np.where(rev == 0, b, b + ln)
    exit_ = np.where(rev == 0, b + ln, b)
    return enter, exit_


# --- chunking (SPEC 1.3) ----------------------------------------------------

def chunk_bounds(a, b, ln, amem, bmem, count, bspan, seqbound):
    """Close a chunk on run count, axis-b span, or a member crossing on a named axis.

    `seqbound` must be 'both' when both axes are seqlocal: a chunk holds runs for exactly one
    member per axis, which is what makes its b-extent a single comparable interval and hence
    what makes the reverse skip work at all. A chunk always takes at least one run, so a run
    longer than `bspan` may exceed it.
    """
    n = len(a)
    bounds = [0]
    held = 0
    lo = hi = 0
    for i in range(n):
        if held == 0:
            lo, hi = b[i], b[i] + ln[i]
            held = 1
            continue
        nlo, nhi = min(lo, b[i]), max(hi, b[i] + ln[i])
        crossed = ((seqbound in ("a", "both") and amem[i] != amem[i - 1]) or
                   (seqbound in ("b", "both") and bmem[i] != bmem[i - 1]))
        if held >= count or (bspan and (nhi - nlo) > bspan) or crossed:
            bounds.append(i)
            lo, hi = b[i], b[i] + ln[i]
            held = 1
        else:
            lo, hi = nlo, nhi
            held += 1
    bounds.append(n)
    return np.array(bounds)


# --- run codec (SPEC 1.2) ---------------------------------------------------

def encode_chunk(a, b, ln, rev, a_enc="uvarint"):
    """Return the chunk's streams, in declaration order: a, b, len, strand."""
    prev_a_exit = np.concatenate(([a[0]], (a + ln)[:-1]))
    ga = a - prev_a_exit

    enter, exit_ = b_enter_exit(b, ln, rev)
    prev_b_exit = np.concatenate(([enter[0]], exit_[:-1]))
    gb = enter - prev_b_exit

    streams = [u_enc(ga) if a_enc == "uvarint" else zz_enc(ga),
               zz_enc(gb),
               u_enc(ln - 1),
               bitmap_enc(rev)]
    return streams


def decode_chunk(streams, n, a0, b_enter0, a_enc="uvarint"):
    """Inverse of encode_chunk.

    `a0` and `b_enter0` are the chunk's bases and come from the chunk directory: every delta
    stream restarts at zero at a chunk boundary, which is what makes chunks independently
    decodable.
    """
    ga, _ = (u_dec(streams[0], n) if a_enc == "uvarint" else zz_dec(streams[0], n))
    gb, _ = zz_dec(streams[1], n)
    lm1, _ = u_dec(streams[2], n)
    ln = lm1 + 1
    rev = bitmap_dec(streams[3], n).astype(np.int64)

    a = np.empty(n, np.int64)
    b = np.empty(n, np.int64)
    pa = a0
    prev_exit = b_enter0
    for i in range(n):
        a[i] = pa + ga[i]
        pa = a[i] + ln[i]
        enter = prev_exit + gb[i]
        b[i] = enter if rev[i] == 0 else enter - ln[i]
        prev_exit = (b[i] + ln[i]) if rev[i] == 0 else b[i]

    return a, b, ln, rev


# --- block compression (SPEC 1.2: each stream compressed independently) -----

def compress_streams(streams, level=9):
    return [zlib.compress(s, level) for s in streams]


def decompress_streams(blobs):
    return [zlib.decompress(s) for s in blobs]
