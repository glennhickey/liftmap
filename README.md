# libintervalmap

A small C library for **indexed pairwise interval maps**: a set of colinear ungapped *runs*
relating intervals on one coordinate axis to intervals on another, chunked, compressed, and
indexed so that lookups from *either* axis are fast.

It is a generalization of the working core of taffy's `.tui` universal column index, extracted
so that more than one format can sit on it. Two profiles are specified:

- **`hal2.edge`** — a child genome's coordinates onto its parent's, one file per node of a
  phylogenetic tree. The intended storage layer for a modular replacement of HAL's HDF5
  backend.
- **`taffy.tui`** — a genome's coordinates onto a universal MAF column axis. What `.tui`
  does today, re-expressed over this library.

The two differ only in axis declarations, run order, chunk constants, and which optional
sections are present. Everything else — codec, chunking, directories, container — is shared.

## Status

**Draft.** `doc/SPEC.md` is the normative specification. `ref/imap.py` is a Python reference
implementation of the run codec and chunking, used to validate the spec, not to be fast.

The C library does not exist yet.

## What is validated

The codec and chunking in `ref/imap.py` have been round-tripped field-by-field against
9,010,509 real runs taken from a HAL edge at deep vertebrate divergence
(`GCA_024256435.1` → `SiluriformesAnc2`, from a 577-way VGP fish alignment): sequence-local
coordinates on both axes, strand-aware b deltas, `len-1`, a strand bitmap, and per-stream
compression. Result: **2.073 bytes per run** at `count 8192 bspan 1e6 seqbound both`, with a
reverse point-lookup fan-out of 1.24 chunks on average (p99 8). Member purity was asserted on
all 11,218 chunks and per-field round-trip checked on a sample: 0 mismatches.

A smaller variant (1.984 B/run) that carries the axis-b member as a per-run stream instead of
cutting chunks at member boundaries was measured and **rejected**: it lets a chunk span 160
parent scaffolds, which destroys the reverse skip (fan-out 318, p99 621). The negative result
is recorded in `doc/SPEC.md` §1.3 so it is not re-litigated.

## Design constraints that are not negotiable

These are the ones where a wrong choice cannot be fixed later:

1. Every stored offset is relative to the file's own byte 0 — so a file can be embedded in a
   container without a format change.
2. All I/O goes through one `pread(handle, buf, off, len)` seam, where the handle carries a
   base offset and a length. Local file, HTTP range, container slice and memory are all the
   same code path. Never a `FILE*`.
3. Open costs two `pread`s — a fixed trailer, then the footer. No temp files, no eager
   per-type index arrays, no scan.
4. Both axes are indexed from the same payload bytes.
5. Extension slots are reserved up front (`aux.0`–`aux.3`), because retrofitting a stream
   into a deployed format is a re-index of every file.

## Licence

TBD.
