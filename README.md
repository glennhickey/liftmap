# liftmap

A small C library for **indexed pairwise coordinate maps between genomes**: a set of
colinear ungapped *runs*, each relating an interval of a sequence on one axis to an interval
of a sequence on another, chunked, compressed, and indexed so that lifting coordinates from
*either* side is fast.

```
run = (seqA, a, seqB, b, len, strand)      [a, a+len) on seqA  <->  [b, b+len) on seqB
```

Lifting from axis a returns at most one hit per base (axis a is a function, unless a file
declares otherwise); lifting from axis b can return several — paralogy is reported, not
resolved.

## Users

- **taffy** — the `.tui` universal column index (a genome onto a universal MAF column axis)
  is a liftmap file, profile `taffy.tui`.
- **hal2** — the planned modular HAL backend stores each tree node's edge (child genome onto
  parent) as a liftmap file, profile `hal2.edge`.

## Status

The C library implements the format in `doc/SPEC.md` (format major 2): writer, reader,
a cursor for key-ordered and point access from either axis, metadata, groups (e.g. a
genome's sequences), application sections, and a `pread` I/O seam over files, memory and
container slices. `ref/liftmap.py` is a Python reference for the run codec and chunking.

liftmap began as `libintervalmap`; its history is carried over.

## What is validated

- Round trip and both query directions against an oracle on 9,010,509 runs from a real HAL
  edge (`GCA_024256435.1` → `SiluriformesAnc2`, 577-way VGP fish alignment): 0 mismatches,
  ~2.07 bytes per run.
- taffy's `.tui` on this library matches the previous ONEcode implementation byte for byte
  on `taffy view -U` and `taffy lift` over rodent and fish-subtree universal alignments, in
  2–3% less space (`doc/SPEC.md` §4).
- Corruption: header, footer, directories and every chunk are CRC-checked; a CRC-repairing
  fuzzer (`tests/crcfuzz.py`) runs mutants under ASan/UBSan.

## Design constraints

1. Every stored offset is relative to the file's own byte 0, so a file can be embedded in a
   container without a format change.
2. All I/O goes through one `pread(handle, buf, off, len)` seam. Local file, HTTP range,
   container slice and memory are the same code path.
3. Open reads the trailer, footer and directories — never the payload.
4. Both axes are indexed from the same payload bytes.
5. Extension slots are reserved up front (`aux.0`–`aux.3`).

## Command line

```
liftmap from-paf   in.paf[.gz]   out.lmap   [--swap] [--allow-overlap] [--group-sep C [--group-fields N]]
liftmap from-chain in.chain[.gz] out.lmap   [--swap] [--allow-overlap] [--group-sep C [--group-fields N]]
liftmap to-paf     in.lmap [out.paf]        [--max-gap N]
liftmap to-chain   in.lmap [out.chain]      [--max-gap N]
liftmap lift       in.lmap in.bed [out.bed] [--from a|b]
liftmap dump | info | verify  in.lmap
```

Axis a is the sequence named first in each record — the PAF query, the chain target — so
a → b is each format's own lift direction. PAF needs `cg:Z:` CIGARs.

Import keeps the aligned base pairs and their strands, as the unique set of maximal runs:
the same base pairs give the same file however they were split into records, and duplicate
records collapse. Scores, mapping qualities and record boundaries are not kept; export
regroups runs into records across gaps of at most `--max-gap` bp and writes placeholder
scores. If a base of axis a aligns more than once (secondary alignments, or `--swap` onto
a sequence with several copies), import refuses unless `--allow-overlap`. `--group-sep`
groups sequences into genomes by the name up to the N-th separator: `.` 1 for `hg38.chr1`,
`.` 2 for `GCA_000001635.9.chr1`, `#` 2 for PanSN `HG002#1#chr1`.

Checked against independent tools on the evolver mammals HAL (8 branches): `lift` agrees
with `halLiftover` on 294,220 single-base mappings in both directions, and with UCSC
`liftOver` (via `to-chain`, and `chainSwap` for the reverse) on 73,604, strands included.

## Building

```
make            # libliftmap.a, bin/liftmap and the test programs
make check      # self-contained: generated data, brute-force oracle
```

## Licence

TBD.
