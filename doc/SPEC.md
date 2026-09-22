# `libintervalmap` — format schema (draft 0.1)

An **interval map** is a set of colinear ungapped *runs* relating intervals on axis **a** to
intervals on axis **b**, chunked, indexed from both axes, in one self-describing file.

Two profiles are defined here: `hal2.edge` (a child genome onto its parent) and `taffy.tui`
(a genome onto a universal column axis). They differ only in the schema text — the container,
codec, chunking and index machinery are identical.

Design rules this draft is required to honour:

- **Every stored offset is relative to the file's own byte 0.** No `SEEK_END`, no absolute
  paths. This is what lets a file be embedded in a zip/tar/pack later without a format change.
- **All I/O through `pread(handle, buf, off, len)`** where `handle = {backend, base, len}`.
  Backends: local fd, `udc2` (HTTP range), container slice, memory.
- **Open is two `pread`s** — the trailer, then the footer. No temp files, no eager per-type
  index arrays, no scan.
- **The schema describes every byte, payload included.** A reader that understands this
  document plus the embedded schema text can decode every run without knowing what the axes
  mean. (This is the thing ONEcode stops short of: its `D R 2 3 INT 6 STRING` declares "an int
  and a blob", which is ~1% of a `.tui`.)

---

## 1. Schema text

Plain ASCII, stored verbatim in the footer, one directive per line, `#` to end of line is a
comment. It is the normative description of the file: **if the schema and the writer disagree,
the writer is wrong.**

### 1.1 Directives

```
imap <major>                          format version of this document
profile <name> <version>              who wrote it and what the axes mean

axis <a|b> <label> coord <seqlocal|global> [extent <n>]
                                      seqlocal: coordinates are (member, offset)
                                      global:   one integer axis, `extent` = its cardinality

field <name> delta <none|prev_end> enc <uvarint|zigzag|bitmap> [min <n>]
                                      one declaration per stream, in stream order

order <a|b>                           run order within a chunk
chunk count <n> [bspan <n>] [seqbound <a|b|both|none>]
block <none|zstd|zlib>

section <id> <kind> [args...]         one per section present
```

### 1.2 Field semantics (normative)

A run is `(a, b, len, strand)`: `[a, a+len)` on axis **a** corresponds to `[b, b+len)` on
axis **b**, in the orientation given by `strand` (`0` = same, `1` = opposite). On
`strand = 1`, `a + i` pairs with `b + len-1 - i`.

| field | stored as | delta base |
|---|---|---|
| `a` | `a - prev_a_exit` | `prev_end` — previous run's `a + len` |
| `b` | `b_enter - prev_b_exit` | `prev_exit` — **strand-aware, see below** |
| `len` | `len - 1` (a run is ≥ 1 long) | `none` |
| `strand` | one bit | `none` |

A chunk holds runs for exactly one member on each axis (see 1.3), so member ids are carried
once per chunk in the directory, never per run.

**`delta prev_exit` — the strand-aware b delta.** Define, per run:

```
enter = (strand == 0) ? b       : b + len
exit  = (strand == 0) ? b + len : b
stored_b[i] = enter[i] - exit[i-1]          exit[-1] := enter[0], so stored_b[0] == 0
```

Decoding needs `len` and `strand` first, then `b = enter` if forward else `enter - len`.
Streams are separately length-delimited, so this imposes no ordering constraint on the file —
only on the reconstruction step.

The point is that a reverse-strand stretch walks *backwards* through axis b. A naive
`b - (prev_b + prev_len)` charges a large negative jump for every run in such a stretch;
following the traversal direction keeps the delta near zero instead. Measured on the fish
edge, 48.4% of whose runs are reverse: the b stream drops from 0.838 to 0.497 B/run, and the
whole file from 2.369 to **1.984 B/run — a 16% saving**. `.tui` does not do this today; it
instead sorts runs by b within a chunk to get small b deltas, which achieves something
similar by a more expensive route (and costs it the non-negative a delta). **Whether tui can
drop its g-sort once it has this delta is an open question worth testing** — if so it loses
a sort and gains `uvarint` on axis a.

Encodings: `uvarint` = LEB128 of a non-negative value. `zigzag` = LEB128 of
`(v << 1) ^ (v >> 63)`. `bitmap` = one bit per run, LSB-first, zero-padded to a byte.

`strand` gets its own stream rather than being packed into `len` (as `.tui` does with
`len<<1|strand`): packing destroys the length stream's delta structure. Measured on real HAL
edges the split is 19–22% smaller overall while the strand bitmap costs 0.002 B/run.

**Each stream is compressed independently**, not concatenated then compressed. Measured gain
is modest but free and grows as chunks shrink: −0.6% at 65536 runs, −4.1% at 8192. Stream
order has no measurable effect on size, so it is fixed by declaration order for determinism.

### 1.3 Chunking (normative)

A chunk is closed when **any** of these trips:

- `count` runs have accumulated;
- adding the run would make the chunk's span on axis **b** exceed `bspan`;
- the run crosses to a different member on an axis named by `seqbound`.

**`seqbound both` is required when both axes are `seqlocal`.** A chunk must hold runs for
exactly one member on each axis. This is not a size optimisation — it is what makes a chunk's
`b_min`/`b_span` a single well-defined interval, and therefore what makes the reverse skip
possible at all. Coordinates from different members are not comparable, so a chunk spanning
several axis-b members has no meaningful b-extent to test a query against.

Cutting at axis-a boundaries additionally makes axis-a deltas provably non-negative within a
chunk, which is what licenses `uvarint` there: measured, 0 negative a-deltas with the cut,
93 without.

**Rejected alternative, recorded so it is not re-litigated:** cut on axis a only and carry
the axis-b member as a per-run delta stream. It is 4% *smaller* (1.984 vs 2.064 B/run) because
the parent here is a Cactus ancestor with 5,575 scaffolds and this avoids shattering the
payload into ~10.8k chunks. It is also unusable: a chunk then spans up to **160 axis-b
members**, and reverse point lookup degrades from a mean of 1.23 chunks to **318** (p99 621).
Size is not the thing to optimise here.

**`bspan` is load-bearing, not a tuning knob.** Axis-a lookups always touch exactly one chunk
(a-intervals are disjoint and sorted), but axis-b lookups fan out over every chunk whose
b-range covers the query. Measured, reverse point lookup:

| edge | chunking | mean chunks touched | p99 |
|---|---|---|---|
| ape internal | `count 65536` | 49.6 | 50 |
| ape internal | `count 65536 bspan 1e6` | **1.99** | 30 |
| fish (deep divergence) | `count 65536` | 113.4 | 132 |
| fish (deep divergence) | `count 65536 bspan 1e6` | **1.48** | 12 |

A writer that ignores `bspan` produces a readable file with catastrophic reverse-lookup cost.
A single run longer than `bspan` is allowed to exceed it — a chunk always takes at least one
run.

`count` matters less than expected on fragmented genomes: with `seqbound both` against a
5,575-scaffold parent, member boundaries dominate the cutting, so 65536 and 8192 produce
10,839 and 11,218 chunks and 2.064 vs 2.073 B/run — a 0.4% difference. It still matters on
genomes with few, long members, where `count` is the only thing closing a chunk.
**Default 8192**, which bounds bytes inflated per random probe without measurable cost.

`order` is declared, not assumed. `order a` is the default and keeps axis-a deltas
non-negative; `order b` is what `.tui` does today.

**Implementation status.** The C writer implements `order a` and `seqbound both` only,
and encodes axis a as `uvarint`; it rejects a chunk whose extent would not fit the u32
span fields rather than letting one wrap. A file declaring `order b` or
`field a ... enc zigzag` is legal per this document but is not yet produced or read.

---

## 2. Container

Little-endian throughout. Sections are 8-byte aligned.

```
HEADER  (64 B, at offset 0 — sized so format sniffers that read 64 bytes see all of it)
  0   u8[8]   magic  = "IMAP\x1A\x0D\x0A\x00"
  8   u32     imap_major
  12  u32     imap_minor
  16  u64     feature_flags      unknown bit set => reader MUST refuse
  24  u8[32]  profile            NUL-padded, e.g. "hal2.edge"
  56  u64     reserved (0)

...sections...

FOOTER  (variable): section table, member directories, chunk directories,
                    schema text, provenance, summary stats

TRAILER (32 B, last 32 bytes of the file)
  0   u64     footer_off         relative to byte 0
  8   u64     footer_len
  16  u32     footer_crc32          CRC-32 (IEEE); computed in bounded steps,
                                 since zlib's uInt truncates a >4 GiB length
  20  u32     reserved (0)
  24  u8[8]   magic (repeated)
```

Open = `pread` the last 32 B, then `pread` the footer. Nothing else is read until queried.

### 2.1 Section table

```
  u32 n_sections
  repeat:
    u8   id_len;  char id[id_len]      e.g. "runs", "dir.a", "mem.b"
    u64  off                           relative to byte 0
    u64  len                           bytes on disk
    u64  raw_len                       uncompressed
    u8   codec                         0 none, 1 zlib  (per-section; see below)
    u8   flags
    u32  crc32
```

Unknown section ids are skipped, not an error — that plus `feature_flags` is the extension
mechanism. **Reserved for future use: section ids `aux.0`–`aux.3` and the `run` field names
`aux0`–`aux3`.** Retrofitting a stream into a deployed `.tui` cost 27 hours and 11.8 TiB of
scratch; the cost of reserving them now is these two lines.

### 2.2 Chunk directory `dir.a` — fixed width, uncompressed, binary-searchable

**64 bytes** per chunk, sorted by `(a_member, a_min)`, little-endian, in this order:

```
   0  u32 a_member       4  u32 b_member
   8  u64 a_min         16  u32 a_span     20  u32 b_span
  24  u64 b_min
  32  u64 b_enter0      <- the first run's axis-b entry point; a decode base
  40  u64 off           <- relative to the start of the `runs` section
  48  u32 clen          52  u32 rawlen     56  u32 n_runs     60  u32 codec
```

`b_enter0` is **not** `b_min`. The deltas begin at the first run, and under `order a` that
is not the run with the smallest b. `a_min` doubles as the axis-a decode base, which holds
only because `order a` makes the first run the smallest — under `order b` it would not.

`codec` is per chunk and not inherited from the section: the `runs` section itself is stored
uncompressed and each chunk blob inside it is compressed individually, so the reader has to
be told which. Inferring it from `rawlen != clen` is wrong for a chunk that happens not to
compress.

A chunk blob is `u32 len[4]` — the four stream lengths in declaration order — followed by the
four streams; `codec` applies to that whole blob.

`a_span`/`b_span` are the *tight* extent of the chunk's runs on each axis, and `b_span` is
what lets a reverse query skip a chunk without decompressing it. **Tightness is a validated
invariant, not a convention** — every skip is silently wrong if a range is loose, and `.tui`
never asserted it. Both are u32: a writer must close a chunk before its extent would exceed
2^32 rather than let the field wrap, which is the same failure in a quieter form.

Left uncompressed so it can be `pread` and binary-searched in place. A 1.3 Gb genome at 8192
runs/chunk is ~3,200 chunks ≈ 205 KB.

### 2.3 Chunk directory `dir.b` — a permutation, not a copy

12 bytes per chunk, sorted by `(b_member, b_min)`:

```
  u32 chunk_id                  index into dir.a
  u64 prefix_max_b_end          running max of (b_min + b_span) over dir.b[0..i]
```

`prefix_max_b_end` is the early-exit bound for the backward walk over overlapping chunks.
Storing a permutation rather than a second copy of the runs was measured: a second
parent-sorted payload costs +2.94 B/run and buys fan-out 1.00 against the cap's 1.48. Not
worth it.

### 2.4 Member directory `mem.a` / `mem.b`

```
  u32 n_members
  repeat: u32 name_off  u32 name_len  u64 length
          u64 axis_min  u64 axis_max          extents on this axis
          u32 first_chunk  u32 n_chunks
  name blob
  u32 by_name[n_members]                      member ids, sorted by name
```

Extents live here on purpose: omitting them is what costs `.tui` 9 s per query at 577-way
with no cheap workaround. `by_name` gives O(log n) name lookup with no build-time structure;
a persisted CHD perfect hash is a later optimisation, not a format change.

Because coordinates are `seqlocal`, **there is no global position → member index anywhere.**
That whole structure — and HAL's `getSequenceBySite`, which `DnaIterator::getSequence()` calls
per base per column and which degrades to an allocating binary search above 1000 sequences
(`hdf5Genome.cpp:37`, `:360-378`) — simply does not exist here.

### 2.5 Optional sections

| id | kind | for |
|---|---|---|
| `anchor` | sparse `key → value` map, both monotone, delta+varint | `.tui`'s column → MAF file offset |
| `group` | `(cluster_id, member, start, len)` records | hal2 paralogy as equivalence classes |
| `stats` | counts, coverage, run-length and gap histograms | O(1) `halStats` |
| `prov` | program, version, command, date — one record per writer | provenance |
| `meta` | free-form `key\tvalue` text | e.g. the `# hal` Newick |

---

## 3. The two profiles

### 3.1 `hal2.edge` — a child genome onto its parent

```
imap 1
profile hal2.edge 1

axis  a  child   coord seqlocal
axis  b  parent  coord seqlocal

field a       delta prev_end   enc uvarint
field b       delta prev_exit  enc zigzag
field len     delta none       enc uvarint  min 1
field strand  delta none       enc bitmap

order   a
chunk   count 8192  bspan 1000000  seqbound both
block   zstd

section mem.a  memdir  axis a
section mem.b  memdir  axis b
section runs   runs
section dir.a  chunkdir   order a
section dir.b  chunkperm  dir.a  order b  prefixmax b_end
section stats  stats
section prov   prov
```

Axis a is strictly sorted and disjoint (a child base has 0 or 1 parent), so `a` is `uvarint`.
Axis b may repeat and reorder — that is paralogy and rearrangement — so `b` is `zigzag`.
**Paralogy is not stored.** The copies of the child at parent position *p* are exactly what a
`dir.b` query at *p* returns; `--noDupes` picks the smallest `a` among them at query time.
Measured against HAL's stored canonical on a real fish edge: agreement on 100.000% of 876,531
duplicated positions, zero differences.

### 3.2 `taffy.tui` — a genome onto the universal column axis

```
imap 1
profile taffy.tui 3

axis  a  genome  coord seqlocal
axis  b  column  coord global  extent 72489835721

field a       delta prev_end   enc zigzag
field b       delta prev_exit  enc zigzag
field len     delta none       enc uvarint  min 1
field strand  delta none       enc bitmap

order   b
chunk   count 65536  bspan 1000000  seqbound a
block   zstd

section mem.a   memdir  axis a
section runs    runs
section dir.a   chunkdir   order a
section dir.b   chunkperm  dir.a  order b  prefixmax b_end
section anchor  anchor  key b  value fileoff  every 10000
section meta    meta
section prov    prov
```

No `mem.b` — axis b is `global`, so it has one implicit member of size `extent`. `order b`
and `count 65536` preserve today's `.tui` behaviour; `seqbound a` only, since axis b has no
members to bound against. The `# hal` Newick and `max_gap` live in `meta`.

That is the whole difference between the two formats: **axis declarations, run order, chunk
constants, and which optional sections are present.**

---

## 4. Status and open questions

**Validated by reference implementation** (`ref.py`, round-tripped against 9,010,509 real
runs from `GCA_024256435.1 -> SiluriformesAnc2`, every 211th chunk decoded and compared
field-by-field):

| configuration | chunks | B/run (zlib-9) | reverse fan-out (mean / p99) | round-trip |
|---|---|---|---|---|
| `count 8192 bspan 1e6 seqbound both` | 11,218 | **2.073** | **1.24 / 8** | OK |
| `count 65536 bspan 1e6 seqbound both` | 10,839 | 2.064 | 1.23 / 8 | OK |
| *(rejected)* `seqbound a` + bmem stream | 656 | 1.984 | 318 / 621 | OK |

All with sequence-local coordinates on both axes, the strand-aware b delta, `len-1`, a strand
bitmap, and per-stream compression. zstd at the same settings measured 5–7% below zlib-9 on
this data in earlier runs, so ~1.93–1.97 B/run is the expected figure.

Per-field round-trip (a, b, len, strand) was checked on every 97th chunk: 0 mismatches.

Still open, to be measured rather than argued:

1. **`order a` vs `order b` for `hal2.edge`.** `.tui` measured −18% from b-ordering on apes,
   but that predates the strand-aware delta, which may subsume most of the benefit. Test on a
   deep and a shallow edge.
2. **Whether `.tui` can drop its within-chunk g-sort** if it adopts `delta prev_exit`. Same
   experiment, run against tui's own data.
3. **zstd level, and whether a per-section trained dictionary earns its keep.** A dictionary
   measured 2.55 vs 2.57 B/run at 65536 runs but 2.82 vs 3.13 at 256 — it matters only if
   chunks end up small.
4. **`b` as `uvarint` under `order b`.** If b-ordered chunks make the enter-delta monotone it
   can drop the zigzag bit.
