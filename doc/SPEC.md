# `liftmap` — format specification

A **liftmap** is a set of colinear ungapped *runs* relating intervals on axis **a** to
intervals on axis **b**, chunked, indexed from both axes, in one self-describing file.

Two profiles are defined here: `hal2.edge` (a child genome onto its parent) and `taffy.tui`
(a genome onto a universal column axis). They differ only in run order, chunk constants,
metadata and application sections — the container, codec, chunking and index machinery are
identical.

Design rules this draft is required to honour:

- **Every stored offset is relative to the file's own byte 0.** No `SEEK_END`, no absolute
  paths. This is what lets a file be embedded in a zip/tar/pack later without a format change.
- **All I/O through `pread(handle, buf, off, len)`** where `handle = {backend, base, len}`.
  Backends: local fd, HTTP(S) range (liftmap's libcurl backend, §5), container slice,
  memory, or any backend an application supplies.
- **Open is two `pread`s** — the trailer, then the footer. No temp files, no eager per-type
  index arrays, no scan.
- **This document describes every byte, payload included.** A reader that implements it can
  decode every run without knowing what the axes mean. (This is the thing ONEcode stops short
  of: its `D R 2 3 INT 6 STRING` declares "an int and a blob", which is ~1% of a `.tui`.)

**History.** liftmap began as `libintervalmap` (prefix `imap_`, magic `IMAP`); its format 4
is liftmap's format 1, renamed. Measurements below that mention `.lmap` files were taken
under the old name. Format 2 replaced a free-text "schema" section with metadata (2.5) and
added groups (2.4a). Format 3 added per-page directory checksums (`dir.crc`, 2.3a) so a
reader loads directories lazily.

---

## 1. Runs

The run model and its codec are fixed by this document and the format major version. An
earlier draft embedded a schema text (`axis`, `field`, `order`, `chunk` directives) meant to
make them configurable per file; nothing ever read it, so format 2 dropped it. What varies
per file is recorded where it is used: run order and the axis-a overlap bit in the header
(2), chunk constants implicitly in the chunk directory, and what the axes mean in metadata
(2.5).

### 1.2 Field semantics (normative)

A run is `(a, b, len, strand)`: `[a, a+len)` on axis **a** corresponds to `[b, b+len)` on
axis **b**, in the orientation given by `strand` (`0` = same, `1` = opposite). On
`strand = 1`, `a + i` pairs with `b + len-1 - i`.
Every run lies inside both of its members: `a + len ≤` the axis-a member's length and
`b + len ≤` the axis-b member's. The writer refuses a run that does not.

**Axis a is a partial function** — each position of an axis-a member maps to at most one
position on axis b — **unless the header sets feature bit 1 (`A_OVERLAP`)**. That bit exists
for maps that are approximate by construction: a chained `.tui` merges runs across gaps, so a
run's length is its span on axis b and no longer its span on axis a, and runs of one sequence
routinely overlap on it. On the rodent universal index, forcing such a map back to a partial
function by trimming dropped 3.8M of its 7.8M runs outright. `A_OVERLAP` requires order b,
because order a stores `a` as an unsigned delta from the previous run's end, which cannot go
backwards. Axis b is a multimap either way.

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
whole file from 2.369 to **1.984 B/run — a 16% saving**.

**The b delta base must match the run order.** `prev_exit` follows the traversal, which is
right when runs are in axis-a order. Under `order b` the runs are sorted by b-start, so the
naive `prev_end` delta is already small and mostly non-negative for both strands, and
following the traversal instead *breaks* that for reverse runs. Measured on the evolver
universal-column index, order b: `prev_end` 2.146 B/run, `prev_exit` 2.285. Hence:
`order a` ⇒ `b delta prev_exit`; `order b` ⇒ `b delta prev_end`. This also settles an earlier
open question: the strand-aware delta does not let `.tui` drop its within-chunk g-sort.

Encodings: `uvarint` = LEB128 of a non-negative value. `zigzag` = LEB128 of
`(v << 1) ^ (v >> 63)`. `bitmap` = one bit per run, LSB-first, zero-padded to a byte.

`strand` gets its own stream rather than being packed into `len` (as `.tui` does with
`len<<1|strand`): packing destroys the length stream's delta structure. Measured on real HAL
edges the split is 19–22% smaller overall while the strand bitmap costs 0.002 B/run.

**Each stream is compressed independently**, not concatenated then compressed, using raw
deflate (no zlib wrapper: four streams per chunk and ~10^4 chunks per file make the 6-byte
header and Adler checksum add up, and each chunk's directory CRC covers integrity). A stream is
deflated only where that makes it smaller; otherwise it is stored verbatim. Measured in the C
implementation against compressing each chunk as one blob: −1.5% on the fish HAL edge in
order a, −2.9% in order b, −1.7% on the evolver `.tui` data. Stream order has no measurable
effect on size, so it is fixed by declaration order for determinism.

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

`order` is declared, not assumed, and **the right choice depends on the data** — measured,
per-stream zlib-9, `seqbound both`:

| data | order a, `prev_exit` | order b, `prev_end` |
|---|---|---|
| HAL edge, fish leaf→ancestor (deep, 48% reverse) | **2.074** | 2.511 |
| HAL edge, ape internal (shallow) | **2.938** | 3.517 |
| universal columns, evolver `.tui` (123,471 runs) | 3.097 | **2.098** |

A HAL child→parent map is mostly colinear, so axis-a order keeps b nearly monotone as well,
and the strand-aware delta absorbs inversions; sorting by b scrambles a instead. On a
universal column axis one sequence's consecutive runs are separated by every other
lineage's novel columns, so a-order produces large b jumps that b-order avoids. Hence
`hal2.edge` is `order a` and `taffy.tui` is `order b`.

**Chunks are cut along the order axis.** Under `order a` runs arrive in axis-a order and are
cut as they come. Under `order b` a whole axis-a member is buffered, sorted by
`(b_member, b, a)`, and cut in that order — which is what `.tui` does (`tui.c` phase 2: a
sequence's runs are colinear-merged in `t` order, the whole sequence is sorted by `g`, then
cut). Axis a is then `zigzag`, and a member's chunks overlap on axis a, which `dir.a.max`
(section 2.2) handles.

An earlier draft said order b cuts in axis-a order and sorts only *within* each chunk. That
is not equivalent, and it was wrong: on a universal column axis a sequence's consecutive runs
jump across columns at every rearrangement, so the `bspan` cap fired almost every run. On the
rodent `.tui` (85M runs, T = 2.73 × 10^9) that produced 44.6M chunks — under two runs each —
and a file 22× the size of the `.tui`. Cutting along b gives 36,709 chunks and a file 4%
*smaller* than the `.tui`.

`bspan` caps a chunk's span on axis b in both orders, but plays a different role in each.
Under order a it bounds the *non-cut* axis, which is what keeps reverse-lookup fan-out near
one. Under order b it bounds the cut axis itself, so one chunk never covers more than `bspan`
columns — the role of `.tui`'s `TUI_CHUNK_G_MAX`, which cut a full-chromosome lift to a
distant target from 196 s to 1.7 s.

**Coarsening.** A coarsened file (`lmap_coarsen`) is an ordinary liftmap whose runs were
chained across gaps: along each member of a key axis, a member of the other axis extends its
open chain when the strand matches and both gaps are within `max_gap`, and the chain's
length becomes its key-axis span. Its extent on the other axis is therefore approximate, it
may overlap there (`A_OVERLAP`), and it is clipped at the member's end. `max_gap` in the
metadata says so. Readers need nothing special; exact-coordinate consumers should refuse a
file whose `max_gap` is non-zero unless told otherwise, as taffy does.

**Implementation status.** The C library implements both orders with `seqbound both`,
recorded in header feature bit 0 (section 2), and `A_OVERLAP` (bit 1) under order b. It rejects a chunk whose extent would not fit
the u32 span fields rather than letting one wrap. Under order b the writer holds one axis-a
member's runs in memory at a time, which is the same bound `.tui`'s own builder has.

**The writer streams the payload.** The `runs` section is written straight after the header
as each chunk is encoded, and every other section follows it at close (the header is
written last, over a placeholder); section order is free, since the footer locates
everything. What the writer holds is the chunk directory (80 bytes per chunk: ~0.8 GB at
10^7 chunks) and, under order b, one axis-a member's runs while they are sorted and cut.
An earlier version buffered the whole compressed payload until close -- ~100 GB for a
577-way `.tui`.

---

## 2. Container

Little-endian throughout. Sections are 8-byte aligned.

```
HEADER  (64 B, at offset 0 — sized so format sniffers that read 64 bytes see all of it)
  0   u8[8]   magic  = "LMAP\x1A\x0D\x0A\x00"
  8   u32     lmap_major         3; a reader requires an exact match
  12  u32     lmap_minor         0
  16  u64     feature_flags      unknown bit set => reader MUST refuse
                                 bit 0: runs within a chunk are in order b
                                 bit 1: A_OVERLAP -- axis a may overlap (1.2);
                                        only with bit 0
  24  u8[32]  profile            NUL-padded, e.g. "hal2.edge"
  56  u64     reserved (0)

...sections...

FOOTER  (variable): the section table

TRAILER (32 B, last 32 bytes of the file)
  0   u64     footer_off         relative to byte 0
  8   u64     footer_len
  16  u32     footer_crc32          CRC-32 (IEEE); computed in bounded steps,
                                 since zlib's uInt truncates a >4 GiB length
  20  u32     header_crc32          CRC-32 of the 64-byte header
  24  u8[8]   magic (repeated)
```

**The header is checksummed, and verified before any field in it is trusted**, because
`feature_flags` decides how every chunk decodes. Before it was covered, flipping bit 0 on a
real order-a file made the reader open it as order b: the per-run containment checks
rejected 11,199 of its 11,218 chunks, and the other 19 decoded without complaint into
different runs.

Open = `pread` the last 32 B, then the footer, then the member directories, metadata,
groups and `dir.crc`. **Neither the `runs` section nor the chunk directories are read at
open**: both scale with the data, and an index is opened to touch a few members. Each chunk
carries its own CRC in its directory entry, checked whenever the chunk is read; each page
of the chunk directories has its own CRC in `dir.crc` (2.3a), checked when the page is
loaded. The whole-section CRCs in the section table are checked only by a full verification
pass (`lmap_verify`, which streams them in bounded reads).

History: reading and checksumming the whole payload at open cost 175 ms per open on the
228 MB rodent index; dropping that brought it to 7 ms. Reading the chunk directories still
cost 21 ms at 288,047 chunks and would reach hundreds of MB per open at 577-way scale;
paging them brought that to 0.2 ms.

### 2.1 Section table

```
  u32 n_sections
  repeat:
    u8   id_len;  char id[id_len]      e.g. "runs", "dir.a", "mem.b"
    u64  off                           relative to byte 0
    u64  len                           bytes on disk
    u64  raw_len                       uncompressed
    u8   codec                         0 none (sections are stored; chunks carry their own)
    u8   flags
    u32  crc32
```

Unknown section ids are skipped, not an error — that plus `feature_flags` is the extension
mechanism. **Reserved for future use: section ids `aux.0`–`aux.3` and the `run` field names
`aux0`–`aux3`.** Retrofitting a stream into a deployed `.tui` cost 27 hours and 11.8 TiB of
scratch; the cost of reserving them now is these two lines.

### 2.2 Chunk directory `dir.a` — fixed width, uncompressed, binary-searchable

**72 bytes** per chunk, sorted by `(a_member, a_min)`. Under `order a` a member's chunks are
also disjoint on axis a, since they were cut along it; under `order b` they generally
overlap. The writer sorts the directory at close (chunks are addressed by payload offset, so
reordering the directory moves nothing), and requires each axis-a member's runs to arrive in
non-decreasing `a` without overlap *across* chunks, not merely within one — axis a is a
partial function at the level of runs whatever the chunk layout (unless `A_OVERLAP`, 1.2). Little-endian, in this
order:

```
   0  u32 a_member       4  u32 b_member
   8  u64 a_min         16  u32 a_span     20  u32 b_span
  24  u64 b_min
  32  u64 base          <- decode base of the NON-ordering axis (below)
  40  u64 off           <- relative to the start of the `runs` section
  48  u32 clen          52  u32 rawlen     56  u32 n_runs     60  u32 codec
  64  u32 crc32         <- CRC-32 of the chunk's clen stored bytes
  68  u32 reserved (0)
```

Decoding a chunk needs the first run's position on both axes. The first run in stored order
is the minimum on the *ordering* axis, so that base is already in the entry; `base` holds the
other one:

| order | axis-a base | axis-b base |
|---|---|---|
| a | `a_min` | `base` = first run's b entry point (`b+len` if reverse, else `b`) |
| b | `base` = first run's a | `b_min` |

Neither stored base is a minimum: under order a the first run is not the one with the
smallest b, and under order b it is not the one with the smallest a.

`codec` is per chunk: `0` stores every stream verbatim, `1` raw-deflates each stream where
that makes it smaller. The `runs` section itself is stored uncompressed.

A chunk blob is

```
  u32 raw_len[4]       uncompressed length of each stream, in declaration order
  u32 stored_len[4]    bytes actually stored for it
  the four stored streams, back to back
```

`stored_len == raw_len` means the stream is stored verbatim; `stored_len < raw_len` means it
is raw deflate that must inflate to exactly `raw_len` bytes, consuming exactly
`stored_len`; `stored_len > raw_len` is invalid. There is no separate flag to disagree with
the lengths. Under `codec 0` every stream must be verbatim. In the directory,
`clen = 32 + Σ stored_len` and `rawlen = 16 + Σ raw_len`.

`a_span`/`b_span` are the *tight* extent of the chunk's runs on each axis, and `b_span` is
what lets a reverse query skip a chunk without decompressing it. **Tightness is a validated
invariant, not a convention** — every skip is silently wrong if a range is loose, and `.tui`
never asserted it. Both are u32: a writer must close a chunk before its extent would exceed
2^32 rather than let the field wrap, which is the same failure in a quieter form.

**Limits and consistency.** `n_runs` sizes a reader's decode buffer, so it must be
answerable to the bytes actually present, not merely declared:

- `1 <= n_runs <= 2^24`, and `rawlen <= 2^30`
- `rawlen >= 16 + 3*n_runs + ceil(n_runs/8)` — every run costs at least one varint byte
  in each of the a, b and len streams plus a strand bit, after the four-length header
- `clen >= 32`; `codec 0`: `clen == rawlen + 16`; `codec 1`: `clen <= rawlen + 16` and
  `rawlen <= 1032*clen + 64`, deflate's maximum expansion

Without these a 256-byte chunk could declare two billion runs and a reader would try to
allocate ~80 GB before discovering the streams are empty — a fuzzer found exactly that.

Every decoded run must lie inside its chunk's declared extents on both axes. That is what
makes skipping a chunk on its extent sound; a reader rejects a chunk that violates it.

Left uncompressed so it can be `pread` and binary-searched in place. A 1.3 Gb genome at 8192
runs/chunk is ~3,200 chunks ≈ 205 KB.

### 2.2a `dir.a.max` — the axis-a early-exit bound

One u64 per `dir.a` position: the running maximum of `a_min + a_span` over the entries at or
before it, **reset at each axis-a member** — the axis-a twin of `dir.b`'s prefix max. It is
what lets `query_a` walk backwards over chunks that overlap on axis a and stop as soon as
nothing earlier can reach the window. Under order a the member's chunks are disjoint and the
walk ends after one step. A reader recomputes it at open and refuses a mismatch, and under
order a additionally refuses overlapping extents within a member.

### 2.3 Chunk directory `dir.b` — a permutation, not a copy

12 bytes per chunk, sorted by `(b_member, b_min)`:

```
  u32 chunk_id                  index into dir.a
  u64 prefix_max_b_end          running max of (b_min + b_span) over dir.b[0..i]
```

`prefix_max_b_end` is the early-exit bound for the backward walk over overlapping chunks.
**It resets at each member.** Members are separate coordinate spaces, so a bound carried
over from the previous member would be meaningless; worse, it would be too large, which
does not produce wrong answers but silently disables the early exit for the next member.
A reader recomputes it (per member, on first use; 2.3a) and refuses a file whose stored
values differ.
Storing a permutation rather than a second copy of the runs was measured: a second
parent-sorted payload costs +2.94 B/run and buys fan-out 1.00 against the cap's 1.48. Not
worth it.

### 2.3a `dir.crc` — per-page checksums, and lazy loading

```
  u32 page_entries                       4096
  u32 reserved (0)
  u64 n_chunks
  u64 total_runs, payload_bytes          the directory's sums, so stats need no scan
  u32 crc[n_pages]  for dir.a            n_pages = ceil(n_chunks / page_entries)
  u32 crc[n_pages]  for dir.a.max
  u32 crc[n_pages]  for dir.b
```

`dir.a`, `dir.a.max` and `dir.b` are read a page at a time: the page is loaded, checked
against its CRC, and every entry in it checked on its own (2.2 limits, member ids, chunk ids
in range), the first time anything needs it. A shared reader installs pages atomically, so
concurrent readers stay safe. The invariants that span entries — a member's range carries
that member, is sorted, has the recomputed running max, and (for `dir.b`) names each chunk
once — are checked for a member's range the first time a cursor uses it, which is what the
early exits need. Together with the ranges tiling each directory exactly (checked at open),
the per-range checks imply the global ones: every chunk id in a `dir.b` range must carry
that range's member, so no id can appear in two ranges. `lmap_verify` checks every range.
A file whose directory is damaged therefore opens, and fails at the first query that
touches the damage.

### 2.4 Member directory `mem.a` / `mem.b`

```
  u32 n_members
  repeat: u32 name_off  u32 name_len  u64 length
          u64 axis_min  u64 axis_max          extents on this axis
          u32 first_chunk  u32 n_chunks
  name blob
  u32 by_name[n_members]                      member ids, sorted by name
```

`first_chunk`/`n_chunks` name a contiguous range **in a directory, not in the payload**:
for `mem.a` a range of `dir.a` positions, for `mem.b` a range of `dir.b` positions. A
member's chunks are contiguous in its own axis's directory and generally not in the
other's — a parent scaffold's chunks are scattered through `dir.a`. Every chunk belongs to
exactly one range per axis; a reader checks that the ranges tile the directory and that
each entry in a range carries that member's id.

Extents live here on purpose: omitting them is what costs `.tui` 9 s per query at 577-way
with no cheap workaround. `by_name` gives O(log n) name lookup with no build-time structure;
a persisted CHD perfect hash is a later optimisation, not a format change.

Because coordinates are `seqlocal`, **there is no global position → member index anywhere.**
That whole structure — and HAL's `getSequenceBySite`, which `DnaIterator::getSequence()` calls
per base per column and which degrades to an allocating binary search above 1000 sequences
(`hdf5Genome.cpp:37`, `:360-378`) — simply does not exist here.

### 2.4a Groups `grp.a` / `grp.b` (optional)

A group is a named set of members of one axis — typically one genome's sequences, where one
axis carries several genomes (`.tui`, multi-genome PAFs). A member is in at most one group.

```
  u32 n_groups                          >= 1; the section is absent when there are none
  u32 name_off[n_groups + 1]            name i is blob[name_off[i] .. name_off[i+1])
  name blob                             names unique, non-empty
  u32 group[n_members]                  group id per member, 0xFFFFFFFF = none
```

Group ids are declaration order, which is how a writer records an order it cares about
(`.tui` keeps its genome order). A reader derives each group's members in name order and
their summed length at open, so a genome's sequences, its size and its sequence count cost
nothing to ask for, and the member ids go straight to a cursor filter (2.6). Groups replace
the `"<genome>."` name-prefix convention, which cannot be parsed in general: `hg38.chr1`
splits at the first `.`, `GCA_000001635.9.chr1` at the second, PanSN `HG002#1#chr1` at the
second `#`. The writer records the answer instead.

Member names are unique per axis, and `mem.*`'s `by_name` must list every member once in
strictly increasing name order; a reader refuses anything else. (A fuzzer found that an
unchecked `by_name` with a repeated id overran the group lists built from it.)

### 2.5 Metadata and optional sections

**Metadata** (`meta`, required, loaded at open): string key → string value.

```
  u32 n
  n x { u32 key_len, key, u32 value_len, value }     sorted by key, keys unique
```

Keys are 1–255 bytes, and neither keys nor values contain NUL; anything else is allowed,
tabs and newlines included. Standard keys — a library-level vocabulary so tools can describe
any liftmap without knowing its profile:

| key | meaning |
|---|---|
| `axis.a`, `axis.b` | what each axis is: `child`/`parent`, `query`/`target`, `sequence`/`universal column` |
| `max_gap` | runs were merged across gaps of up to this many bp; coordinates are approximate. Absent or `0` = exact |
| `tree` | Newick tree over the groups' names |
| `source`, `program`, `description` | provenance |

Applications prefix their own keys (`tui.format`, `tui.anchors`).


**Application sections.** Any section whose id begins `x.` belongs to the application: the
library stores it verbatim, checksums it like every other section, and never interprets it.
It is loaded only when asked for, so a large one costs nothing at open. Ids are unique; a
reader refuses a duplicate. This is for bulk profile data that does not belong in metadata:
`.tui` keeps its column → MAF file-offset anchors, used by `view -U`, as `x.tui.anchor`.

Sections the library itself may define later, not yet implemented:

| id | kind | for |
|---|---|---|
| `paralogy` | `(cluster_id, member, start, len)` records | hal2 paralogy as equivalence classes |
| `stats` | counts, coverage, run-length and gap histograms | O(1) `halStats` |

Member names can also be looked up by prefix (the members whose name begins with a string
form one contiguous range in name order), but a group is the reliable way to name a
genome's sequences (2.4a).

---

## 2.6 Queries (normative semantics)

A query names one member on one axis and a half-open window `[lo, hi)` on it, and returns
every run overlapping the window **clipped to it**, with the other axis mapped through the
run's strand. For a window of offsets `[o1, o2)` into a run of length `len`, the other axis
gets `[o1, o2)` on the forward strand and `[len-o2, len-o1)` on the reverse strand.

- **axis a** (`query_a`): binary search the member's `dir.a` range for the last entry with
  `a_min < hi`; walk backwards, decoding each chunk whose extent overlaps the window, and stop
  at the first entry whose `dir.a.max <= lo`. Results are sorted by `(a, b, len, strand)`
  and never overlap, because axis a is a partial function — unless the file sets
  `A_OVERLAP`, when they may. Under order a a point query decodes exactly one
  chunk; under order b it decodes the few whose axis-a extents overlap (4.3 on average on the
  fish edge, against 1.18 for an axis-b point query — the trade order b makes).
- **axis b** (`query_b`): binary search the member's `dir.b` range for the last entry with
  `b_min < hi`; walk backwards, decoding each chunk whose extent overlaps the window, and
  stop at the first entry whose `prefix_max_b_end <= lo`. Results are sorted by
  `(b, a_member, a)`. They may overlap on axis b: that is paralogy, and a query reports it
  rather than resolving it. Choosing a representative — the smallest `a` among the copies
  covering a position — is a caller's policy, not the format's.

Measured on the fish edge (9,010,509 runs, 5,575 parent scaffolds), mean chunks decoded per
query: axis a 1.00 at a 2 bp window, 3.0 at 200 kb; axis b 1.35 at 2 bp, 2.3 at 200 kb.

**Cursor.** Both queries are one-shot wrappers over a cursor, which is the general access
path: one member of a key axis, optionally restricted to a set of members of the other axis
(tui restricts the column axis to one genome's sequences). It returns whole runs in key
order, ties broken by (other member, other position, len, strand), *across* chunks that
overlap: chunks are opened in `key_min` order only as the merge reaches them, so a scan holds
just the chunks overlapping its position. Decoded chunks are kept under a byte budget,
least recently used first out; budget 0 streams. Point lookups use the same cache and search
forward from the previous answer, so a per-position sweep costs O(1) per position. A cursor
is single-threaded; the file is not, so concurrent readers each open their own.

Measured through taffy on the rodent universal index against its previous ONEcode reader:
per-column reverse lift 0.94–0.98× the time, range lift at parity, output identical except
the order of paralogous matches at a column, which was unspecified and layout-dependent.

## 3. The two profiles

### 3.1 `hal2.edge` — a child genome onto its parent

| | |
|---|---|
| axes | a = child sequences, b = parent sequences (`axis.a` = `child`, `axis.b` = `parent`) |
| order | a |
| chunks | count 8192, bspan 10^6 |
| groups | none needed: one genome per axis |

Axis a is strictly sorted and disjoint (a child base has 0 or 1 parent), so `a` is `uvarint`.
Axis b may repeat and reorder — that is paralogy and rearrangement — so `b` is `zigzag`.
**Paralogy is not stored.** The copies of the child at parent position *p* are exactly what a
`dir.b` query at *p* returns; `--noDupes` picks the smallest `a` among them at query time.
Measured against HAL's stored canonical on a real fish edge: agreement on 100.000% of 876,531
duplicated positions, zero differences.

### 3.2 `taffy.tui` — a genome onto the universal column axis

| | |
|---|---|
| axes | a = `genome.sequence` of every genome, b = one member `column` of length T |
| order | b |
| chunks | count 65536, bspan 10^6 (`TUI_CHUNK_RUNS`, `TUI_CHUNK_G_MAX`) |
| groups | axis a, one per resolved genome, in genome order |
| metadata | `tui.format`, `tree` (the `# hal` Newick), `max_gap`, `tui.anchors`, `axis.*`, provenance |
| sections | `x.tui.anchor`: the column → MAF file-offset index |

A chained index (`max_gap > 0`) sets `A_OVERLAP`; a base index does not. T is the column
member's length rather than a separate field.

That is the whole difference between the two profiles: **run order, chunk constants,
metadata, groups, and which application sections are present.**

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

**Validated against taffy's `.tui` reader** (`tests/lmap_transcode.c` in the taffy
checkout, up to commit 09ced8e; retired once `.tui` itself moved onto this library, below). Every run of `tests/tui/evolverMammals.uni.maf.gz.tui` (16 sequences, 123,471
runs, T = 708,019 columns) was read with tui's own reader and written through this library:

- run identity: every sequence's runs come back from `query_a` exactly
- forward lift: `tui_query` vs `query_a`, base by base — 5,041,395 bases over whole sequences
  plus 30,433,565 over 2,000 random windows, 0 mismatches
- reverse lift: `tui_genome_lift_column` vs `query_b` over every column and all 9 genomes —
  5,041,395 matches, 0 mismatches

Written in `order a` that transcode came out 43% larger than the `.tui` (3.12 vs 2.18
B/run), exactly as the order table above predicts. **Written in `order b` — the `taffy.tui`
profile — it is 267,804 bytes against the `.tui`'s 268,906: parity, 0.4% smaller**, with the
same three checks at zero mismatches. With per-stream compression it is **263,492 bytes,
2.0% smaller than the `.tui`**.

The same check on the rodent universal index (`vgp-577way-v1-MuridaeAnc3`, regenerated in
tui format 0.3: 11,698 sequences, 84,970,744 runs, T = 2,731,506,489 columns):

| | `.tui` | `.lmap`, order b |
|---|---|---|
| size | 235,525,448 bytes | **226,116,332 bytes (−4.0%)** |
| forward lift, 11,243 whole sequences (1.22 × 10^9 bases) | 0.361 s | 0.050 s |
| forward lift, 2,000 random windows | 0.113 s | 0.056 s |
| reverse lift, 8.02M columns × 5 genomes | 5.57 s | 0.455 s |
| every run of every sequence | 9.13 s | 14.4 s |

Run identity, forward lift and reverse lift all agree exactly. Transcoding took 42 s and
0.58 GB peak memory.

And on the fish-subtree universal index (`vgp-577way-v1.RayFinnedFishesAnc64`, regenerated
in tui format 0.3: 57 genomes, 120,805 sequences, 846,503,338 runs, T = 2,655,269,561):

| | `.tui` | `.lmap`, order b |
|---|---|---|
| chunks | 354,369 | **354,369** — identical, so the layout matches rather than approximates |
| size | 1,956,422,670 bytes | **1,911,774,892 bytes (−2.3%)** |
| forward lift, 118,931 whole sequences (2.25 × 10^9 bases) | 7.13 s | 2.76 s |
| forward lift, 2,000 random windows | 0.383 s | 0.355 s |
| reverse lift, 8.04M columns × 57 genomes | 40.3 s | 3.58 s |
| every run of every sequence | 100.5 s | 150.4 s |

All three checks agree exactly: 846,503,338 runs, 2,245,605,736 bases of forward lift and
50,889,391 reverse matches. Transcoding took 766 s at 4.8 GB peak, most of it the writer
holding the ~1.9 GB payload (see the known limit above).

**taffy's `.tui` now runs on this library** (tui format 0.4, taffy.intervalmap `876ab73`),
checked against the ONEcode build of taffy on the same inputs:

| | rodent | fish subtree |
|---|---|---|
| `.tui` size, ONEcode → liftmap | 235.5 → 227.9 MB | 1,956 → 1,917 MB |
| `taffy index -u` wall | 2:46 → 2:46 | 37:50 → 39:00 |
| open (format 3: directories paged) | 2.9 ms | 8.4 ms (354,369 chunks, 120,805 sequences) |
| `view -U`, 30 regions × 3 modes | 90/90 identical, 5.5 → 5.4 s | 90/90 identical, 35.3 → 35.7 s |
| `lift`, 300 BED intervals × 4 modes | 5 genomes, 19/20 identical* | 4 genomes, 16/16 identical |
| `lift -F --chainOverlapFrac 0` | 1.3–1.6 → 0.7–1.3 s | 1.4–2.9 → 0.9–1.5 s |

\* The one difference: where paralogous matches meet at a column, taffy's default lift
splits rows by the order matches arrive, which the ONEcode reader left to the chunk layout
and the cursor fixes; coverage is identical.

Before format 3 every chunk directory was read at open: 72 ms on fish, repeated per
invocation, which made `view` 14% slower there. Paging the directories removed it; what
remains at open is mostly the member directory (120,805 sequence names on fish), which is
also where the next gain at 577-way scale would come from.

Still open, to be measured rather than argued:

1. *(settled: order a for HAL edges, order b for universal columns — table in 1.3)*
2. *(settled: no — the strand-aware delta hurts under order b)*
3. **zstd level, and whether a per-section trained dictionary earns its keep.** A dictionary
   measured 2.55 vs 2.57 B/run at 65536 runs but 2.82 vs 3.13 at 256 — it matters only if
   chunks end up small.
4. **`b` as `uvarint` under `order b`.** If b-ordered chunks make the enter-delta monotone it
   can drop the zigzag bit.

---

## 5. Remote access

A liftmap file is meant to be read where it sits: on a web server, in a bucket, behind a
browser (UCSC blockViz and track hubs). ONEcode could not do this at all, since it is
stdio-only, and that was one of the reasons to leave it. Nothing in the format has to
change for it. Every stored offset is file-relative (constraint 1), all reads go through
`lmap_pread` (constraint 2), and open touches only the trailer, footer and member
directories (constraint 3). What remains is transport, and making each round trip
count.

### 5.1 Design

**The seam and the cache live in liftmap; transports mostly live in the client.**

- `lmap_io_open_backend(be, ctx, size)` turns any `read(ctx, buf, off, len)` into a
  handle. That is how a client plugs in the remote stack it already has: taffy uses
  htslib's hFILE for the schemes only it speaks. liftmap depends on none of them.
- `lmap_io_cache(inner, block, capacity, readahead, own)` goes in front of any handle.
  - Blocks are 64 KiB, with 64 MiB of them held least-recently-used.
  - Adjacent misses are coalesced into one request.
  - Sequential misses read ahead, doubling up to 4 MiB.
  - A read bigger than a quarter of the cache bypasses it.
  - It is thread-safe, and `lmap_io_get_cache_stats` reports requests, bytes, hits and
    misses.
  - Directory pages (4 KiB, §2.3a) and chunks (tens of KiB) are both smaller than a block,
    so one request usually serves several of them.
- `lmap_io_open_url` is an optional libcurl backend (`make HTTP=1`). Without it, the core
  still needs only zlib and pthreads. It makes bounded range requests
  (`Range: bytes=off-end`) on one reused connection, and it gets the size from `HEAD`.
  `lmap_open` uses it, behind the cache, for any path containing `://`, so `bin/liftmap`
  reads URLs directly.

**Use bounded ranges.** A transport that asks for an open-ended range (`bytes=off-`) and
hangs up once it has enough makes the server send far more than is read. The server keeps
writing into socket buffers until the connection drops, and object stores bill for those
bytes.

Measured on the rodent universal index over a local range server, with `taffy lift`
on 20 BED intervals:

| transport | requests | open-ended | bytes sent by the server | output |
|---|---|---|---|---|
| htslib hFILE (open-ended ranges) | 71 | 69 | ~150 MB | identical |
| liftmap libcurl (bounded ranges) | 85 | 0 | **8.8 MB** | identical |

Wall time was the same on loopback. On a real network, the extra ~140 MB is the whole
difference. taffy therefore uses liftmap's backend for `http://` and `https://`, and
hFILE only for the schemes it alone speaks (`s3://`, `gs://`, ...).

hFILE also has a correctness bug, found here. In htslib 1.17, a forward seek of less than
1 MB is served by reading and discarding data (`hfile_libcurl.c`, `libcurl_read`). When a
network read ends exactly where the skip ends, `to_skip` is never reset, so the next read
skips the same amount again and returns bytes from the wrong offset, with no error.
liftmap's per-chunk CRCs caught it (§2.2): the error surfaced as "corrupt run data"
rather than as wrong lift output. taffy's hFILE backend works around the bug by doing
short forward skips itself. This is a second reason to validate payload checksums on
every remote read.

**Measured costs** (loopback range server, default cache):

| | requests | bytes |
|---|---|---|
| open, rodent (11,698 sequences) | 5 | 1.2 MB |
| open, fish subtree (120,805 sequences) | 6 | 11.7 MB |
| `taffy lift`, 20 intervals, rodent | 70–85 | 7.9–8.8 MB |
| `taffy view -U`, 2 kb region, rodent | 15 | ~11 MB, mostly the `.maf.gz` blocks, not the `.tui` |

In `tests/remote.c` on rodent, 500 random queries plus cursor sweeps took 1,039
uncached reads and 613 through the cache. With 20 ms of latency per read, that is 7.4 s
against 4.0 s. A repeated workload makes 0 new reads. The cache trades bytes for round trips. It fetches whole
64 KiB blocks, so that workload moved 39.7 MB against 13.5 MB uncached. On any network
where latency dominates, that is the right trade, and `block_bytes` is the knob.

Correctness is checked in two places. `tests/test_tools.py` serves files from an
in-process HTTP server with Range support and requires `dump`, `info` and `lift` over
`http://` to match the local results. `tests/remote.c` runs the same workload directly,
through a custom backend and through the cache, and requires identical hits.

### 5.2 Roadmap

| phase | what | status |
|---|---|---|
| R1 | `lmap_io_open_backend` seam; `lmap_io_cache` block cache with coalescing, read-ahead and stats | **done** |
| R2 | optional libcurl backend, `lmap_open` on URLs, `bin/liftmap` over http(s) | **done** |
| R3 | taffy: remote `.tui` for `view -U` and `lift`; libcurl for http(s), hFILE for other schemes | **done**, identical output to local on rodent |
| R4 | HAL format 3: `.hgn` genome files, and `.2bit` DNA through twobit64, both over liftmap's libcurl backend and `lmap_io_cache`, not udc. That keeps UCSC kent code out of the dependency tree, and gives bounded range requests. The manifest resolves relative locators against its own URL (as `LodManager::resolvePath` does); blockViz reads a remote `.hal`. twobit64 already reads URLs this way. | planned |
| R5 | performance at 577-way and beyond, below | planned |

**R5 items**, in rough order of expected payoff:

1. **Page the member directory.** It is most of the open cost: 11.7 MB for 120,805 fish
   sequences, and it grows with sequence count, not with what a query touches. Page
   `mem.a`/`mem.b` the way the chunk directories already are (§2.3a), with a small top
   level holding the first name of each page. Open then costs a fixed few KB, and a
   lookup by name costs one page. This needs a format minor bump, and a reader keeps
   accepting the unpaged layout.
2. **Parallel fetches.** Today the cache issues one request at a time. The next step is
   the curl multi interface, with HTTP/2 where the server offers it, so that a lift over
   many intervals, or a `view` needing several chunks, overlaps its round trips. The
   cursor already knows which chunks it will need next (the heap merge of overlapping
   chunks), so it can prefetch them.
3. **Open in one round trip.** The trailer is 32 B at the end, and the footer and member
   directories sit just before it. So open could fetch the last N KB speculatively with a
   suffix range (`bytes=-N`), sized from a hint in the header, and return from a single
   request when everything fits.
4. **Persistent disk cache**, content-addressed by URL, `ETag` or `Last-Modified`, and
   block index, so a restarted process or a second tool reuses what was fetched. Browser
   use of HAL format 3 needs this: blockViz under UCSC has relied on udc's disk cache,
   and this replaces it. So it moves ahead of items 2 and 3 once R4 starts.
5. **Remote container slices.** A manifest may point into a tarball or zip64 archive at
   an offset (constraint 1). `lmap_io_slice` over a remote handle already covers the
   read side. What remains is resolving a member name to an offset through the
   container's index, and fetching that index once.
6. **Validation against real object stores** (S3, GCS, a CDN) over real latency.
   Loopback numbers bound request counts, not wall time.

The format commits to none of these. R5.1 is the only item that touches the file layout.
