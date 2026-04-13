# 07 — Implementation Phases

This document proposes a layered rollout of the architecture described
in [06-final-architecture](./06-final-architecture.md). Each phase
delivers working functionality that builds on the previous; the layers
can be shipped independently and each provides value on its own.

## Phase 0 — Preparatory refactors (optional)

Small, independent cleanups that simplify later phases:

- Make `FeatureReference` an extensible enum (already a sum type;
  no change needed, but pattern-matching sites should be reviewed).
- Extract the per-tile `write_inputs` logic (`src/main.rs:171-238`)
  from the export loop so it can be replaced with per-cell
  concatenation later.
- Introduce a typed `Format` enum that replaces the current string
  handling in geof-subprocess parameters.

None of these are required; they make subsequent phases smaller.

## Phase 1 — Geof Rust library + CLI

**Goal:** replace the geof subprocess with a native Rust crate. Tyler
continues to use tile-grained conversion (no streaming yet).

**Deliverables:**

- `geof` crate with `convert_tile()` (layer 1 API).
- `geof` standalone CLI (single input → single output).
- Initial format support: 3D Tiles (GLB). CityJSON / OBJ / GPKG can be
  added incrementally.
- Tyler's `src/main.rs` tile-export loop uses `geof::convert_tile()`
  directly instead of spawning subprocesses.
- Remove `subprocess` crate dependency.

**Value delivered:** eliminates subprocess spawn overhead per tile.
Simpler deployment (no external `geof` binary). Typed error handling.

**Does not require:** stdin support, morton-sort, progressive
assembly. Tyler's pipeline is otherwise unchanged.

**Acceptance tests:**
- Run tyler on an existing cjindex dataset; compare output tileset
  to the subprocess-based baseline for a golden-file comparison.
- Performance comparison: wall-clock time and max RSS.
- CLI: `geof input.city.jsonl output.glb` produces a valid GLB.

## Phase 2 — Stdin ingestion with per-cell spill files (single-pass)

**Goal:** add `--stdin` mode with pre-declared extent. Per-cell spill
files replace in-memory geometry. Assembly still at EOF (no
progressive).

**Deliverables:**

- New `InputSource::Stdin` variant.
- CLI flags: `--stdin`, `--extent x1,y1,z1,x2,y2,z2`.
- Header parsing: the first CityJSONSeq line populates CRS / transform;
  if `metadata.geographicalExtent` is present, it populates the extent.
- Per-cell spill file writer. Format: raw CityJSON bytes (for Phase 2;
  format-specific intermediates come in Phase 4).
- `FeatureReference::SpillRef { cell_id, byte_offset, byte_len }`.
- Modified `write_inputs` path that concatenates cell spill files into
  per-tile NDJSON.

**Value delivered:** tyler works with `roofer | tyler --stdin
--extent ...`. Memory footprint drops to grid + feature index +
bounded buffers, independent of input size.

**Does not require:** morton-sort, progressive assembly.

**Acceptance tests:**
- `cat existing.city.jsonl | tyler --stdin --extent ... -o out/`
  produces the same tileset as running tyler on the file-based dataset.
- Memory ceiling verified with a large input.

## Phase 3 — Two-pass fallback for missing extent

**Goal:** allow tyler to work without a supplied extent, by running an
automatic two-pass mode.

**Deliverables:**

- When `--stdin` is used without `--extent` and no header extent is
  present: spill raw stdin lines to a temp file, compute extent
  incrementally, then re-process the temp file via the Phase-2 path.
- CLI flag: `--stdin-temp-dir` to override the default temp location.

**Value delivered:** ergonomic safety net. Users can `roofer | tyler
--stdin` without thinking about extents.

**Acceptance tests:**
- Same output as Phase 2 given equivalent input.
- Automatic cleanup of the temp file on success.

## Phase 4 — Format-specific per-feature encoders (progressive layer 1)

**Goal:** move per-feature conversion work (triangulation) into the
streaming loop. Spill files contain format-specific intermediates.

**Deliverables:**

- `FeatureEncoder` and `TileAssembler` traits in the geof crate (layer
  2 API).
- Implementations for GLB, CityJSON, OBJ (and GPKG as a stretch goal).
- Tyler's ingestion loop calls `encode_feature()` inline and writes the
  format-specific intermediate to the cell spill file.
- Tyler's assembly calls `assemble_tile()` with iterators over
  pre-encoded intermediates.
- The `convert_tile()` facade is re-implemented on top of the two-trait
  API.

**Value delivered:** triangulation happens progressively during
streaming. Post-EOF assembly becomes a fast buffer-packing operation.
For non-triangulated formats, the per-feature step is a cheap
passthrough.

**Acceptance tests:**
- Each encoder produces byte-identical output to the previous
  monolithic conversion (golden files).
- Per-feature / per-tile split does not change the output.
- Throughput comparison: ingestion time, total wall-clock time.

## Phase 5 — Live quadtree + progressive assembly

**Goal:** trigger assembly for sealed leaves as soon as their cells
complete. Output files appear during ingestion, not only at EOF.

**Deliverables:**

- `LiveQuadtree` mutable structure tracking per-node vertex counts and
  sealed flags.
- Capacity-sealing logic: when `sum > capacity`, mark group sealed.
- Rayon assembly worker pool and a channel from the ingestion thread.
- At EOF: finalize all non-sealed quadrants, assemble remaining tiles,
  write `tileset.json`.

**Dependency on Phase 6:** without morton ordering, progressive
assembly works only for sealed leaves whose feature set happens to be
complete. In practice: for dense cells that exceed capacity early AND
have all their features by the time they're sealed. Partial benefit.

**Value delivered:** some assembly work moves earlier. Most benefit
materializes with Phase 6.

**Acceptance tests:**
- Tile output files match the batch pipeline byte-for-byte.
- Progressive output observable (tiles appear on disk during ingestion).

## Phase 6 — Morton-sort tool and low-water mark

**Goal:** add the standalone `morton-sort` tool and tyler's low-water
mark logic. Full progressive assembly is now enabled.

**Deliverables:**

- `morton-sort` crate (separate binary in the workspace).
- Reorder buffer implementation with configurable size.
- Header + control-line passthrough.
- Tyler's streaming loop advances an LWM per feature; the
  `LiveQuadtree` supports `on_lwm_advanced()` to detect completed
  quadrants.

**Value delivered:** tiles are assembled progressively as the morton
frontier advances. For spatially coherent input, the vast majority of
assembly work overlaps with ingestion. EOF work is limited to sparse
regions and `tileset.json`.

**Acceptance tests:**
- `morton-sort` preserves semantic content (every input feature is
  emitted exactly once).
- With a buffer size of ~10K on a spatially coherent input, the output
  is in strict morton order.
- Tyler's behavior is unchanged whether or not `morton-sort` is in the
  pipe (progressive benefit differs; correctness does not).

## Phase 7 — Optional control signals

**Goal:** support explicit `+RegionComplete` signals for producers
that can emit them.

**Deliverables:**

- CityJSONSeq parser recognizes `{"type": "+RegionComplete", "bbox":
  [...]}` lines.
- `morton-sort` treats control lines as barriers (drain buffer before
  passing through).
- Tyler's streaming loop handles `+RegionComplete`: mark all cells
  fully contained in the bbox as drained, trigger assembly for newly
  complete quadrants.
- Documentation on the control-line protocol (forward-compatible for
  future signals).

**Value delivered:** producers that can cheaply declare regions
complete (e.g., a tile-at-a-time processor) get progressive behavior
even without morton ordering.

**Acceptance tests:**
- Signals correctly trigger assembly.
- Unknown `+`-prefixed types are ignored (forward compatibility).
- Output is identical whether signals are present or not (only timing
  differs).

## Phase 8 — Non-triangulated formats

**Goal:** complete support for CityJSON, OBJ, GPKG output.

**Deliverables:**

- Full `FeatureEncoder` / `TileAssembler` implementations for each
  format.
- CLI support: `--format cityjson`, `--format obj`, `--format gpkg`.
- Per-format output layout decisions (e.g., a single monolithic
  CityJSON vs per-tile CityJSON; GPKG schema choices).

**Value delivered:** tyler becomes a general-purpose streaming
tile-generation tool, not specific to 3D Tiles.

**Acceptance tests:**
- Each format's output validates against its spec.
- Spot-check visual inspection (a few tiles loaded in a viewer).

## Proposed ordering summary

| Phase | Theme | Blocked by |
|---|---|---|
| 0 | Preparatory refactors (optional) | — |
| 1 | Geof Rust library + CLI | (optional) 0 |
| 2 | Stdin ingestion + per-cell spill | 1 |
| 3 | Two-pass fallback | 2 |
| 4 | Format-specific per-feature encoders | 1, 2 |
| 5 | Live quadtree + progressive assembly | 4 |
| 6 | Morton-sort + LWM | 5 |
| 7 | Control signals (optional) | 5 |
| 8 | Non-triangulated formats | 4 |

Phases 1–3 are the critical path to basic streaming support.
Phases 4–6 deliver the full progressive architecture.
Phase 7 is optional enhancement.
Phase 8 is independent format expansion.

## Risk areas and mitigations

**Risk: geof library triangulation correctness.**
Mitigation: extensive golden-file testing in Phase 1 against the
existing subprocess output. Keep the CLI available as a fallback.

**Risk: per-cell spill file count grows unbounded in pathological
datasets.**
Mitigation: the grid has at most `length²` cells. For the Netherlands
this is ~1.5M cells but only a subset will have features; empty cells
have no spill files. File descriptor exhaustion can be avoided with a
bounded LRU of open `BufWriter`s.

**Risk: live quadtree complexity bugs.**
Mitigation: extensive unit tests on sealing transitions. Run in
parallel with the batch pipeline for a period and compare outputs.

**Risk: `morton-sort` buffer size insufficient for real roofer
output.**
Mitigation: `--on-late` handling, metrics to detect buffer pressure,
configurable size, fallback to `+RegionComplete` signals from roofer
if sorting is impractical.

**Risk: cross-feature optimizations inside the GLB packer would make
per-feature encoding incorrect.**
Mitigation: verify during Phase 4 that the triangulation /
attribute-encoding of feature A is independent of feature B. If not,
demote GLB to per-tile encoding (still works; loses Phase-4 progressive
benefit for GLB only).

## Non-phases (work we are deliberately not doing)

- Replacing the `SquareGrid` with a sparse data structure. Orthogonal
  optimization; not required for streaming.
- Replacing morton order with Hilbert order. Morton matches the
  existing quadtree; switching costs exceed benefits.
- Unbounded / infinite stream support. Tyler is designed for bounded
  inputs; EOF is required for `tileset.json`.
- Distributed / multi-machine execution. Out of scope.
