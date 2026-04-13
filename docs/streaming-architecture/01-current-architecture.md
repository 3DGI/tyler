# 01 — Current Architecture Analysis

This document describes tyler's current pipeline and identifies the
specific constraints that make streaming input non-trivial.

## Current pipeline (simplified)

```
Extent computation  →  Grid allocation  →  Feature indexing  →  Quadtree  →  Tileset  →  Tile export
    (all features)     (needs extent)      (needs grid)        (needs all   (needs      (needs
                                                                 cells)      quadtree)   features)
```

The pipeline is strictly sequential. Each stage must complete before the
next begins.

## Input sources

Today tyler accepts two input source types (`src/parser.rs:77-86`):

```rust
pub enum InputSource {
    LegacyFeatureFiles { features_root: PathBuf },
    CjIndexDataset { dataset_root, index_path, layout },
}
```

The `CjIndexDataset` variant covers three cjindex layouts: `Ndjson`,
`CityJson`, and `FeatureFiles` (`src/parser.rs:88-93`). All current
input paths require random-access file storage; there is no stdin support.

## Phase-by-phase analysis

### Phase 1: Extent computation (`src/parser.rs:160-258`)

Walks every feature (or every bbox page in cjindex) to compute the
spatial extent. Takes ~45% of total runtime on a Netherlands-scale dataset
per `docs/performance/2026-03-30-three-ref-time-comparison.md`.

**Why this is a streaming bottleneck:** the extent is required before the
grid can be allocated. With an unknown extent, the grid size is unknown.
With stdin streaming, all features must arrive before the extent is known.

### Phase 2: Grid allocation (`src/spatial_structs.rs:495-573`)

Creates a **dense** 2D grid: `Vec<Vec<Cell>>` sized from the extent:

```rust
let mut row: Vec<Vec<Cell>> = Vec::with_capacity(d_cells);
row.resize_with(d_cells, || {
    let mut column: Vec<Cell> = Vec::with_capacity(d_cells);
    column.resize(d_cells, Cell { feature_ids: Vec::new(), nr_vertices: 0 });
    column
});
```

Grid dimensions must be a power of 2 (`d_cells = 2^n`) so the quadtree
has `4^n` cells. For the Netherlands at 250m cellsize: ~1250×1250 cells,
estimated 2.6GB with 10M features (per `docs/design_document.md`).

### Phase 3: Feature indexing (`src/parser.rs:596-646`)

For each feature:
1. Parse the CityJSONFeature
2. Extract bbox, centroid, per-cell vertex counts
3. Assign feature to one or more cells
4. Store lightweight `Feature { centroid: [f64;2], bbox: [f64;6], reference }`
   in `World.features` (the `FeatureSet` = `Vec<Feature>`)

**Key observation:** tyler already does not hold full geometry data in
memory during indexing. Only a ~72-byte `Feature` struct per feature is
retained. The full CityJSON is parsed for stats and then discarded.

The feature's `FeatureReference` (`src/parser.rs:824-828`) is either a
`LegacyPath(PathBuf)` (relative path to the feature file) or a
`CjIndexId(String)` — both of which enable re-reading the geometry from
disk during tile export.

This is the correct design for memory efficiency. The streaming
architecture should preserve this property.

Takes ~46.6% of total runtime (the second heavy pass over the dataset).

### Phase 4: Quadtree construction (`src/spatial_structs.rs:40-148`)

Builds a morton-ordered quadtree bottom-up from the grid. Groups of 4
sibling cells are tested for merging:

```rust
if sum_items <= limit {
    // MERGE: 4 children collapse into one leaf tile
    QuadTree { children: vec![], cells: all_cells, nr_items: sum_items }
} else {
    // SPLIT: 4 children remain separate tiles
    QuadTree { children: tiles, cells: vec![], nr_items: sum_items }
}
```

This is the critical merge decision. It requires **complete** vertex
counts per cell — the decision cannot be made partially.

**The monotonicity property** (important for streaming): vertex counts
only increase as features arrive. Once `sum_items > capacity`, the merge
decision is permanent (the group will never merge). This property enables
progressive sealing, detailed in
[03-progressive-tile-generation](./03-progressive-tile-generation.md).

### Phase 5: Tileset generation (`src/main.rs:543-629`)

Converts the quadtree into a 3D Tiles tileset structure. Fast (<5% of
runtime). Produces `tileset.json` and optionally subtree files for
implicit tiling.

### Phase 6: Tile export (`src/main.rs:648-850+`)

For each leaf tile, in parallel (Rayon `into_par_iter`):

1. Collect feature IDs from the tile's cells
2. Read each feature's geometry from disk (file or cjindex)
3. Write a per-tile NDJSON file (`write_inputs`, `src/main.rs:171-238`)
4. Spawn a `geof` subprocess with the NDJSON path and tile parameters
5. Receive the produced GLB file

This is the most expensive phase in a full run. The geof subprocess
performs CityJSON parsing, triangulation, and GLB packing.

## Current memory profile

| Data | Size | Notes |
|---|---|---|
| Grid (`Vec<Vec<Cell>>`) | ~O(cells) | Dense 2D array; ~2.6GB for the Netherlands at 250m |
| Feature set (`Vec<Feature>`) | N × ~72 bytes | Lightweight index; ~720MB for 10M features |
| Geometry data | **0 bytes in RAM** | Held on disk, re-read during tile export |

The on-disk requirement for geometry is the constraint the streaming
design must address.

## Output formats

Today tyler supports only 3D Tiles output. CityJSON output is stubbed
but unimplemented (`src/main.rs:371` panics with "cityjson output is not
supported"). The design must accommodate future OBJ, CityJSON, and GPKG
outputs — some of which do not involve triangulation.

## What tyler gets right (to preserve)

1. Lightweight per-feature memory footprint during indexing.
2. Parallel tile export via Rayon.
3. The grid + quadtree spatial indexing scheme.
4. Morton-ordered cell iteration.
5. Clear separation of concerns between tyler (tiling) and the
   conversion engine (currently geof subprocess).

## What needs to change for streaming

1. **Extent discovery without a pre-pass.** The extent must be provided
   upfront (CLI argument or CityJSONSeq header) or inferred via a
   fallback two-pass mode.
2. **Ephemeral input bytes must be persisted somewhere.** Since stdin
   data cannot be re-read, features must be written to disk during
   ingestion. The current `FeatureReference::LegacyPath` and
   `CjIndexId` mechanisms assume pre-existing persistent storage.
3. **Conversion must decouple from subprocess.** To perform per-feature
   work progressively during streaming, the conversion engine must be
   callable at feature granularity, not only at tile granularity. A
   subprocess-per-feature is prohibitively expensive; an in-process
   library is required.
4. **The pipeline must support progressive assembly.** The current
   batch assumption (build everything, then export tiles) must be
   generalized to handle tiles that can be assembled as soon as their
   inputs are complete.
