# 02 — Alternatives Considered for Streaming Input

This document enumerates the streaming-input alternatives that were
evaluated and explains why each was accepted or rejected.

## Evaluation criteria

Alternatives were evaluated against four criteria, in order of priority:

1. **Memory efficiency** — footprint bounded by grid + per-feature index,
   independent of geometry size.
2. **Cleanest, most elegant architecture** — minimum mechanism; streaming
   concern isolated from tiling concern.
3. **Performance** — single-pass ingestion, zero re-serialization,
   parallelism preserved.
4. **Compatibility with existing pipeline** — the `SquareGrid` +
   `QuadTree` + `Tileset` machinery should be preserved unchanged.

## Alternative 1 — Pre-declared extent + single-pass stdin (ACCEPTED as base)

User (or CityJSONSeq header) provides the spatial extent. Tyler allocates
the grid immediately, then reads CityJSONSeq from stdin line by line.
Each feature is parsed for stats, indexed into grid cells, and its raw
bytes are appended to a per-cell spill file on disk. When stdin reaches
EOF, tyler builds the quadtree and generates tiles.

**Pros:**
- True single-pass streaming. Bounded memory for geometry (zero bytes in
  RAM; spilled to disk).
- Minimal change: the grid/quadtree/tileset/export pipeline stays as-is.
- Works naturally with `roofer | tyler --extent ...` pipelines.
- Per-cell spill files double as the input to the tile-assembly phase —
  no re-serialization, no re-parsing. Tile export becomes file
  concatenation instead of a parse/serialize roundtrip.

**Cons:**
- User must supply an extent (or it must be embedded in the CityJSONSeq
  header).
- Features outside the declared extent must be handled (warn or error).

**Verdict:** accepted as the base streaming architecture. The per-cell
spill file design elegantly unifies two roles: persistent storage for
ephemeral stdin data, and pre-grouped input for tile assembly.

## Alternative 2 — Two-pass stdin with spill file (FALLBACK)

First pass: read stdin, write every feature to a single spill file,
compute extent incrementally. Second pass: run the existing batch
pipeline against the spill file.

**Pros:**
- No user-supplied extent needed.
- Core tiling pipeline is completely unchanged.

**Cons:**
- Not truly streaming — tile generation waits for all input to arrive.
- Extra disk I/O for the monolithic spill file.
- Random access during export requires seeking to arbitrary byte offsets.
- The single-file design cannot be directly fed to the converter without
  a second rewrite step.

**Verdict:** accepted as an automatic fallback when no extent is
provided. Ergonomic safety net; not the designed-for path.

## Alternative 3 — Stdin → cjindex SQLite database (REJECTED)

Stream features from stdin into a cjindex SQLite database. Once stdin
closes, tyler reads via the existing `CjIndexDataset` code path.

**Pros:**
- Leverages existing cjindex infrastructure.
- Per-feature bboxes stored in SQLite enable fast extent queries.

**Cons:**
- SQLite transactions, WAL journaling, and B-tree insertions add
  overhead on the ingestion hot path that buys indexed random access
  (which isn't needed during linear ingestion).
- Adds SQLite as a required component even for simple piped workflows.
- Same EOF-blocked generation as Alternative 2.

**Verdict:** rejected. The per-cell spill file design achieves the same
persistence goal with strictly less mechanism. SQLite is the right
abstraction for a *database*; append-only files are the right abstraction
for a *spill buffer*.

## Alternative 4 — Adaptive / growing grid (REJECTED)

Replace the dense `Vec<Vec<Cell>>` grid with a `HashMap<CellId, Cell>`,
or start with a small grid and resize as features arrive outside its
bounds.

**Pros:**
- No extent required upfront.
- Could save memory when the dataset is spatially sparse.

**Cons:**
- Replacing the dense grid with a hash map slows down `locate_point`
  from O(1) array indexing to hashed lookup.
- Resizing while maintaining the power-of-4 cell count invariant is
  complex.
- The memory cost of empty cells in the dense grid is modest (~48 bytes
  per `Cell`); the justification for a sparse structure is weak.

**Verdict:** rejected as the base streaming mechanism. Could be
revisited as an orthogonal optimization for pathologically sparse
datasets, but should not be coupled with the streaming design.

## Alternative 5 — Hash-based spatial partitioning (REJECTED)

Replace the grid + quadtree with a spatial hash (geohash, S2, Hilbert
curve). Each feature assigned to a hash bucket on arrival; buckets form
the tile hierarchy.

**Pros:**
- No extent needed; hashes work on global coordinates.
- O(1) per-feature assignment.

**Cons:**
- Incompatible with tyler's projected-coordinate grid and morton-order
  quadtree. The entire tileset generation and 3D Tiles export pipeline
  assumes a square grid with `4^n` cells.
- Requires rewriting the quadtree, tileset builder, and implicit tiling
  support.
- No clear benefit for projected-coordinate datasets.

**Verdict:** rejected. A ground-up redesign with no measurable advantage
over the grid approach.

## Alternative 6 — Progressive tile generation with regeneration (SUPERSEDED)

Generate preliminary tile content as features arrive. If new features
arrive in a tile's cells later, mark the tile dirty and regenerate at
EOF.

**Pros:**
- Maximum parallelism — tile generation overlaps with ingestion.
- Useful when upstream has spatial locality.

**Cons:**
- Speculative execution adds invalidation tracking and dirty-flag
  management.
- Tiles might be generated multiple times (wasteful).
- Cannot guarantee tile-content correctness without either a completion
  signal or a guarantee about input ordering.

**Verdict:** superseded by the morton-ordered input approach
([04-morton-ordered-input](./04-morton-ordered-input.md)), which
provides a clean, non-speculative completion signal via spatial
ordering.

## Alternative 7 — Progressive per-feature conversion (ACCEPTED as core)

Decompose tile conversion into per-feature work (triangulation, format
encoding) and per-tile assembly (GLB packing, merging). The per-feature
work is independent across features and can happen during streaming.
Per-tile assembly waits for the tile's feature set to be complete.

**Pros:**
- The expensive per-feature work (triangulation) happens progressively
  during streaming, not at EOF.
- Post-EOF (or post-quadrant-completion) assembly is fast — mostly
  buffer concatenation and index rebasing.
- For non-triangulated formats, per-feature work is trivial (parse +
  spill bytes) — the architecture degenerates gracefully to the base
  case.

**Cons:**
- Requires an in-process Rust library for geometry conversion (replaces
  the geof subprocess).
- The per-feature intermediate format is format-specific.

**Verdict:** accepted as the core progressive mechanism. Requires the
native Rust geof library ([05-geof-library](./05-geof-library.md)).

## Accepted set

The final architecture combines three accepted alternatives:

1. **Alternative 1** (pre-declared extent + per-cell spill files) as the
   base streaming mechanism.
2. **Alternative 2** (two-pass with monolithic spill) as an automatic
   fallback when no extent is supplied.
3. **Alternative 7** (progressive per-feature conversion) as the core
   mechanism for pushing expensive work earlier.

These are then composed with morton-ordered input
([04-morton-ordered-input](./04-morton-ordered-input.md)) and the geof
library ([05-geof-library](./05-geof-library.md)) to produce the
complete architecture
([06-final-architecture](./06-final-architecture.md)).
