# 03 — Progressive Tile Generation

This document analyzes how much of tyler's conversion work can be
performed progressively — i.e., before the end of the input stream — and
establishes the theoretical basis for the architecture's progressive
behavior.

## The question

Tyler's current pipeline pushes all tile conversion to the end: every
feature must be indexed before the quadtree is built, the quadtree must
be built before the tileset is generated, and the tileset must exist
before tiles are exported. Can any of this work be moved before EOF?

## The monotonicity property

The quadtree merge decision (`src/spatial_structs.rs:111-123`) is:

```rust
let sum_items: usize = tiles.iter().map(|t| t.nr_items).sum();
if sum_items <= limit {
    // MERGE: 4 children collapse into one leaf tile
} else {
    // SPLIT: 4 children remain separate tiles
}
```

Vertex counts (or object counts, depending on `QuadTreeCriteria`) only
**increase** as features arrive. This gives us a key invariant:

> **Monotonicity:** once `sum_items > capacity` for any group of 4
> siblings at any quadtree level, that decision is permanent. The group
> will never merge; its children's tile boundaries are final.

Cascading upward: once a group is permanently split, the parent group's
sum includes this group's sum, so the parent is also >= capacity and
also permanently split. Sealing propagates monotonically up the tree.

## The three layers of conversion work

Conversion is not monolithic. It decomposes into three layers with
different finalization requirements:

| Layer | Operation | Can finalize before EOF? |
|---|---|---|
| 1. Per-feature | Parse, triangulate, compute normals, extract geometry | **Yes** — each feature is independent |
| 2. Per-tile boundary | Determine which cells form which tile | **Partially** — sealed groups are final |
| 3. Per-tile assembly | Pack features into output format | **Only after tile's feature set is complete** |

### Layer 1 is fully progressive

Every feature is triangulated (or format-encoded) independently. There
is no dependency between features — feature A's triangulation does not
depend on feature B. This means:

- During streaming, each feature can be processed immediately on arrival.
- The expensive computational-geometry work (triangulation, normal
  computation) happens once, distributed over the full runtime, with no
  batching required.
- The result is stored in a per-cell spill file (in a format-specific
  intermediate representation) and the raw CityJSON can be discarded.

### Layer 2 is partially progressive via sealing

A leaf tile's boundary becomes final when:

1. **Capacity seal:** the tile's cell (or the sum of 4 sibling cells)
   exceeds capacity. At this point, the tile will not be merged with
   siblings; its boundary is fixed.
2. **Quadrant completion:** all 4 cells in the tile's quadrant have
   been fully populated (no more features will arrive). The merge
   decision is final.

Condition 1 can be detected during streaming. Condition 2 requires
either EOF or an ordering guarantee on the input
([04-morton-ordered-input](./04-morton-ordered-input.md)).

### Layer 3 blocks on both (a) boundary finalization and (b) feature-set completeness

Even when a tile's boundary is known, assembly (packing all features
into the output format) cannot begin until no more features will be
added to the tile's cells. Without ordering or signals, this is only
guaranteed at EOF.

**However:** once a tile's feature set is complete, assembly is cheap
if layer 1 has been done progressively. Assembly for GLB becomes:
- Read pre-triangulated meshes from cell spill files
- Concatenate vertex buffers
- Rebase index buffers (offset by cumulative vertex count)
- Write glTF JSON header + binary buffer

This is O(n) in the total vertex count, with no computational geometry —
just memory layout work.

## The progressive work budget

For a Netherlands-scale run (10M buildings, 250m cells, 42K vertex
capacity):

- **Layer 1 (progressive):** triangulation of ~10M buildings is the
  dominant CPU cost. Happens during streaming. **Push earlier.**
- **Layer 2 (mostly progressive):** most dense urban cells individually
  exceed 42K vertices, so they seal immediately when the threshold is
  crossed. Sparse rural cells remain undecided until EOF (or until the
  morton frontier passes).
- **Layer 3 (fast, post-seal):** once a leaf is sealed and its feature
  set is complete, assembly is a memory-layout operation. For dense
  cells this can happen well before EOF (if ordering is guaranteed).

The runtime profile shifts from "ingest fast, convert slow at end" to
"ingest and convert continuously, finalize sparse regions at end."

## Progressive sealing algorithm

The sealing logic maintains per-quadtree-node "sealed" flags. When a
feature is added to a cell:

```
on_feature_added(cell_id, added_vertices):
    cell.nr_vertices += added_vertices
    // Walk up the tree, sealing groups that newly exceed the limit
    for level in (0..max_level):
        group_id = cell_to_group(cell_id, level)
        group_sum = sum_of_cells_in_group(group_id, level)
        if group_sum > capacity and !group.sealed:
            group.sealed = true
            // All direct children of this group have fixed boundaries
            for child in group.children:
                if child.is_leaf_tile():
                    submit_to_assembly_queue(child)
```

Once a group is sealed, the propagation can stop at the next level that
is already sealed (since monotonicity guarantees sealing cascades
upward).

## When a sealed leaf is also assemblable

A sealed leaf has a final boundary, but its feature set may still grow.
Assembly can begin only when we know no more features will arrive. The
three ways to know this:

1. **EOF** — always works, but forces all assembly to the end.
2. **Morton-ordered input + low-water mark** — the stream's spatial
   position guarantees that cells to the "left" of the frontier are
   complete. Clean, implicit.
3. **Explicit `+RegionComplete` control signal** — upstream declares a
   bbox complete. Explicit, optional.

The recommended design uses (2) for progressive assembly, falls back to
(1) for sparse regions or when ordering is imperfect, and accepts (3)
as an optional future enhancement.

## What cannot be made progressive

- `tileset.json` — describes the complete tile hierarchy. Requires the
  final quadtree. **Must be written at EOF.** Small JSON file; fast to
  generate.
- Sparse-region tiles where the quadrant sum stays below capacity —
  their merge decision depends on the final vertex counts. Cannot be
  finalized until EOF (or until the morton frontier passes all 4
  siblings).
- Cross-feature optimizations in the output format (e.g., cross-feature
  vertex deduplication in GLB) — if the converter performs these, per-
  feature encoding is not possible. The design assumes per-feature work
  is truly independent.

## Implications for the design

1. The streaming ingestion loop performs layer-1 work (per-feature
   encoding) inline.
2. A per-cell spill file holds the format-specific intermediate
   representation, not the raw CityJSON bytes (for triangulated formats;
   for non-triangulated formats, the spill may be raw bytes).
3. A live quadtree tracks the sealing state; newly-sealed leaves trigger
   assembly.
4. Assembly workers (Rayon pool) consume sealed leaves, perform layer-3
   work, and write output files continuously during ingestion.
5. At EOF: finalize remaining quadrants, assemble remaining tiles, write
   `tileset.json`.
