# Drop Empty Primitives After Clipping

## Status

Accepted

## Context

When tyler builds an implicit-tiled 3D Tiles output that includes non-unique
CityObject types (e.g. `Road`, `Railway`, `WaterBody`, `TINRelief`), some
content tiles failed GLB conversion with:

```
Tile X/Y/Z conversion failed: primitive bounds missing for non-empty mesh
```

The failed tiles were then dropped from the tileset by the failed-tile
pruning in `src/main.rs`, leaving holes in the output. The same input did
not fail with explicit tiling.

The bug is not in the grid-assignment / tile-routing logic — it lives
downstream in the GLB writer.

### How the empty primitive was produced

In `cityjson-convert/src/gltf_writer.rs`:

1. `MeshCollector` aggregates triangles into one primitive **per
   `feature_type`** (`primitives: BTreeMap<String, RawPrimitiveMesh>`).
2. `ProcessedScene::from_collector` calls `RawPrimitiveMesh::build()`, which
   clips each triangle against the active `ClipVolume`. A primitive can come
   out with zero vertices and zero indices when every triangle of that
   feature_type fell entirely outside the clip region.
3. All built primitives — including empties — were pushed into
   `ProcessedScene::primitives`.
4. The early-out in `write_glb` (`all primitives empty → write empty GLB`)
   only fires when every primitive is empty. If at least one feature_type
   still has geometry inside the tile and another has been fully clipped
   away, encoding proceeds.
5. Per-primitive encoding then calls `Bounds::from_vertices` in
   `ProcessedPrimitiveMesh::quantize` and
   `ProcessedScene::encode_primitives_raw`. Both return `None` for an empty
   vertex slice and bail with `primitive bounds missing for non-empty mesh`.

### Why only implicit tiling

`clip_geographic_region` is set only for implicit content tiles in
`src/main.rs`. Explicit tiling uses `clip_bbox` derived from the qtree node
bbox, which is fitted to cells that actually contain feature vertices, so
per-feature-type primitives almost never become empty after clipping.

### Why only non-unique CityObject types

Buildings/BuildingParts use unique-assignment: each feature is centroid-routed
to a single content tile that is guaranteed to contain real geometry of that
feature. Non-unique types use bbox-overlap assignment, which routes a feature
to every content tile its bbox overlaps. A feature's bbox can overlap a tile
that the actual triangles never touch — typical for diagonal linear features
(roads, railways), whose axis-aligned bbox is much larger than the geometry
itself.

The failing combination is: implicit tiling + non-unique types + a tile that
receives a feature by bbox-overlap whose triangles all sit outside the tile,
while another feature_type in the same tile still has surviving triangles.

## Decision

Drop empty primitives in `ProcessedScene::from_collector` immediately after
`build()`, so they never reach the encoding path. The existing all-empty
early-out in `write_glb` still applies — if filtering removes everything, an
empty GLB is written as before.

This is a single, localised change in the GLB writer. Grid-assignment and
tile-routing semantics are untouched.

Concretely, in `cityjson-convert/src/gltf_writer.rs`, after the per-feature
`build()` call produces `BTreeMap<String, BuiltPrimitiveMesh>`:

```rust
primitives.retain(|_, primitive| {
    !(primitive.vertices.is_empty() && primitive.indices.is_empty())
});
```

`reorder_features_by_type` tolerates dropped feature_types — it sorts
features and remaps vertex `feature_id`s by index, with no requirement that
every feature_type appear as a primitive.

## Consequences

Good:

- implicit-tiled output no longer drops content tiles that mix fully-clipped
  and surviving feature types
- the GLB writer enforces an invariant the encoder already assumed: every
  primitive in `ProcessedScene::primitives` has at least one vertex
- no change to grid-assignment, tile-routing, or clip-volume logic; the fix
  is local to the GLB writer

Trade-offs:

- a feature_type that is wholly clipped away is silently absent from the
  tile's GLB, with no warning. This matches the existing behaviour for
  fully-empty tiles, but it does mean clipping-induced data loss is no
  longer surfaced through the failed-tile log line. If we want visibility,
  a separate debug log when a primitive is dropped would be the
  follow-up.

## Validation

- new regression test
  `convert_to_glb_drops_feature_types_fully_outside_clip_region` in
  `cityjson-convert/tests/gltf_writer_geometry.rs`: a `TINRelief` triangle
  inside the clip region and a `Road` triangle entirely outside; asserts
  the conversion succeeds and the output mesh contains exactly one
  primitive
- full `cityjson-convert` test suite: 16/16 passing
- full workspace `cargo test --workspace`: all suites green
