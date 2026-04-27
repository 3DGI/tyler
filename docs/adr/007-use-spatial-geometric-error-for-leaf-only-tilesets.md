# Use Spatial Geometric Error for Leaf-Only Tilesets

## Status

Accepted

## Context

Tyler currently writes full-detail GLB content only on leaf tiles. Internal
tiles are traversal and culling nodes, not simplified renderable versions of the
same source geometry. In this model, `geometricError` does not measure mesh
decimation error on internal tiles; it controls when the viewer should refine
toward the detailed leaves.

The previous explicit-tileset calculation mixed two heuristics:

- internal tiles used a configured above-leaf error and grid cell size
- leaf tiles used a small fraction of tile bbox diagonal, with a positive
  minimum

This could produce non-monotonic trees where a child had larger
`geometricError` than its parent. It also gave full-detail leaves positive error,
even though there are no children below them to refine to.

The quadtree split criterion is based on vertex count. That is useful for
controlling tile cost, but vertex count is not itself a spatial error metric. A
small dense tile can have many vertices without needing to refine earlier than a
larger sparse tile.

## Decision

For leaf-only 3D Tiles output, Tyler will compute internal-tile
`geometricError` from tile footprint size:

```text
geometricError = tile_width * geometric_error_factor
```

Leaf tiles use:

```text
geometricError = 0.0
```

The CLI option is `--geometric-error-factor` and keeps `-e` as the short option.

The default factor is `0.024`, which preserves roughly the old default behavior
near the deepest parent for the default 250 m grid cell size:

```text
500 m parent tile * 0.024 = 12 m geometricError
```

The top-level tileset `geometricError` is set to the root tile
`geometricError`.

## Consequences

Good:

- explicit trees have a simple monotonic error ladder based on spatial size
- full-detail leaves use the standard `0.0` geometric error
- vertex count continues to control tile subdivision and tile cost, while
  `geometricError` controls spatial refinement
- implicit tilesets inherit a root error that is tied to root footprint size,
  so Cesium's per-level implicit error reduction starts from a spatial value

Trade-offs:

- the factor is still empirical and may need dataset/viewer tuning
- dense urban tiles no longer refine earlier solely because they contain more
  vertices; if this is desired later, it should be a separate cost-aware loading
  policy rather than the primary geometric-error metric
- when Tyler adds true simplified parent content, geometric error should be
  revisited and based on actual simplification error, such as removed geometry,
  mesh decimation, or texture error

## Validation

- CLI test verifies `--geometric-error-factor` parsing
- 3D Tiles format test verifies generated geometric errors decrease toward
  zero-valued leaves
- full `cargo test` passes
