# Use Spatial Traversal Errors and Requested Content Allowances

## Status

Accepted; updated by geometric-error plans 02 and 04.

## Context

Vertex count controls subdivision and tile cost, but does not measure spatial
refinement or geometry simplification. A tile is a patch of the dataset; some
patches contain geometry and others direct the viewer toward smaller children.
Larger error numbers generally cause entry and refinement farther away with
the same viewer settings.

## Decision

Explicit tiling uses three independent settings:

| Setting | Default | Responsibility |
| --- | ---: | --- |
| `--3dtiles-max-content-width` | 4000 metres | Largest internal patch eligible for simplified content |
| `--3dtiles-max-geometric-error` | 5 metres | First eligible content level's simplification allowance |
| `--3dtiles-geometric-error-factor` | 0.024 | Spatial traversal and tileset entry error multiplier |

Non-empty internal patches at or below the width limit receive content.
Original leaves always receive full-detail geometry, even when wider than the
limit. The first content level receives the configured maximum allowance;
each deeper level receives half, and original leaves use zero. Moving the
width limit restarts this sequence at the newly selected first content level.
Widths usually halve, so different width limits can select the same level.

```text
Starting content error = requested simplification allowance
Starting traversal error = factor × patch width in metres
Final tile error = max(starting error, largest surviving child error)
Empty tile error = 0

Non-empty tileset error = max(factor × enclosing dataset width, 2 × final root error)
Empty tileset error = 0
```

The tile patch width is measured before bounds expand around complete objects.
The enclosing dataset width includes that expansion. External tilesets use the
same separate entry rule. A single full-detail root has zero refinement error
but retains a size-based tileset entry error. The multiplier 2 remains fixed;
both spatial estimates are viewing rules rather than measured geometry errors.

Store starting traversal and representation errors separately so failed exports
can be removed and errors recalculated. Failed parent content falls back to
spatial error when successful children survive. Successful simplified content
keeps its allowance even when all children fail. Store enclosing widths and
pass the same configured factor to generation, pruning and external splitting.
Internal estimates and the factor are excluded from serialized tileset JSON.

The factor applies to explicit, implicit, debug and tileset-only output. Implicit
tiling retains its spatial calculation and viewer-driven per-level reduction;
content width and simplification allowance remain explicit-only settings.

The CLI accepts finite, non-negative factors, including zero for experiments,
and rejects negative values, NaN and infinity before planning or exporting.
Calculated errors that overflow to infinity return errors. Zero removes the
size contribution, but child errors can remain positive; an entirely zero-error
dataset may fail to display or reach its content in the viewer.

## Consequences

- Parents retain errors at least as large as surviving children.
- Geometry placement, requested simplification and spatial viewing policy can
  be adjusted independently with the same executable and separate outputs.
- Default 0.024 retains plan 02's calculations; custom factors affect only
  spatial estimates. Compare outputs with identical camera and viewer settings.
- The converter currently exports full-detail geometry for requested simplified
  representations. Allowances support viewing experiments now; mesh reduction
  and measured representation errors remain future converter work.

## Validation

Existing active tests cover CLI validation, custom factors, content schedules,
failed exports, external entry errors, full-detail roots, implicit and debug
output, overflow rejection and two-factor CLI comparisons. `just ci` is the
acceptance check, without additional Clippy allowances.
