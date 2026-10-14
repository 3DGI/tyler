# Use Rust Voxel Remeshing for Optional GLB Proxies

## Status

Proposed. The converter implementation and evaluations described here are
implemented on the `grid-remesh-simp` branch, including remeshed explicit parent
content. Clipped geometry support and validated representation errors remain
future work.

## Related Commits

- `b412366`: add optional meshoptimizer voxel remeshing and simplification to
  `cjconvert` GLB output.
- `a496d7f`: replace the capped cubic grid with a rectangular grid and explicit
  voxel size.
- Subsequent branch work ports reconstruction to Rust, reduces memory and
  repeated work, and introduces shared-grid processing interfaces.
- Later branch work moves the Rust code into `cityjson-convert`, removes the
  vendored crate, and switches native algorithms to upstream master with
  meshoptimizer v1.2.

## Context

Tyler needs compact geometry for future 3D Tiles parent content. Simplifying
each building independently leaves a minimum geometry cost per object and
cannot combine nearby, disconnected buildings through ordinary edge collapse.
Combining their buffers does not create edges between them. Reconstruction can
instead replace the original connectivity and aggregate surfaces at a chosen
spatial resolution before simplification.

Meshoptimizer v1.3 introduced an experimental voxel remesher with optimized
vertex placement. Its initial integration exposed two constraints: a capped
cubic grid could silently coarsen the requested resolution, and finer grids
required substantial memory and produced large intermediate meshes. Shallow
building datasets also wasted cells when all three axes used the longest
extent.

The first supported case is complete objects, such as buildings assigned to
one tile, without clipping. Implementation starts in `cityjson-convert` so
geometry quality, metadata transfer, cost, and error reporting can be evaluated
for Tyler's remeshed parent content. The spatial traversal-error policy in
[ADR 007](007-use-spatial-geometric-error-for-leaf-only-tilesets.md) continues to
apply to Tyler's existing tilesets.

## Decision

### Optional converter API and CLI

Add `ExportOptions.voxel_remesh: Option<VoxelRemeshOptions>`, disabled by
default. `VoxelRemeshOptions.geometric_error` supplies a positive, finite
quality budget in Cartesian metres. Library callers prepare the geometry and
select source LoDs before export; remeshing combines all supplied triangles,
including different object types, in one reconstruction.

`cjconvert` exposes `--remesh-error <METRES>` for GLB output and an optional
`--lod <LOD>` for every output format. LoD selection matches one exact source
LoD and retains object relationships and geometry-less parents. Without `--lod`,
all source LoDs are exported, including when remeshing. Selection fails only if
no geometry matches across the input, including across CityJSONSeq features.
The CLI uses source-coordinate placement;
callers must supply metric Cartesian geometry rather than geographic degrees.

`convert_to_glb` and the GLB writer return `GlbExportReport` with output bytes
and an optional `RemeshReport`. The latter records requested error, effective
voxel size, three padded grid dimensions, error components, and input,
reconstructed, and final triangle counts. The CLI prints these diagnostics to
stderr.

Both Cartesian and geographic clipping options are rejected with remeshing.
Tyler maps its explicit content plan directly into this option: full-detail jobs
use `None`; simplified internal jobs use `Some(VoxelRemeshOptions)` with the
planned maximum error. `--3dtiles-max-geometric-error` sets the first internal
content level's budget, halved at each deeper level. Leaves and implicit content
remain full detail. Explicit content includes complete objects without clipping.

Remove `GeometryRepresentation` and its field from converter export options.
`voxel_remesh` is the single converter control for this geometry reduction.
Tyler retains its internal representation policy for planning, converting its
`f64` allowance to `f32` toward finer quality. Non-finite, non-positive,
overflowing, or underflowing allowances return errors instead of being clamped.

### Surface reconstruction followed by simplification

Use the meshoptimizer v1.3 reconstruction algorithm, ported to safe Rust:

1. Collect the prepared triangles in a local Cartesian frame, ordering them
   consistently by source object ID. Build one padded grid for the complete
   export.
2. Sample triangle interiors and edges. Sampling density depends on each
   original triangle's longest edge and the voxel size. The first pass marks
   sampled surface cells; the second accumulates weighted sample positions,
   triangle-plane quadrics, and half-cell occupancy flags.
3. Classify enclosed empty cells as solid interiors using outside propagation.
   The converter uses solid mode; interior cells influence connectivity without
   receiving output vertices.
4. Fit one candidate vertex per sampled surface cell. A regularized quadric
   solve keeps it near the contributing planes and sample average, with the
   solution clamped to its cell. Unsuitable solves retain the sample average.
5. Connect fitted vertices using lookup rules for overlapping `2 x 2 x 2` cell
   neighbourhoods. Half-cell occupancy and quadric costs resolve ambiguous
   cases. The result is indexed geometry with new connectivity, independent of
   source object ownership.
6. Weld coincident positions, then apply meshoptimizer's native
   position-updating quadric simplifier with absolute error.
   The target is at most four triangles for the whole export; error and topology
   constraints normally stop simplification much earlier.

Different objects can contribute to the same cell or become connected through
reconstruction. This permits later edge collapse across former object
boundaries. Components that remain disconnected are still separate for edge
collapse. Reconstruction is resolution-dependent and is not an exact solid
union; narrow gaps can close, interiors can fill, and small components can
disappear.

After simplification, assign each output triangle the feature ID of the source
triangle nearest its centroid, using exact point-to-triangle distance through a
BVH and source-order tie breaking. All three corners share that categorical ID.
Retain original metadata rows, including owners with no surviving geometry.
Regenerate normals and group materials by the attributed feature type.
Remeshed materials use backface culling. Existing colour, smooth-normal,
quantization, and compression options remain available.

### Explicit voxel size and compact rectangular storage

Accept a voxel edge length directly, with each allocation axis sized from its
own extent and empty padding. Remove the original 256-cell cap and packed
10-bit coordinate restriction; never increase the voxel size to fit a fixed
resolution. Invalid sizes, unrepresentable coordinates, index overflows, and
allocation failures return errors. Float half-cell precision and local axis
dimensions still impose a `2^23`-cell limit.

Keep the reconstruction algorithm while reducing its overhead:

- Build occupancy once and sample the source twice, removing repeated grid
  preparation for separate counting and output calls.
- Pack surface and solid-interior occupancy into 64-bit words, with 32-bit
  rank prefixes for surface-vertex lookup. These arrays use approximately
  20 bytes per 64 cells in solid mode, with X rounded to a word boundary.
- Propagate outside state across words and skip uniform spans during
  polygonization.
- Emit shared positions and triangle indices directly. Release occupancy,
  quadrics, cell provenance, and remapping buffers when their work is complete.

The core modules forbid unsafe Rust and check sizes before fallible allocations.
Storage and scans remain dense in the bounding volume, approximately
`O(volume / h^3)`. Sampling work remains approximately
`O(sum((longest_triangle_edge / h)^2))`, and intermediate surface geometry
typically grows with area divided by `h^2`. Sparse storage and a different
triangle voxelizer are deferred.

### Geometric-error accounting

For requested budget `E`, choose:

```text
h = E / (2 * sqrt(3))
reconstruction allowance = sqrt(3) * h
simplification tolerance = E / 2
estimated geometric error = reconstruction allowance
                          + simplifier-reported error
                          + float-rounding and quantization allowances
```

Round the voxel size toward finer quality when converting to `f32`. This keeps
the voxel-diagonal allowance within half the requested budget independently of
dataset extent. Encoding allowances can raise the final estimate above `E`.

This report is a quality-selection estimate, not a certified maximum
surface-distance bound. It does not account adequately for disappearing
components or all changes in source coverage. In the largest v1.3 5 m evaluation,
the estimate is 5.815 m while the sampled source-to-proxy maximum is 44.288 m.
Sparse datasets at coarser settings lose buildings hundreds of metres from the
remaining proxy. Missing owner IDs alone do not prove geometry loss because
attribution is approximate; measured source-to-proxy distances establish it.

Parent-tile error validation must address component coverage and establish a
conservative error against selected source geometry, including encoding. Area
samples in both directions are useful diagnostics but do not certify a maximum.
Propagate validated errors monotonically (`parent >= children`), use `REPLACE`
refinement, and generate proxies from the selected source geometry. Do not use
the current report directly as a parent tile's `geometricError`.

### Upstream native dependency and local Rust ownership

The initial implementation used a vendored `meshopt` 0.6.2 fork with native
meshoptimizer v1.3. Reconstruction was ported from v1.3 commit
`9e1f07b159d3cb777f1c67ed31fc11fd117986f4`; the native remesher and its bindings
were then removed. That fork is now removed entirely.

Keep the reconstruction and its existing checks in
`cityjson_convert::voxel_remesh`, together with the original
[MIT licence](../../cityjson-convert/src/voxel_remesh/LICENSE.md). This public
module owns `VoxelGrid`, `VoxelExtent`, `RemeshOptions`, `RemeshResult`, and
`remesh`; its error type is the converter's existing `anyhow::Result`.
The grid, reconstruction, and triangulation-table modules forbid unsafe code.
The same module provides the checked native position-update wrapper with
optional per-vertex lock flags.

Use upstream `meshopt` master as a git dependency pinned at
`3c8c9126cd98e70f34c6d9c54d02f41af6749406` (crate version 0.7.0). Its native
meshoptimizer submodule is v1.2, commit
`9d9890c73011d75920af614485296d1e03e95448`. Native simplification, remapping,
optimization, and compression remain C++; there is no local patch, native-source
copy, or custom binding generation. Cargo and Nix pin the upstream source.

v1.2 already exposes position-updating simplification and vertex locks. It lacks
`PreserveFolds`, introduced in
[v1.3](https://github.com/zeux/meshoptimizer/releases/tag/v1.3), so remove that
flag and use `ErrorAbsolute` alone. This changes simplification behaviour and
can reduce preservation of sharp creases in double-sided sheets. The earlier
v1.3 measurements remain historical; evaluate the v1.2 outputs separately rather
than claiming byte equivalence with them. See the
[upstream v1.2 evaluation](../performance/voxel-remeshing/upstream-meshopt-1.2-evaluation.md).

### Groundwork for clipped tile boundaries

Expose `VoxelGrid` with an explicit origin and voxel size, and `VoxelExtent`
with a half-open global processing-cell range. `RemeshResult` returns indexed
geometry, grid dimensions, and queryable global cell identities. Processing
extent is separate from any eventual output clipping boundary.

Sampling depends on the original triangle and lattice, not processing extent.
Out-of-domain samples are skipped rather than clamped. Shared regions require
identical original triangles and consistent accumulation order. Regional shell
reconstruction is supported; cropped solid domains are rejected because the
current interior fill assumes its domain boundary is outside. Solid mode
requires complete input coverage with padding.

The converter currently uses a centred grid in its local export frame and
passes no simplifier locks. It does not yet reconstruct adjacent tiles on a
shared lattice. Future integration needs a common Cartesian frame and origin,
buffered original geometry with matching boundary-cell contributions, valid
interior classification, consistent clipping intersections, locked boundary
vertices, and coordinated GLB quantization and placement. Different-resolution
neighbours need a transition policy as well. A shared grid or halo alone does
not guarantee matching boundaries: interior classification is nonlocal, and
simplification and encoding can change otherwise matching vertices.

## Alternatives Considered

These are candidates for evaluation, not implemented alternatives or measured
quality rankings on the branch datasets. QEM is an error metric; restricting
contractions to existing edges is a separate algorithm choice.

| Approach | Cross-object simplification | Benefits and limitations |
|---|---|---|
| Grid snapping, welding, and triangle cleanup | Vertices from different objects merge when they snap to the same grid position. | Cheap processing of existing vertices; snapping alone does not reduce triangle count. Welding and cleanup can remove triangles, but introduce grid distortion and depend strongly on source vertex placement. |
| Quadric vertex clustering | Vertices in one cell share a representative regardless of source object. | A sparse map of occupied cells and a quadric solve can provide approximately linear expected work over vertices and triangles, without volume filling or surface sampling. Placement can fit source planes, but topology can become non-manifold and large nearby faces with distant corners may remain separate. |
| QEM with non-edge pair contractions | Nearby vertices without a connecting edge become eligible for contraction. | Supports adaptive aggregation of disconnected components. Requires spatial candidate searches, collapse bookkeeping, and non-manifold handling; dense proximity candidates can be quadratic. It does not provide an exact union or a certified surface-distance bound. |
| Meshoptimizer `simplifySloppy` | Spatial clustering merges topologically disconnected features. | Available as a baseline without a new simplifier. It selects an existing representative vertex per cell using quadric error, rather than solving a new position. Quality is generally lower than regular simplification; its internally normalized grid has a 1024-position-per-axis limit and does not expose our shared origin and arbitrary voxel-size interface. |

[Open3D's vertex-clustering API](https://www.open3d.org/docs/release/python_api/open3d.geometry.TriangleMesh.html#open3d.geometry.TriangleMesh.simplify_vertex_clustering)
illustrates average and quadric representative placement. The
[original Garland-Heckbert QEM paper](https://www.cs.cmu.edu/~garland/Papers/quadrics.pdf)
explicitly supports aggregation through non-edge pair contraction, including
closely spaced disconnected cubes. Meshoptimizer documents the trade-offs of
[`simplifySloppy`](https://github.com/zeux/meshoptimizer/blob/v1.3/README.md#basic-simplification);
its grid and representative selection are visible in the
[implementation](https://github.com/zeux/meshoptimizer/blob/v1.3/src/simplifier.cpp).

Vertex clustering and pair contraction merge existing vertices rather than
sampling complete surfaces. They can miss nearby large walls when the input
has no nearby vertices, and do not automatically remove all internal or
overlapping faces. Nearest-grid snapping bounds vertex movement by
`sqrt(3) * h / 2`, but this does not bound source coverage after degenerate
triangles or entire components are removed.

Reconstruction remains the implemented baseline because it observes triangle
interiors and rebuilds connectivity across objects. Further comparisons should
measure runtime, peak memory, triangle count, GLB size, bidirectional surface
deviation, and disappearing components at comparable measured quality.

## Consequences

Good:

- GLB export can create compact proxies while retaining categorical metadata
  and the existing encoding options.
- Requested resolution is independent of a fixed grid cap, and rectangular
  sizing suits shallow building datasets.
- The Rust implementation reduces repeated work and occupancy/output storage
  while establishing explicit interfaces for future regional processing.

Trade-offs:

- Whole-object source coverage, certified parent geometric error, and clipped
  boundary preservation remain unresolved; the feature stays experimental.
- Fine cells still require dense volume scans and can create substantially more
  triangles than the source before simplification.
- Geometry ownership after merging is approximate, and retained metadata
  imposes a GLB size floor.
- Upstream native algorithms still require C++ build tooling. The v1.2 backend
  lacks fold-preserving simplification, while the local Rust port still needs
  maintenance. Existing source-coordinate float conversion before recentering
  also contributes precision loss; moving that conversion is separate work.

## Validation and Measurements

Workspace acceptance is `just ci`, including the relocated reconstruction and
vertex-lock checks. Strict workspace linting now covers those modules as well;
no new clippy allowances are added. Checks cover invalid inputs, wide grids
beyond the old packed-index limits, negative global cell coordinates, matching cell
positions across overlapping shell domains, cropped-solid rejection, vertex
locks, feature IDs, metadata, normals, encoding combinations, and clipping
rejection.

The existing content-plan check also verifies that full-detail jobs disable
remeshing and internal jobs enable it with their scheduled budgets. Fractional
`f64` budgets round toward finer quality; overflowing, underflowing, and NaN
budgets fail conversion. No new test functions are added for this wiring.
A local Tyler export of eight complete building families wrote four remeshed
internal tiles and four full-detail leaves without failures. The root proxy has
120 triangles at a 5 m allowance; deeper internal tiles use 2.5 m.

Before switching the native simplifier to v1.2, a local comparison against the
preserved native implementation exercised 48 elongated-box cases with different
axes, voxel sizes, and solid/shell/solve
options. Expanded triangle positions and ordering matched exactly. All 16
compared source/proxy GLBs from the building datasets are byte-identical across
the native rectangular-grid implementation and Rust port. This validates those
comparisons without claiming that every possible input is byte-identical or
that reconstruction error is bounded.

After the v1.2 switch, the reconstruction comparison still matches all 48 cases
exactly. All four unremeshed source GLBs are byte-identical, while the 12 proxies
change. The largest 5 m proxy now has 54,932 triangles and occupies 3,935,048
bytes, with 989.6 MiB peak RSS. Its 5.810 m estimate remains below the sampled
43.943 m source-to-proxy maximum; removing the fork does not resolve coverage
loss. Updated exports and bidirectional distance checks are recorded separately.

The original Rust-port measurements use whole objects at LoD 2.2 from
`/home/ravi/sync/tyler-lods/CitySJON-Gebouwen`, with quantization, compression,
flat normals, and retained metadata. Single-run timings and RSS are indicative:

- Largest 5 m export: peak RSS falls from 1490.8 MiB to 989.9 MiB (33.6%);
  runtime changes from 7.989 s to 7.759 s. It reconstructs 5,270,688 triangles
  from 775,314 and simplifies to 55,292, producing a 3,938,252-byte GLB.
- Medium 5 m export: RSS falls from 390.0 MiB to 78.1 MiB and runtime from
  0.581 s to 0.236 s.
- Sparse eight-family subset at 2 m: RSS falls from 606.4 MiB to 85.7 MiB.

The full reports, JSON measurements, and comparison image are consolidated in
Tyler's `docs/performance/voxel-remeshing/` directory:

- [Original resolution-based evaluation](../performance/voxel-remeshing/voxel-remeshing-evaluation.md)
- [Native rectangular-grid evaluation and surface distances](../performance/voxel-remeshing/rectangular-voxel-grid-evaluation.md)
- [Rust port evaluation and performance comparison](../performance/voxel-remeshing/rust-voxel-remeshing-evaluation.md)
- [Upstream v1.2 evaluation](../performance/voxel-remeshing/upstream-meshopt-1.2-evaluation.md)
