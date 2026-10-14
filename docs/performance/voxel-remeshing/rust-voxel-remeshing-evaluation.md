# Rust voxel remesher evaluation

These measurements describe the Rust reconstruction with the former native
meshoptimizer v1.3 backend at commit `1439eaa`. Reconstruction now lives in
`cityjson-convert`, and native algorithms use upstream `meshopt` master with
meshoptimizer v1.2, without `PreserveFolds`. See the
[upstream v1.2 evaluation](upstream-meshopt-1.2-evaluation.md) for current outputs.
The byte-equivalence and performance comparisons below describe the former
backends, not the subsequent v1.2 switch.

Evaluated on 2026-10-12. The native comparison uses the rectangular-grid backend
at commit `a496d7f`; the current backend ports its reconstruction to safe Rust.
Both use release `cjconvert`, whole objects selected at LoD 2.2, flat normals,
quantization, meshopt compression, and retained source metadata. Inputs and subset
selection match the [native evaluation](rectangular-voxel-grid-evaluation.md).
No input geometry was clipped.

The port implements the selected focused optimisations: one prepared grid,
packed occupancy, indexed output, and interfaces for consistent regional grids.
Sparse storage and a different triangle voxelizer remain separate work.

## Export measurements

Each row is one local run including parsing, triangulation, reconstruction,
simplification, attribution, metadata, and encoding. C++ exports ran first; the
final Rust exports ran after other validation processes had finished. These are
indicative wall times and child-process peak RSS, rather than repeated statistical
benchmarks. RSS is MiB. Source exports have comparable time and memory in both
versions; the largest source is about 577 MiB and 1.58 s.

| Dataset suffix | Requested m | C++ seconds | Rust seconds | C++ RSS MiB | Rust RSS MiB | RSS reduction |
|---|---:|---:|---:|---:|---:|---:|
| eight-buildings | 2 | 1.221 | 0.754 | 606.4 | 85.7 | 85.9% |
| eight-buildings | 5 | 0.122 | 0.111 | 62.9 | 24.3 | 61.4% |
| eight-buildings | 10 | 0.032 | 0.035 | 23.8 | 16.7 | 30.0% |
| 110000_492000 | 5 | 0.347 | 0.186 | 202.4 | 38.9 | 80.8% |
| 110000_492000 | 20 | 0.024 | 0.023 | 22.1 | 17.5 | 20.6% |
| 110000_492000 | 40 | 0.014 | 0.015 | 18.6 | 17.7 | 4.8% |
| 174000_440000 | 5 | 0.581 | 0.236 | 390.0 | 78.1 | 80.0% |
| 174000_440000 | 20 | 0.106 | 0.102 | 59.8 | 58.4 | 2.3% |
| 174000_440000 | 40 | 0.098 | 0.095 | 57.3 | 58.2 | -1.7% |
| 120000_484000 | 5 | 7.989 | 7.759 | 1490.8 | 989.9 | 33.6% |
| 120000_484000 | 20 | 2.189 | 2.261 | 830.2 | 780.4 | 6.0% |
| 120000_484000 | 40 | 2.036 | 2.069 | 805.3 | 768.3 | 4.6% |

The largest 5 m export falls from 1490.8 MiB to 989.9 MiB peak RSS (33.6% lower).
Time is approximately unchanged, at 7.989 s versus 7.759 s. The medium 5 m export
uses 78.1 MiB instead of 390.0 MiB, and takes 0.236 s instead of 0.581 s.
The sparse eight-family subset at 2 m drops from 606.4 MiB to 85.7 MiB.
Coarser exports have smaller grids, and much of their memory/time is elsewhere.
Small percentage increases on some rows should not be interpreted as stable
regressions from single runs.

All 16 source/proxy GLBs are byte-identical between the two backends on these
inputs. The largest 5 m proxy still contains 55,292 final triangles, occupies
3,938,252 bytes, and has a 5.815 m estimate. Its intermediate reconstruction has
5,270,688 triangles. The previous surface-distance observations therefore also
apply to these outputs: the estimate remains unsuitable as a certified parent-tile
`geometricError`, and component/source-coverage loss is still unresolved.

## Where memory and work were removed

- The native wrapper invoked reconstruction separately for counting and output.
  Rust builds occupancy once, then accumulates samples, solves positions, and
  emits an indexed mesh. Source sampling runs twice instead of three times.
- Surface/interior occupancy is packed into 64-bit row words, with 32-bit rank
  prefixes for surface vertex lookup. Solid grids use about 20 bytes per 64 cells,
  with X rounded up to a word, instead of four bytes per cell. On the largest
  fine grid this replaces roughly 547 MiB of grid entries with 43 MiB of masks
  and rank data, excluding row work buffers and surface records.
- Flood propagation operates on whole words; polygonization skips uniform word
  spans and only expands cubes that cross the surface.
- Output has one position per occupied surface cell and three indices per
  triangle. It does not materialize three duplicated positions per triangle.
  The converter welds only the indexed positions and releases source/remapping
  buffers before simplification. Compact provenance is released through
  `RemeshResult::into_mesh` when it is not needed.

The algorithm remains dense: volume storage/word scans still scale with the
bounding volume divided by `h^3`, although the constants are much smaller.
Sampling still costs approximately `sum((longest_triangle_edge / h)^2)`;
intermediate surface geometry typically grows with area divided by `h^2`.

A CPU-clock profile of the initial Rust version (734 samples at 99 Hz) attributed
25.9% of exclusive samples to the main native simplification routine, 8.7% to
collapse ranking, and 4.4% to adjacency updates. Sampling and half-cell coordinate
conversion accounted for approximately 23.3%. The port reduces grid overhead but
cannot avoid simplifying millions of fine-grid triangles with the current
reconstruction strategy. The final version inlines the hot coordinate conversion
and uses the checked 32-bit range for its float-to-integer operation.

## Groundwork for shared clipped boundaries

`VoxelGrid` supplies an explicit origin and voxel size in the input Cartesian
frame. `VoxelExtent` supplies a half-open global processing cell range with empty
padding. A smaller processing allocation does not shift the lattice or change
triangle sampling density. Unlike the native implementation, sampling rate is
not capped by the processing grid dimensions; it is determined by the original
triangle and voxel size. Out-of-domain samples are skipped rather than clamped.

`RemeshResult` exposes an indexed mesh and each vertex's global cell identity.
Regional callers must use the same original triangles, grid frame, and consistent
accumulation order for shared contributions. Reconstruction must occur in that
shared frame before per-GLB recentering. The centred convenience constructor is
for current standalone whole-object exports, not independent tiled alignment.

Solid reconstruction still assumes its processing boundary is outside. It rejects
extents that do not cover the full supplied input with padding. Regional shell
reconstruction is supported, providing the local processing interface without
silently treating cropped solids as valid. The native simplification wrapper now
accepts per-vertex lock flags so later clipping can preserve boundary vertices.

Clipped GLB remeshing remains rejected. Buffered input retrieval, solid boundary
classification, shared clipping intersections, locked simplification boundaries,
coordinated GLB quantization/placement, and different-resolution transitions are
future integration work. Shared cell identities alone do not solve those steps.

## Validation

`just ci` passes. The standalone fork's strict library clippy check and all 11
tests also pass. The existing remeshing test was refactored for indexed output,
wide rows and coordinates beyond 1023 along each axis, negative global cell
coordinates, exact shared cell positions across cropped shell domains, invalid
grid ranges, cropped-solid rejection, and locked simplification vertices.
Existing GLB checks still cover metadata, one owner per triangle, encoding modes,
normal generation, and clipping rejection. Crate packaging includes the new nested
Rust modules and excludes the removed C++ remesher.

A local comparison against the preserved native implementation exercised 48
elongated-box cases: three long axes, voxel sizes 0.2/0.7/1.43/3.0, and all four
solid/shell/solve option combinations. Triangle counts, order, and positions
matched exactly; inspected solid meshes had positive signed volume. These are
localized port-validation experiments, not a retained native fallback or tests
for the superseded API.

Full paired runs are in [performance JSON](rust-voxel-remeshing-performance.json).
Input SHA-256 checksums are recorded in the original evaluation.
