# Rectangular voxel grid evaluation

These measurements describe the native backend at commit `a496d7f`. The current
backend is Rust; updated comparisons are in the
[Rust remesher evaluation](rust-voxel-remeshing-evaluation.md).

Evaluated on 2026-10-11 using the release `cjconvert`, with whole objects at LoD
2.2 from `/home/ravi/sync/tyler-lods/CitySJON-Gebouwen`. Inputs and encoding settings
match the [original evaluation](voxel-remeshing-evaluation.md): flat normals,
quantization, meshopt compression, retained source attributes, and no clipping.
The earlier measurements at commit `b412366` remain the resolution-based baseline.

The native fork now accepts a voxel edge length in model units. Each axis uses
`max(4, ceil(axis_extent / voxel_size) + 3)` cells, including empty padding.
Grid entries are 32-bit indices within occupied rows; row references use `size_t`,
and voxel coordinates are stored separately. Neither the old 256-cell resolution
limit nor the 10-bit coordinate packing limits the requested quality.

`cjconvert` uses `voxel_size = requested_error / (2 * sqrt(3))`, rounding to a
float toward finer quality. Voxel reconstruction therefore has half the requested
error as its diagonal allowance, independent of tile extent. Absolute
simplification uses the other half. Encoding can increase the final estimate.
Grids that exceed float-coordinate precision or byte-addressing limits fail
instead of increasing the voxel size. Available memory also limits dense grids.

## Export measurements

Each row is one run, including parsing, triangulation, remeshing, simplification,
feature attribution, metadata, and GLB encoding. Wall time and child-process peak
RSS are indicative. MB is decimal; memory is MiB.

| Dataset suffix | Requested m | Grid X×Y×Z | Voxel m | Triangles | GLB MB | Seconds | Peak RSS MiB | Estimated error m |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 110000_492000 | source | — | — | 3,526 | 0.075 | 0.011 | 17.5 | — |
| 110000_492000 | 5 | 1421×1327×24 | 1.443 | 474 | 0.038 | 0.337 | 202.4 | 5.828 |
| 110000_492000 | 20 | 358×334×9 | 5.774 | 60 | 0.030 | 0.024 | 22.1 | 20.767 |
| 110000_492000 | 40 | 181×169×6 | 11.547 | 28 | 0.029 | 0.015 | 18.1 | 40.307 |
| 174000_440000 | source | — | — | 39,229 | 0.692 | 0.076 | 44.0 | — |
| 174000_440000 | 5 | 1405×1400×44 | 1.443 | 3,762 | 0.245 | 0.592 | 391.4 | 5.733 |
| 174000_440000 | 20 | 354×353×14 | 5.774 | 218 | 0.171 | 0.107 | 59.3 | 20.634 |
| 174000_440000 | 40 | 179×178×9 | 11.547 | 28 | 0.167 | 0.093 | 58.0 | 39.863 |
| 120000_484000 | source | — | — | 775,314 | 13.802 | 1.664 | 577.2 | — |
| 120000_484000 | 5 | 1408×1567×65 | 1.443 | 55,292 | 3.938 | 8.266 | 1490.8 | 5.815 |
| 120000_484000 | 20 | 355×394×19 | 5.774 | 7,408 | 2.915 | 2.201 | 830.2 | 20.793 |
| 120000_484000 | 40 | 179×199×11 | 11.547 | 1,596 | 2.778 | 2.011 | 804.7 | 40.541 |

On the largest sample, 5 m now uses a `1408×1567×65` grid with 1.443 m voxels.
The old implementation clamped to `256³`, using approximately 8.887 m voxels and
reporting an 18.685 m estimate. The new estimate is 5.815 m, including encoding.
Peak RSS rises from 814.9 MiB to 1490.8 MiB, and runtime from 2.313 s to 8.266 s:
the finer requested quality creates 5,270,688 intermediate triangles before
simplification. The dense grid itself uses approximately 547 MiB.

At 20 m, the same sample uses a shallow `355×394×19` grid and runs in 2.201 s,
with 830.2 MiB peak RSS. Its 7,408 triangles produce a 2.915 MB GLB, compared with
775,314 triangles and 13.802 MB for the source. Rectangular sizing avoids spending
cells on empty height, while finer voxel sizes can still need substantial memory.

The sparse eight-family subset illustrates that cost: a 2 m budget creates a
`1944×2468×30` dense grid and uses 606.7 MiB despite having only 943 source
triangles. Grouping fewer whole objects or a future sparse grid could reduce
this cost.

## Surface-distance observations

Both directions use 1,024 area-weighted samples with NumPy `default_rng(271828)`
and exact point-to-triangle distance via a BVH. Both meshes are decoded from their
actual compressed, quantized GLBs with node transforms applied. The reference is
the selected source GLB, including internal interfaces, rather than unencoded
CityJSON. Sampled maxima are observations, not certified Hausdorff bounds.

| Dataset suffix | Requested m | Source → proxy max / RMS m | Proxy → source max / RMS m | Represented owners / retained rows |
|---|---:|---:|---:|---:|
| 110000_492000 | 5 | 14.807 / 1.321 | 9.597 / 0.834 | 43 / 47 |
| 110000_492000 | 20 | 892.404 / 82.579 | 6.538 / 2.448 | 11 / 47 |
| 110000_492000 | 40 | 1958.334 / 567.340 | 11.477 / 4.389 | 5 / 47 |
| 174000_440000 | 5 | 24.305 / 1.690 | 2.870 / 0.613 | 501 / 630 |
| 174000_440000 | 20 | 1192.135 / 227.674 | 10.114 / 2.266 | 98 / 630 |
| 174000_440000 | 40 | 2370.985 / 652.026 | 14.629 / 3.817 | 17 / 630 |
| 120000_484000 | 5 | 44.288 / 2.289 | 3.309 / 0.560 | 9422 / 11297 |
| 120000_484000 | 20 | 80.093 / 5.925 | 11.143 / 2.125 | 4128 / 11297 |
| 120000_484000 | 40 | 159.483 / 15.310 | 17.493 / 3.954 | 961 / 11297 |

Every decoded triangle retains a single owner across its three corners, all IDs
are in range, and metadata row counts match the source. Missing owner IDs alone
do not establish missing geometry, since attribution is approximate near merged
surfaces. Large forward distances do establish loss of source coverage.

The finer voxel size removes the clamping problem but does not guarantee smaller
sampled maxima. On the largest 5 m export the sampled forward maximum is 44.288 m
against a 5.815 m estimate. At 20 m, sparse datasets still lose isolated buildings,
with errors hundreds of metres beyond the estimate. This report must therefore
remain a quality-selection estimate: parent-tile `geometricError` still needs
component preservation or conservative surface-distance validation, followed by
monotonic propagation to parents.

The existing dense eight-family subset retains all eight owners at 0.5 m, using
`74×130×51` cells and producing 38 triangles from 96. Its sampled forward / reverse
maxima are 1.695 / 0.068 m. At 5 m it produces four triangles with four owners and
a 9.735 m forward maximum. Full subset measurements are included in the JSON.

## Validation and reproduction

`just ci` passes, as do the standalone vendored crate's strict library clippy
check and all 11 tests. Existing wrapper checks now exercise fully occupied rows
beyond byte-index capacity and coordinates beyond 1023 along each axis, for both
solid and shell reconstruction. Existing GLB checks cover widely separated boxes
with a grid axis above 1024, all encoding combinations, and regenerated normals.
A local C++ experiment with AddressSanitizer, UndefinedBehaviorSanitizer and float
cast checks passed for elongated closed boxes along each axis, all four option
combinations, and injected allocation failures with allocation cleanup.

Reproduce an export from the repository root:

```sh
cargo build --release -p cityjson-convert --bin cjconvert
target/release/cjconvert \
  /home/ravi/sync/tyler-lods/CitySJON-Gebouwen/buildings_2025_120000_484000.city.json \
  --lod 2.2 --remesh-error 5 --output /tmp/proxy.glb
```

Full measurements are in [performance JSON](rectangular-voxel-grid-performance.json)
and [distance JSON](rectangular-voxel-grid-distances.json). Dense-subset RSS includes
the NumPy launcher footprint and is not comparable to the full-file measurements.
Input SHA-256 checksums are recorded in the original evaluation.
