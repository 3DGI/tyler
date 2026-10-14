# Upstream meshopt v1.2 evaluation

Evaluated on 2026-10-12 after removing the vendored crate. Native algorithms use
upstream `meshopt` master, pinned at `3c8c9126cd98e70f34c6d9c54d02f41af6749406`
(crate version 0.7.0). Its meshoptimizer submodule is v1.2,
`9d9890c73011d75920af614485296d1e03e95448`.

The Rust reconstruction, explicit voxel size, rectangular packed grid, indexed
output, processing extents, and global cell identities now live in
`cityjson_convert::voxel_remesh`. Position-updating simplification still uses
upstream FFI with optional vertex locks. v1.2 lacks `PreserveFolds`, so converter
simplification uses `ErrorAbsolute` alone. This is a change in the native
simplifier and its options; the reconstruction algorithm is unchanged.

## Export comparison

Inputs and subset selection match the [previous Rust evaluation](rust-voxel-remeshing-evaluation.md):
complete building families at LoD 2.2, Cartesian source coordinates, flat normals,
quantization, compression, and retained metadata. The baseline consists of the
saved v1.3 outputs measured at commit `1439eaa`. All four unremeshed source GLBs
are byte-identical, while all 12 remeshed GLBs differ. Input and intermediate
reconstruction triangle counts match for all cases.

Current timings and peak child RSS are single-run observations, including
parsing through encoding. They are not statistical benchmarks; older v1.3
timings are available in the previous report.

| Dataset suffix | Requested m | v1.3 final triangles | v1.2 final triangles | v1.2 GLB bytes | v1.2 seconds | v1.2 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| eight-buildings | 2 | 280 | 272 | 24,344 | 0.723 | 85.2 |
| eight-buildings | 5 | 114 | 124 | 21,948 | 0.107 | 23.8 |
| eight-buildings | 10 | 58 | 58 | 20,784 | 0.031 | 17.0 |
| 110000_492000 | 5 | 474 | 488 | 38,468 | 0.177 | 38.7 |
| 110000_492000 | 20 | 60 | 62 | 29,972 | 0.022 | 18.0 |
| 110000_492000 | 40 | 28 | 26 | 29,192 | 0.016 | 17.2 |
| 174000_440000 | 5 | 3,762 | 3,756 | 245,068 | 0.230 | 77.5 |
| 174000_440000 | 20 | 218 | 208 | 170,564 | 0.102 | 57.9 |
| 174000_440000 | 40 | 28 | 18 | 166,396 | 0.098 | 58.1 |
| 120000_484000 | 5 | 55,292 | 54,932 | 3,935,048 | 7.521 | 989.6 |
| 120000_484000 | 20 | 7,408 | 7,200 | 2,910,024 | 2.235 | 780.3 |
| 120000_484000 | 40 | 1,596 | 1,572 | 2,777,352 | 2.046 | 767.8 |

The largest 5 m export produces 54,932 triangles and a 3,935,048-byte GLB,
compared with 55,292 triangles and 3,938,252 bytes previously. It takes 7.521 s
and 989.6 MiB, compared with the earlier 7.759 s and 989.9 MiB. These timings do
not establish a stable speed improvement. The 5,270,688-triangle reconstruction
intermediate is unchanged, so the main fine-grid cost remains.

## Surface-distance checks

Use the previous methodology: 1,024 area-weighted samples per direction,
NumPy `default_rng(271828)`, and exact nearest point-to-triangle distances through
a BVH. Measurements use decoded, quantized/compressed GLBs with node transforms
applied; reference geometry includes source internal interfaces. Triangle IDs
are categorical, in range, and identical across each triangle's corners.
Metadata row counts match every source export. Represented owner counts do not
by themselves measure lost geometry.

The largest dataset illustrates the change:

| Requested m | v1.3 source-to-proxy max m | v1.2 source-to-proxy max m | v1.3 proxy-to-source max m | v1.2 proxy-to-source max m |
|---|---:|---:|---:|---:|
| 5 | 44.288 | 43.943 | 3.309 | 4.993 |
| 20 | 80.093 | 82.585 | 11.143 | 11.032 |
| 40 | 159.483 | 184.541 | 17.493 | 15.713 |

The switch is not uniformly better or worse in these observations. For example,
the 40 m source-to-proxy maximum rises from 159.483 m to 184.541 m, while the
reverse maximum falls from 17.493 m to 15.713 m. At 5 m, the current reported
estimate is 5.810 m against a sampled source-to-proxy maximum of 43.943 m.
Sparse datasets still lose isolated components, with sampled distances exceeding
one kilometre at coarse settings. These sampled maxima are not certified
Hausdorff bounds, and the current estimate remains unsuitable for direct use as
a parent tile's `geometricError`. No new guarantee of thin-sheet crease
preservation or clipped boundaries follows from these checks.

## Reconstruction and workspace validation

The relocated Rust reconstruction was compared again against the preserved
native rectangular-grid implementation on the same 48 elongated-box cases:
three axes, four voxel sizes, and all solid/shell/solve option combinations.
Triangle counts, order, and positions match exactly, with maximum coordinate
difference zero; inspected solid meshes have positive signed volume. This local
experiment confirms the checked numeric-conversion changes made for workspace
linting on those inputs.

`just ci` passes with 175 workspace tests, including the existing remesher test
moved from the removed crate. Its checks still cover invalid inputs, wide grids,
regional shell consistency, cropped-solid rejection, and simplifier locks.
Strict workspace clippy covers the relocated code without new allowances.
The MIT licence and upstream attribution are retained beside the Rust modules.

The upstream meshopt Nix source hash was verified with its meshoptimizer
submodule. During this switch, Nix Cargo-dependency evaluation exposed a stale
`cityjson-index-0.9.0` output-hash entry while the lockfile used 0.10.0. A
follow-up packaging fix updates that entry and the CityJSON source hash to the
locked revision, and derives the package version from `Cargo.toml`.
`nix build --no-link .#tyler` now succeeds. Both packaged CLIs pass version and
help checks, and packaged `cjconvert` exports the eight-building 5 m proxy with
all eight metadata rows retained. `just ci` still passes with 175 tests.

Full exports, sampled distances, byte comparisons, and reconstruction comparisons
are recorded in [measurements JSON](upstream-meshopt-1.2-measurements.json).
