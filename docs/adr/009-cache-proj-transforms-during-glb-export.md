# Cache PROJ Transforms During GLB Export

## Status

Proposed

## Context

After the grid indexing and feature-preparation optimizations, large
`BuildingPart` exports still underused CPU during the GLB tile writing phase.
Rayon had enough tile jobs, and `perf sched` showed that Linux scheduler delay
was negligible:

- `tyler` accumulated about 328 seconds of CPU runtime across 9 threads
- average scheduler delay was about 0.010 ms
- maximum scheduler delay was about 1.476 ms

This ruled out tile-job starvation and OS scheduling as the main cause. The
process still only consumed about 3.5 to 3.6 CPU cores on average with
`RAYON_NUM_THREADS=8`, so the missing time had to be blocked/off-CPU time inside
the tile export work.

`strace -f -c` then showed that the missing time was mostly lock wait:

```text
85.73%  269.944s  futex
4.54%    14.305s  lseek
4.34%    13.678s  read
2.10%     6.602s  pread64
```

A futex profile showed almost all futex traffic on one address, but the first
capture grouped by futex address rather than by owner. A focused userspace
stack capture of glibc lock waits identified the owner as PROJ's internal
SQLite-backed authority database:

```text
__lll_lock_wait
sqlite3_column_text
osgeo::proj::io::DatabaseContext::Private::run
osgeo::proj::io::AuthorityFactory::...
osgeo::proj::operation::CoordinateOperationFactory::...
```

The repeated construction happened inside the GLB writer:

- `CoordinateTransform::from_model` created a source-to-ECEF transform with
  `Proj::new_known_crs(source_crs, "EPSG:4978", None)`
- `ClipVolume::geographic_region` optionally created a source-to-geographic
  transform with `Proj::new_known_crs(source_crs, "EPSG:4979", None)`

Both were called per tile. For 598 tiles, this caused many parallel requests for
the same CRS operations. PROJ responded by repeatedly querying its internal
authority database and serializing on its own SQLite/glibc locks.

## Decision

`cityjson-convert` will cache constructed PROJ transforms per worker thread
during GLB export.

The GLB writer now stores lightweight transform specifications in tile-local
state rather than owning `Proj` values directly:

- `CachedProjTransform` stores the `(source_crs, target_crs)` key.
- A thread-local `HashMap<(String, String), Proj>` stores the actual PROJ
  objects.
- `CachedProjTransform::convert` creates the PROJ object on first use in that
  worker thread, then reuses it for subsequent vertex conversions.

This applies to:

- source CRS to ECEF (`EPSG:4978`) for ENU / ECEF-relative GLB placement
- source CRS to geographic (`EPSG:4979`) for geographic clipping

The cache is thread-local rather than global because `Proj` wraps raw PROJ
context/object pointers and should not be shared across threads without a
separate thread-safety audit. Thread-local caching matches the Rayon execution
model: each worker pays the expensive PROJ operation construction once per CRS
pair, then reuses it for many tiles.

No CLI behavior or output semantics change.

## Consequences

Good:

- repeated per-tile PROJ operation construction is avoided
- PROJ's internal SQLite authority database should no longer dominate futex
  wait time after the first few tiles per worker
- the fix is local to `cityjson-convert` and preserves the existing
  `convert_to_glb(&CityModel, path, &ExportOptions)` API
- transform construction errors still surface at conversion time with source
  and target CRS context

Trade-offs:

- transform creation is now lazy, so an invalid CRS may fail on first vertex
  conversion instead of at `CoordinateTransform::from_model`
- each worker keeps its own copy of each used transform, so memory usage grows
  with `RAYON_NUM_THREADS * number_of_CRS_pairs`
- a long-lived worker thread may retain cached PROJ objects until the thread
  exits
- this does not remove the cost of `proj_trans` itself; it only removes repeated
  `proj_create_crs_to_crs` / authority lookup overhead

## Rejected Alternatives

- Use a single global `Mutex<HashMap<...>>` of `Proj` values.
  This would likely move the bottleneck from PROJ's internal locks to Tyler's
  own lock, and sharing raw PROJ objects across threads would need an explicit
  safety audit.

- Pre-create transforms once in `src/main.rs` and pass them through
  `ExportOptions`.
  This would require widening public API and threading lifetime/ownership
  concerns through the caller. It is heavier than needed because the problem is
  local to GLB conversion.

- Create one transform per tile but serialize transform creation.
  This would reduce PROJ lock contention but still perform the expensive
  authority lookup hundreds of times.

- Replace PROJ with a hand-coded transform for the observed CRS pair.
  That may be faster for one dataset, but it would be brittle and would bypass
  PROJ's CRS semantics, grid handling, and future datasets.

## Validation Plan

Functional validation:

- run `cargo fmt --all --check`
- run `cargo check --workspace --all-targets --all-features`
- run `cargo test --workspace --all-targets --all-features`
- keep existing GLB writer tests for:
  - ENU placement
  - ECEF reprojection
  - geographic clipping
  - geographic clip intersection solving

Performance validation:

- rerun the Amsterdam subset export with `RAYON_NUM_THREADS=8`
- compare the `Converting and optimizing ... tiles` wall time before and after
- rerun `strace -f -c` and confirm futex time is materially reduced
- rerun the glibc lock-wait stack capture and confirm PROJ authority database
  stacks are no longer dominant
- compare the GLB tile export wall time before and after the cache change

Expected profiling change:

- far fewer stacks under `osgeo::proj::io::DatabaseContext::Private::run`
- far fewer `sqlite3_column_text` / PROJ authority factory lock waits
- no meaningful change to CityJSON decode, feature preparation, or meshopt
  costs

## Notes

The first measurement after this change should still show a small amount of
PROJ authority database work: each Rayon worker must build its own cached
source-to-ECEF transform, and clipping may add a cached source-to-geographic
transform. The goal is to reduce that cost from "once per tile" to "once per
worker per CRS pair".
