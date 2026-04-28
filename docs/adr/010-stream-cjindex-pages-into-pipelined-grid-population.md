# Stream cjindex Pages into a Pipelined Grid Population

## Status

Proposed

## Related Commits

- `36b2071` Parallelise in chunks for extent calc and vertex counting

## Context

After ADR 005 reduced full extent scans and ADR 008 cut allocations in the
grid-population path, the extent and vertex-counting phases were further
parallelised by chunking each `cjindex` page across rayon workers (commit
`36b2071`). On a 32-core EPYC 7542 running with 20 rayon workers, that change
roughly halved the wall time of both phases, but `cpustat` showed a pronounced
oscillation between roughly 100 % and 2000 % CPU usage throughout the phase
rather than sustained ~2000 % saturation.

The cause was a fan-out/fan-in shape inside `World::index_with_grid`:

```text
for page in cjindex_pages {                  // (A) serial: load page index
    let chunks = page.par_chunks(N)          // (B) parallel: rayon workers
        .map(process_chunk)
        .collect::<Vec<_>>();
    for chunk in chunks {
        for fic in chunk?.into_iter().flatten() {
            self.integrate_feature_in_cells(fic);  // (C) serial: needs &mut self
        }
    }
}
```

Each 65 536-feature page hit a hard barrier sequence A → B → C → A → B → C.
The 100 % troughs in the CPU graph correspond to phases A (page load) and C
(integration into the dense grid). C is non-trivial: for every retained
feature, the integrator pushes into `World::features` and mutates
`grid.cell_mut(...)` for each cell the feature touches. While C runs, all
rayon workers idle.

`extent_from_cjindex_features` had the same per-page barrier shape, although
its serial tail (`ExtentStats::merge`) is much cheaper than grid integration
and only contributed a brief trough.

A first attempt removed the per-page barrier by collecting all pages into a
single flat `Vec` of chunk slices and running one big `par_iter()`. That
flattened the per-page A↔B oscillation but introduced a new visible 100 %
phase at the start of each phase: the entire `cjindex` page index had to be
materialised serially before any rayon worker could begin.

## Decision

Tyler will run `index_with_grid` and `extent_from_cjindex_features` as
streaming pipelines under `std::thread::scope`, so that page loading,
parallel chunk processing, and the serial tail run concurrently.

For `index_with_grid` this is a 3-stage pipeline:

1. **Page loader thread** streams `cjindex` pages from
   `iter_all_feature_ref_pages(CJINDEX_PAGE_SIZE)`, splits each page into
   `CJINDEX_PARALLEL_CHUNK_SIZE` slices, clones each slice into an owned
   `Vec<IndexedFeatureRef>`, and pushes it to a bounded `sync_channel`.
2. **Worker stage** runs `chunk_rx.into_iter().par_bridge()` on a second
   scope thread. Each chunk parses its features, counts vertices in cells,
   and pushes its `Vec<Option<FeatureInGridCells>>` to a second
   `sync_channel`.
3. **Integrator (main thread)** drains the result channel and applies
   mutations to `features` and `grid`.

`extent_from_cjindex_features` uses the same page-loader stage but collapses
the worker and reduce stages: the main thread runs
`chunk_rx.into_iter().par_bridge().map(process).try_reduce(ExtentStats::default, merge)`
because `ExtentStats::merge` is associative and cheap enough that no
dedicated integrator thread is needed.

To make this representable to the borrow checker without `Mutex`-wrapping
the grid, two supporting refactors were required:

- A new `GridLayout` snapshot struct in `crate::spatial_structs` exposes the
  read-only spatial queries used by workers (`locate_point`,
  `intersect_bbox_ranges`, `intersect_bbox`). It is a `Copy` projection of
  `SquareGrid`'s metadata, so workers can hold it by value while the
  integrator keeps `&mut SquareGrid` for cell mutation.
- `World::index_feature_model` and `World::index_cjindex_feature_refs_chunk`
  no longer take `&self`. They take `&InputSource`, `&GridLayout`, and
  `Option<&Vec<CityObjectType>>` directly, so `index_with_grid` can split
  the `&mut self` borrow into disjoint field borrows.

The bounded `sync_channel(64)` between stages provides back-pressure: if
the integrator falls behind, the worker stage blocks, which in turn blocks
the page loader. Memory use is therefore bounded by the channel capacity
rather than by the number of pages on disk.

The redundant `if !grid_cell.feature_ids.contains(&fid)` guard inside
`integrate_feature_in_cells` was also removed. `fid` is freshly allocated
per feature and `feature_in_cells.cells` already contains unique cellids
(produced by `compress_vertex_counts`), so the guard could never trip.

## Consequences

Good:

- workers can begin processing as soon as the first chunk is loaded; there
  is no longer a serial pre-collection phase
- grid integration now overlaps with worker chunk processing, removing the
  per-page C trough
- page index I/O runs on a dedicated OS thread, separate from the rayon
  CPU pool, so I/O does not steal a worker slot
- back-pressure through `sync_channel(64)` keeps memory bounded; we do not
  buffer the whole dataset's `IndexedFeatureRef`s in memory
- `GridLayout` is a small, `Copy` value, so passing it to workers is
  allocation-free per call

Trade-offs:

- the page loader now clones each chunk into an owned
  `Vec<IndexedFeatureRef>` so it can cross the channel; with
  `CJINDEX_PARALLEL_CHUNK_SIZE = 2048` this is a small per-chunk allocation
  and string clone, but is non-zero
- `par_bridge` adds a mutex-guarded `next()` call per chunk; with chunk
  sizes in the thousands of features the per-chunk overhead is negligible
  relative to feature parsing
- the integrator path is now a free function `integrate_feature_in_cells`
  taking `&mut FeatureSet` and `&mut SquareGrid` rather than a method on
  `World`; this is a small API shape change
- error handling now flows through `std::io::Error` boxes from the
  page-loader and worker channels, rather than propagating native
  `cityjson_lib::Error` and `cityjson_index` errors directly

Neutral:

- `extent_from_cjindex_features` no longer accepts a `&CityIndex`; it opens
  its own index via `InputSource::open_index` so the index can be moved
  into the page-loader thread without the `Sync` requirement
  (`rusqlite::Connection` is `Send` but not `Sync`)

## Rejected Alternatives

- Wrap the cell data in atomics or a sharded mutex so the integrator can
  also be parallelised. The dense grid is up to 1024 × 1024 ≈ 1 M cells,
  each holding `nr_vertices: usize` and `feature_ids: Vec<usize>`. Atomics
  could serve `nr_vertices`, but `feature_ids` would still need a mutex
  per shard. The pipeline change here is structurally simpler and already
  removes the dominant C trough; sharding can be revisited if profiling
  later shows the integrator as the next bottleneck.

- Pre-collect all `IndexedFeatureRef` pages into one flat `Vec` and run
  `par_iter()` on it. This was the first iteration of the change. It
  removes the per-page A↔B↔C oscillation but introduces a long serial
  page-collection prefix at the start of each phase, which was visible as
  a ~100 % CPU phase before the parallel work started.

- Use `rayon::scope` instead of `std::thread::scope` for the spawned
  stages. `rayon::scope` blocks the calling task until all child tasks
  complete, which would defeat the goal of overlapping page loading with
  chunk processing. `std::thread::scope` runs the stages on dedicated OS
  threads while still allowing borrows from the enclosing function.

- Use `rayon::spawn` to fire-and-forget chunk tasks. `rayon::spawn`
  requires `'static` closures, which would force `Arc` wrapping of
  `InputSource`, `GridLayout`, and the type filter. The
  `std::thread::scope` formulation lets these be ordinary borrows.

## Validation Plan

Functional validation:

- run `cargo fmt --all --check`
- run `cargo check --workspace --all-targets --all-features`
- run `cargo test --workspace --all-targets --all-features`
- keep existing parser tests for:
  - unique-assignment cell selection
  - bbox-spanning vs. vertex-bearing cell counts
  - `count_vertices_in_grid` parity with the reference `BTreeMap`
    implementation

Performance validation:

- rerun the same dataset that produced the 100 %/2000 % oscillation
  (Amsterdam-scale `cjindex` with `--cityobject-types Building`)
- confirm with `cpustat` that the extent and vertex-counting phases
  sustain ~2000 % CPU usage across the rayon pool, with only a brief
  startup blip while the first chunk is being loaded
- confirm wall-time is at least as good as the per-page parallel
  baseline from `36b2071`

Expected profiling change:

- the per-page CPU trough during integration disappears
- the prefix CPU trough caused by serial page pre-collection disappears
- there is no measurable change in CityJSON decode time or PROJ behavior
  (those are addressed by ADRs 005, 008, and 009)

## Notes

The integrator stage still runs on a single thread, by design: the dense
grid mutations and the `features` `Vec` push are serialised. The win here
is overlapping that serial work with chunk processing, not eliminating it.
If a future profile shows the integrator as the next dominant cost, the
cell-storage refactor described under "Rejected Alternatives" would be the
next step.
