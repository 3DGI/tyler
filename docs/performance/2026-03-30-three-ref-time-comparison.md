# Three-ref Tyler time comparison

Date: 2026-03-30

Dataset:
`/home/balazs/Data/3DBAG_3dtiles_test/input`

Refs:

- `master` -> `9085570d7e9ca4e15c1cdcc7b46ac7d01524d0f1`
- `623ef236ceec4d2b210c4aa7fcc7cdf1c670d1f1`
- `6219a556d522b016dd82f802815e558a28a0c190`
- patched working tree on `8560b6d` with local changes to:
  - use a new `cjindex::query_iter_with_ids()` path instead of JSON roundtrips for feature ids
  - count only the selected vertices per feature instead of scanning all vertices once per selected cityobject
- current head on `HEAD` after switching Tyler to `cjindex::iter_all*()` for full corpus scans

Method:

- Each ref was built in release mode with `cargo build --release --locked`.
- The timed binary was the ref-local `target/release/tyler`.
- Timing used `/usr/bin/time -v`.
- Dataset cold start was kept, including cjindex sidecar creation/rebuild when required.
- `--3dtiles-tileset-only` was used to keep the workload focused on Tyler.
- A minimal stub `geof` executable was passed because older refs still validate `--exe-geof` even in tileset-only mode.
- `master` and `623e` use the legacy CLI shape:
  `tyler -m <metadata> -f <features> -o <output> --3dtiles-tileset-only --3dtiles-metadata-class benchmark`
- `6219` uses the single-input CLI shape:
  `tyler <input-root> -o <output> --3dtiles-tileset-only --3dtiles-metadata-class benchmark`
- `master` and `623e` were built against the legacy `cjlib` API commit `1e5c814`; `6219` was built against the current sibling `cjlib` checkout. This was necessary to make the historical refs build in the current workspace layout.

Results:

| Ref | Elapsed | User | System | CPU | Max RSS | FS inputs | FS outputs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `master` | `0:05.50` | `8.95s` | `1.09s` | `182%` | `49,188 KB` | `0` | `16` |
| `623ef236ceec4d2b210c4aa7fcc7cdf1c670d1f1` | `0:30.40` | `37.94s` | `1.20s` | `128%` | `63,536 KB` | `0` | `360` |
| `6219a556d522b016dd82f802815e558a28a0c190` | `4:06.44` | `195.78s` | `50.34s` | `99%` | `188,560 KB` | `528` | `34,551,384` |
| patched working tree on `8560b6d` | `3:55.08` | `186.87s` | `47.59s` | `99%` | `188,284 KB` | `15,552` | `34,551,312` |
| current head after `iter_all*()` | `0:47.49` | `44.11s` | `3.33s` | `99%` | `188,776 KB` | `8` | `156,960` |

Relative elapsed time:

- `623e` is `5.53x` slower than `master`.
- `6219` is `44.81x` slower than `master`.
- `6219` is `8.11x` slower than `623e`.
- patched working tree on `8560b6d` is `42.74x` slower than `master`.
- patched working tree on `8560b6d` is `7.73x` slower than `623e`.
- patched working tree on `8560b6d` is `1.05x` faster than `6219` (`11.36s`, about `4.6%` lower elapsed time).
- current head after `iter_all*()` is `8.64x` slower than `master`.
- current head after `iter_all*()` is `1.61x` faster than `623e`.
- current head after `iter_all*()` is `5.18x` faster than `6219` (`218.95s`, about `88.0%` lower elapsed time).

Notes:

- The `6219` run includes the cjindex cold-start cost, which is part of the current architecture and therefore intentionally included.
- The `6219` run is also much more I/O-heavy than the older refs, based on the `time -v` filesystem counters.
- The patched run kept the same cold-start condition by temporarily moving aside the existing `/home/balazs/Data/3DBAG_3dtiles_test/input/.cjindex.sqlite` sidecar before timing and then restoring it afterward.
- The patch removes two real Tyler-side inefficiencies, but the overall improvement is small. The remaining dominant cost is still the broader cjindex path architecture, especially the extra full dataset pass after extent computation and the cold-start reindex itself.
- Switching Tyler to `cjindex::iter_all*()` removes the spatial lookup bottleneck from the full scan. The new cold run spends about `19s` in extent construction, about `19s` in grid indexing, and about `9s` rebuilding the sidecar.

Instrumented cold-start phase breakdown on the patched working tree:

- `cjindex_reindex`: `9.332s`
- `world_from_cjindex`: `111.019s`
- `world_index_with_grid`: `114.873s`
- total elapsed for the instrumented run: `4:06.73`

Phase share of total elapsed (`246.73s`):

- `cjindex_reindex`: about `3.8%`
- `world_from_cjindex`: about `45.0%`
- `world_index_with_grid`: about `46.6%`
- everything else combined: about `4.7%`

This phase split shows that the dominant regression is not the cold-start reindex by itself. The two full cjindex-backed feature reconstruction passes inside Tyler account for roughly `225.9s` of the `246.7s` run, while reindexing accounts for only `9.3s`.
