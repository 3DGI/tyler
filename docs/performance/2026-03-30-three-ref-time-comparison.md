# Three-ref Tyler time comparison

Date: 2026-03-30

Dataset:
`/home/balazs/Data/3DBAG_3dtiles_test/input`

Refs:

- `master` -> `9085570d7e9ca4e15c1cdcc7b46ac7d01524d0f1`
- `623ef236ceec4d2b210c4aa7fcc7cdf1c670d1f1`
- `6219a556d522b016dd82f802815e558a28a0c190`

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

Relative elapsed time:

- `623e` is `5.53x` slower than `master`.
- `6219` is `44.81x` slower than `master`.
- `6219` is `8.11x` slower than `623e`.

Notes:

- The `6219` run includes the cjindex cold-start cost, which is part of the current architecture and therefore intentionally included.
- The `6219` run is also much more I/O-heavy than the older refs, based on the `time -v` filesystem counters.
