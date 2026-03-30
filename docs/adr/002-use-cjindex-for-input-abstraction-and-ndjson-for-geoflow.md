# Use cjindex for Input Abstraction and NDJSON for Geoflow Exchange

## Status

Accepted

## Related Commits

- `93c0f0e` Integrate cjindex inputs with NDJSON tile export

## Context

Tyler originally assumed one narrow input model:

- `--metadata` points to a single base `.city.json` document
- `--features` points to a directory tree of standalone `.city.jsonl` feature files
- tile export passes those original feature file paths to `geoflow`

That assumption became too restrictive once `cjindex` was introduced as the
dataset abstraction layer. `cjindex` supports three storage layouts:

- `feature-files`
- `ndjson`
- `cityjson`

The conversion boundary between Tyler and `geoflow` also needed to stay
simple. Although `geoflow` can read CityJSON sequence data directly, Tyler
should not grow multiple export contracts for downstream conversion. The
requirement for this integration is:

- Tyler may accept multiple storage layouts as input
- Tyler exports only NDJSON / `CityJSONSeq` to `geoflow`

This means the old "pass through original source file paths" design is no
longer sufficient. Tyler needs an internal feature handle that is independent
of how the input dataset is physically stored.

## Decision

Tyler will use `cjindex` as the input and storage-layout abstraction, and it
will use NDJSON as the only exchange format passed to `geoflow`.

Concretely:

1. Tyler's CLI will take one dataset-root input path.
2. That input path may refer either to:
   - a legacy dataset root that contains `metadata.city.json` and standalone
     `.city.jsonl` feature files under that root
   - a `cjindex` dataset root in `feature-files`, `ndjson`, or `cityjson`
     layout
3. Tyler will detect `cjindex` datasets with `cjindex::resolve_dataset()`.
4. For `cjindex` inputs, Tyler will ensure a usable sidecar index exists and
   will derive one shared base metadata document for internal processing.
5. Tyler's world/grid/quadtree state will store backend-agnostic feature
   references instead of only relative feature file paths:
   - legacy relative feature path
   - `cjindex` feature id
6. For each tile, Tyler will resolve the selected feature references and write
   a tile-local NDJSON file.
7. Tyler will continue to invoke `geoflow` through
   `--path_features_input_file`, but that file will now always point at
   NDJSON created by Tyler for the tile.

## Consequences

Good:

- Tyler can ingest all three storage backends supported by `cjindex`
- the `tyler -> geoflow` boundary is stable and format-specific
- tile export no longer depends on original source file layout
- CityJSON input can be normalized into NDJSON without changing downstream
  tooling
- legacy feature-file input keeps working

Trade-offs:

- Tyler now owns tile-local NDJSON materialization
- `cjindex` datasets must expose a consistent shared metadata document for
  Tyler's current processing model
- export now includes an extra write step even when the source is already
  NDJSON
- internal indexing code is more abstract because it must support both path-
  based and id-based feature references

## Rejected Alternatives

- Pass original source files directly to `geoflow` for every layout.
  This keeps Tyler tightly coupled to storage layout details and does not give
  one stable downstream contract.

- Support multiple Tyler export modes for `geoflow`, such as both CityJSON and
  NDJSON.
  This widens the boundary without a clear need and makes downstream behavior
  harder to reason about and test.

- Restrict `cjindex` integration to NDJSON only.
  This would leave `cityjson` and `feature-files` datasets outside the new
  abstraction, defeating the purpose of integrating `cjindex`.

- Keep storing only feature file paths in Tyler internals.
  That works only for the legacy layout and breaks as soon as features no
  longer correspond to standalone files on disk.

## Notes

This ADR is about the Tyler input/export contract, not about changing
`geoflow`. `geoflow` still receives paths through
`--path_features_input_file`; the decision is that those paths now always
resolve to Tyler-generated NDJSON for the current tile.
