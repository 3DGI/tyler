# cjindex NDJSON Plan

## Goal

Make `tyler` accept the three source layouts covered by `cjindex`, but keep the export contract
simple:

- `tyler` only exports NDJSON / `CityJSONSeq`
- `geoflow` reads those NDJSON files for conversion

## Revised Design

1. Use `cjindex` only as the input/layout abstraction.
2. Stop treating the input path as "directory of standalone feature files" only.
3. Normalize every non-legacy input into NDJSON before the `geoflow` step.
4. Keep the `tyler` -> `geoflow` boundary fixed: NDJSON file paths in
   `--path_features_input_file`.

## Input Modes

### Legacy feature files

Keep the current path working. `tyler` already knows how to index standalone feature files.

### cjindex NDJSON

Use `cjindex` to detect the dataset and reuse NDJSON directly where possible.

### cjindex CityJSON

Use `cjindex` to read the dataset and emit NDJSON that `tyler` can pass to `geoflow`.

## Required Refactor in tyler

`tyler` currently stores feature file paths in the world/grid structures. That is too narrow.

Replace that with a backend-agnostic feature reference:

- legacy feature-file path
- `cjindex` feature id or equivalent handle

That lets tile export resolve selected features into NDJSON regardless of the original storage
layout.

## Export Path

For each tile:

1. Collect the selected feature references from the quadtree/grid.
2. Resolve those references through the active backend.
3. Write NDJSON for that tile.
4. Pass the NDJSON path to `geoflow`.

The key rule is that `tyler` does not need to export feature files anymore for the `geoflow`
boundary. NDJSON is the only transport format.

## Implementation Steps

1. Add `cjindex` as a dependency in `tyler`.
2. Detect whether the single input path is:
   - legacy standalone feature files
   - `cjindex` NDJSON
   - `cjindex` CityJSON
3. Refactor world/index data so it is not hard-wired to feature file paths.
4. Add a `cjindex`-backed extent/indexing path.
5. Add tile-local NDJSON export for `cjindex` inputs.
6. Keep the `geoflow` invocation unchanged except that it now always receives NDJSON paths.
7. Add tests for legacy feature files, `cjindex` NDJSON, and `cjindex` CityJSON.

## Immediate Next Step

Define the `cjindex` API that `tyler` needs for tile-local NDJSON emission, because that is the
new integration boundary.
