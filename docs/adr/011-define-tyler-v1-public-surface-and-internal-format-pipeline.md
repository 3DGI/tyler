# Define Tyler v1.0 Public Surface and Internal Format Pipeline

## Status

Proposed

## Related Issues

- Multi-format Tyler surface: #110, #111, #112, #113, #119
- Format-agnostic orchestration: #109
- `cityjson-convert` expansion: #103, #106, #107, #108
- Shared converter internals: #104, #105
- GPKG typing and nullable attributes: #52
- Mesh coloring by semantic surface: #49
- Runtime and operations: #5, #6, #8, #37
- Dependency and security boundary: #116, #120

## Context

Tyler started as a 3D Tiles exporter for CityJSON-compatible inputs. Through
v0.3.14, Tyler prepared tiled CityJSONSeq-style exchange data and delegated
format conversion to the external C++ Geoflow software through subprocesses. In
v0.4.1, Tyler removed that external runtime dependency and introduced the
in-repository `cityjson-convert` implementation for native Rust GLB conversion.
That release is also the optimized baseline for Tyler's GLB / 3D Tiles path:
the v1.0.0 work broadens the output surface, but it must maintain or improve
the performance characteristics of v0.4.1 for equivalent GLB / 3D Tiles runs.

That history leaves two boundaries to make explicit for v1.0.0. Tyler owns the
tiling pipeline and the preparation of per-tile `CityModel` values.
`cityjson-convert` owns pure format conversion and serialization of those
prepared models. The v1.0.0 milestone expands Tyler into a multi-format tiler
and extends `cityjson-convert` into the shared format conversion and
serialization layer for the repository.

The milestone issues describe several connected changes:

- Tyler needs a public CLI surface that can produce 3D Tiles, CityJSON,
  CityJSONSeq, OBJ, and GPKG outputs.
- Some runs should emit multiple output formats from one input scan and one tile
  plan.
- Format-specific options must be validated as configuration, before input
  indexing or tile materialization begins.
- Tyler should prepare the per-tile `CityModel` boundary, including feature
  selection by tile membership, filtering, attribute selection, and other
  tiling-driven model shaping before conversion begins. Mesh-level geometric
  clipping at tile boundaries is not part of this stage; it remains a
  post-triangulation step inside `cityjson-convert`.
- `cityjson-convert` should own pure format conversion and serialization from a
  prepared `CityModel` into the requested output format, including the
  mesh-level seam clipping required for seamless 3D Tiles content.
- Operational controls such as file logging and explicit job control are part
  of the v1.0 surface, not profiling-only developer switches.

Without an explicit target architecture, the implementation risks growing one
format at a time around the current single-output GLB path. That would make
multi-format output slower than necessary and would blur responsibility between
Tyler's tiling pipeline and `cityjson-convert`'s format conversion pipeline.

## Decision

Tyler v1.0 will expose a format-agnostic tiling CLI backed by an internal
pipeline that separates orchestration, tiling, tile `CityModel` materialization,
and output conversion.

The v1.0 pipeline must preserve the optimized v0.4.1 GLB / 3D Tiles path as a
performance baseline. Adding output backends, backend validation, and
multi-format dispatch must not make equivalent single-format `3dtiles` runs
slower without a measured and accepted trade-off. Where the new architecture
touches shared input scanning, grid/quadtree construction, tile materialization,
or GLB conversion dispatch, the intended direction is equal or better
throughput and memory behavior than v0.4.1.

### Tyler v1.0 Public Surface

Tyler is a command-line tiler for CityJSON-compatible datasets resolved through
`cityjson-index`. A run selects an input dataset, constructs a tiling plan, and
writes one or more tiled outputs.

The `--format` option selects output backends. Tyler v1.0 supports:

- `3dtiles`
- `cityjson`
- `cityjsonseq`
- `obj`
- `gpkg`
- `tsv`

`--format` may be passed more than once in a single run when the selected
formats share the same tiling scheme. Tyler builds one extent/grid/quadtree plan
for the run and reuses one tile `CityModel` materialization path where possible.
Each selected backend receives the same tile model boundary and writes its own
format-specific output.

Common CLI options cover:

- input and output locations
- CityObject type filtering
- LoD selection
- attribute selection and attribute preservation
- grid size, quadtree, and extent controls
- logging, including logging to a file
- explicit job control and parallelism

Format-specific options are only valid for compatible formats.
Format-specific options are prefixed with the selected format name, following the existing `--3dtiles-*` convention.
Candidate prefixes for v1.0 include: `--cityjson-*`, `--cityjsonseq-*`, `--obj-*`, and `--gpkg-*`.
Tyler validates them before input processing starts. Invalid combinations, missing required
format options, ignored options, and conflicting options fail during
configuration validation rather than after index resolution or tile work.

The public behavior is format-oriented rather than implementation-oriented:
users select datasets, filters, tiling controls, operational controls, and output
formats. Tyler decides which intermediate work can be shared by the selected
formats.

### Examples

Short examples of the intended v1.0 CLI shape. Existing option names are kept
where v0.4.1 already has them; repeated `--format`, `--format obj`,
`--format gpkg`, `--jobs`, and `--log-file` are part of the future v1.0
surface described by this ADR.

Common CLI parameters have short versions for convenience:

- `-v` for `--version`,
- `-h` for `--help`,
- `-i` for `--input`,
- `-o` for `--output`,
- `-f` for `--format`,
- `-j` for `--jobs`,
- `-l` for `--log-file`

```sh
tyler data/amsterdam.city.jsonl \
  --output out/amsterdam-tiles \
  --format 3dtiles \
  --3dtiles-implicit \
  --3dtiles-content-clip-to-tile-bounds \
  --3dtiles-content-add-bv \
  --3dtiles-metadata-class cityobject \
  --object-type Building \
  --object-type BuildingPart \
  --lod-building-part 2.2 \
  --object-attributes measuredHeight:float,function:string \
  --include-parent-attributes \
  --smooth-normals \
  --jobs 16 \
  --log-file out/amsterdam-tiles/tyler.log
```

```sh
tyler data/amsterdam.city.jsonl \
  --output out/amsterdam-review \
  --format cityjson \
  --format cityjsonseq \
  --object-type Building \
  --object-type Road \
  --lod-building 2.0 \
  --lod-road 1.0 \
  --grid-cellsize 250 \
  --grid-minz -10 \
  --grid-maxz 120
```

```sh
tyler data/tile-index.cjindex \
  --output out/shared-mesh-and-review \
  --format 3dtiles \
  --format obj \
  --object-type Building \
  --lod-building 2 \
  --color-building '#d8b365' \
  --3dtiles-geometric-error-factor 0.018
```

```sh
tyler data/region.city.json \
  --output out/region-export \
  --format gpkg \
  --object-type Building \
  --object-type Road \
  --object-attributes class:string,height:float,owner:string
```

```sh
tyler data/debug.city.jsonl \
  --output out/debug-cityjson \
  --format cityjson \
  --format gpkg \
  --object-type Building \
  --object-attributes measuredHeight:float,status:string \
  --debug-dump-data
```

This configuration is rejected during validation because a 3D Tiles-only option
is used without selecting a compatible output backend:

```sh
tyler data/region.city.json \
  --output out/invalid-cityjson \
  --format cityjson \
  --3dtiles-content-clip-to-tile-bounds
```

#### 3DBAG

The current commands for generating the 3DBAG data.

3D Tiles:

```sh
RUST_LOG=error \
RAYON_NUM_THREADS=20 \
tyler \
--output /data/3dtiles_test_geof \
--metadata /data/metadata.json \
--features /data/bouwlagen_features \
--object-type Building \
--object-type BuildingPart \
--lod-building-part 2.2 \
--lod-building 2.2 \
--qtree-capacity=280000 \
--grid-minz="-50" \
--grid-maxz="400" \
--3dtiles-metadata-class building \
--object-attribute=b3_bouwlagen:int --object-attribute=b3_dak_type:string --object-attribute=b3_extrusie:string --object-attribute=b3_h_maaiveld:float --object-attribute=b3_h_nok:float --object-attribute=b3_is_glas_dak:bool --object-attribute=b3_kas_warenhuis:bool --object-attribute=b3_mutatie_ahn3_ahn4:bool --object-attribute=b3_mutatie_ahn4_ahn5:bool --object-attribute=b3_n_nok:int --object-attribute=b3_n_vlakken:int --object-attribute=b3_nodata_fractie_ahn3:float --object-attribute=b3_nodata_fractie_ahn4:float --object-attribute=b3_nodata_fractie_ahn5:float --object-attribute=b3_nodata_radius_ahn3:float --object-attribute=b3_nodata_radius_ahn4:float --object-attribute=b3_nodata_radius_ahn5:float --object-attribute=b3_opp_buitenmuur:float --object-attribute=b3_opp_dak_plat:float --object-attribute=b3_opp_dak_schuin:float --object-attribute=b3_opp_grond:float --object-attribute=b3_opp_scheidingsmuur:float --object-attribute=b3_puntdichtheid_ahn3:float --object-attribute=b3_puntdichtheid_ahn4:float --object-attribute=b3_puntdichtheid_ahn5:float --object-attribute=b3_pw_bron:string --object-attribute=b3_pw_datum:string --object-attribute=b3_pw_onvoldoende:bool --object-attribute=b3_pw_selectie_reden:string --object-attribute=b3_rmse_lod12:float --object-attribute=b3_rmse_lod13:float --object-attribute=b3_rmse_lod22:float --object-attribute=b3_t_run:int --object-attribute=b3_val3dity_lod12:string --object-attribute=b3_val3dity_lod13:string --object-attribute=b3_val3dity_lod22:string --object-attribute=b3_volume_lod12:float --object-attribute=b3_volume_lod13:float --object-attribute=b3_volume_lod22:float --object-attribute=identificatie:string --object-attribute=oorspronkelijkbouwjaar:int --object-attribute=status:string \
> tyler_geof.log 2>&1
```

The Tyler v1.0 commands for generating the 3DBAG data.

3DTiles and OBJ formats are created per LoD 1.2, 1.3, and 2.2.
The OBJ format does not support attributes, so the `--object-attributes` option is ignored.
We only include the `BuildingPart` type in both formats, and use the `--include-parent-attributes` option to include the attributes of the `Building` type, in case of the 3DTiles format.

```sh
tyler \
--jobs 20 \
--log-file tyler_geof.log \
--log-level error \
--format 3dtiles \
--format obj \
--output /data/3dtiles_test_geof \
--input /data/bouwlagen_features \
--object-type BuildingPart \
--lod-building-part 2.2 \
--qtree-capacity=280000 \
--grid-minz="-50" \
--grid-maxz="400" \
--include-parent-attributes \
--3dtiles-metadata-class building \
--object-attributes=b3_bouwlagen:int,b3_dak_type:string,b3_extrusie:string,b3_h_maaiveld:float,b3_h_nok:float,b3_is_glas_dak:bool,b3_kas_warenhuis:bool,b3_mutatie_ahn3_ahn4:bool,b3_mutatie_ahn4_ahn5:bool,b3_n_nok:int,b3_n_vlakken:int,b3_nodata_fractie_ahn3:float,b3_nodata_fractie_ahn4:float,b3_nodata_fractie_ahn5:float,b3_nodata_radius_ahn3:float,b3_nodata_radius_ahn4:float,b3_nodata_radius_ahn5:float,b3_opp_buitenmuur:float,b3_opp_dak_plat:float,b3_opp_dak_schuin:float,b3_opp_grond:float,b3_opp_scheidingsmuur:float,b3_puntdichtheid_ahn3:float,b3_puntdichtheid_ahn4:float,b3_puntdichtheid_ahn5:float,b3_pw_bron:string,b3_pw_datum:string,b3_pw_onvoldoende:bool,b3_pw_selectie_reden:string,b3_rmse_lod12:float,b3_rmse_lod13:float,b3_rmse_lod22:float,b3_t_run:int,b3_val3dity_lod12:string,b3_val3dity_lod13:string,b3_val3dity_lod22:string,b3_volume_lod12:float,b3_volume_lod13:float,b3_volume_lod22:float,identificatie:string,oorspronkelijkbouwjaar:int,status:string
```

CityJSON and GPKG and TSV formats are created with all LoD-s within one file.

```sh
tyler \
--jobs 20 \
--log-file tyler_geof.log \
--log-level error \
--format cityjson \
--format gpkg \
--format tsv \
--output /data/3dtiles_test_geof \
--input /data/bouwlagen_features \
--qtree-capacity=280000 \
--gpkg-split-semantics \
--gpkg-split-lod \
--tsv-omit-null-rows \
--object-attributes=b3_bouwlagen:int,b3_dak_type:string,b3_extrusie:string,b3_h_maaiveld:float,b3_h_nok:float,b3_is_glas_dak:bool,b3_kas_warenhuis:bool,b3_mutatie_ahn3_ahn4:bool,b3_mutatie_ahn4_ahn5:bool,b3_n_nok:int,b3_n_vlakken:int,b3_nodata_fractie_ahn3:float,b3_nodata_fractie_ahn4:float,b3_nodata_fractie_ahn5:float,b3_nodata_radius_ahn3:float,b3_nodata_radius_ahn4:float,b3_nodata_radius_ahn5:float,b3_opp_buitenmuur:float,b3_opp_dak_plat:float,b3_opp_dak_schuin:float,b3_opp_grond:float,b3_opp_scheidingsmuur:float,b3_puntdichtheid_ahn3:float,b3_puntdichtheid_ahn4:float,b3_puntdichtheid_ahn5:float,b3_pw_bron:string,b3_pw_datum:string,b3_pw_onvoldoende:bool,b3_pw_selectie_reden:string,b3_rmse_lod12:float,b3_rmse_lod13:float,b3_rmse_lod22:float,b3_t_run:int,b3_val3dity_lod12:string,b3_val3dity_lod13:string,b3_val3dity_lod22:string,b3_volume_lod12:float,b3_volume_lod13:float,b3_volume_lod22:float,identificatie:string,oorspronkelijkbouwjaar:int,status:string
```

### Mapping of CityModel to Output Formats

#### Comma Separated Values (CSV) or Tab Separated Values (TSV)

Geometry-less output uses the shared CityJSON tabular projection schema defined
in [ADR 013](013-define-shared-cityjson-tabular-projection-schema.md) for
CityModel, CityObject, and Semantic attributes in the input.
The tabular output enables easy processing by spreadsheets and other standard data analysis tools.

The geometry-less tabular output of the attributes could be powerful in combination with a generic parquet format for enabling efficient analytics on large datasets.
Could be the last step in a piped data processing chain of citymodel reshaping and subsetting in the command line, where the result piped into a parquet file that can be queried with generic data
analysis tools.

CityArrow and CityParquet already provide the reference projection direction;
Tyler's CSV, TSV, and GeoPackage outputs should align with that logical schema.
Physical formats may encode that logical projection differently.

#### GeoPackage

The GeoPackage mapping can be generalized to other GIS formats that store geometries as Simple Features in a tabular format.

GeoPackage specs: https://www.geopackage.org/spec/

The detailed CityJSON-to-GeoPackage schema mapping is defined in
[ADR 012](012-define-cityjson-to-geopackage-schema-mapping.md), using the shared
tabular projection from
[ADR 013](013-define-shared-cityjson-tabular-projection-schema.md). In summary,
GeoPackage output is GIS-first: coordinates are written as real-world XYZ
coordinates in one declared CRS, feature layers are split into homogeneous
CityObject type and geometry-family layers, CityJSON geometry templates are
resolved, solids are exported as boundary `MultiPolygonZ` features, and
CityObject hierarchy is preserved in a separate relation table.

By default, appearance and metadata are dropped except for CRS metadata.
Semantic primitives are dropped unless `--gpkg-split-semantics` is set, in
which case they are exported into separate feature layers with a semantic
relation table. CityObject, semantic, address, and supported extension
attributes are encoded from the shared logical tabular projection, with JSON
text used only for fields projected as explicit `Json` or for physical cells a
format cannot otherwise represent.

### CLI parameter matrix

CLI parameter matrix of functionality per output format.

### Internal Architecture

`main` will be refactored into orchestration. It will coordinate these phases:

1. Parse CLI arguments into raw configuration.
2. Validate the configuration.
3. Resolve input through `cityjson-index`.
4. Build the extent, grid, and quadtree plan.
5. Materialize per-tile `CityModel` values.
6. Dispatch those tile models to selected output backends.
7. Report progress, logging, and operational errors consistently.

Selected formats will be represented as a list of output backends instead of a
single enum value. Each backend declares its compatibility constraints, required
options, unsupported options, tiling-scheme requirements, and writer entrypoint.
The orchestration layer uses those declarations to decide whether a backend set
can share one run.

Configuration validation is an explicit phase between parsing and input
processing. It checks:

- selected format compatibility
- required options for each format
- options that are ignored by the selected formats
- mutually exclusive options
- output path shape and overwrite behavior
- operational controls such as job count and log-file configuration

Tyler remains responsible for tiling. Its ownership includes:

- interpreting CLI tiling and filtering intent
- resolving and indexing input datasets
- selecting features
- calculating extents
- constructing grids and quadtrees
- assigning features to tiles
- assembling per-tile `CityModel` values
- preparing those tile models for output, including feature-level tile
  membership, LoD filtering, object-type filtering, attribute selection, and
  semantic, material, or color annotations needed by the output configuration
- coordinating shared work across output backends

Tyler does not perform per-vertex geometric clipping against tile bounds. It
selects whole features whose extent overlaps the tile and hands the resulting
`CityModel` to `cityjson-convert` in source CRS.

`cityjson-convert` is responsible for pure conversion and serialization of a
prepared `CityModel` into output formats. Its ownership includes:

- GLB and 3D Tiles content conversion
- CityJSON output conversion
- CityJSONSeq output conversion
- OBJ output conversion
- GPKG output conversion
- format-specific encoding of already-selected attributes, including GPKG typing
  and nullable attributes
- format-local mesh construction required by a target format, such as
  triangulating CityJSON surfaces for mesh-based serialization
- post-triangulation mesh-level clipping against an optional clip volume
  (source-CRS bounding box or geographic region), used by the 3D Tiles path to
  produce seamless tile content. The geographic-region mode is CRS-aware and
  uses bisection against EPSG:4979 planes; it is only invoked when the 3D Tiles
  backend requests it. Other writers never receive a clip volume and never pull
  the geographic clipper into their path.
- serializing color or material information that is already present in the
  prepared model

Shared format-conversion primitives will be moved out of format-specific writers
into reusable `cityjson-convert` modules when they are pure serialization
concerns. In particular, mesh construction and attribute encoding should not be
private implementation details of the GLB writer when OBJ, GPKG, or future
writers need the same format-level behavior. Tiling-driven preparation that
operates on whole features — selection, filtering, attribute shaping — stays in
Tyler before the `CityModel` reaches `cityjson-convert`. Mesh-level operations
such as triangulation and tile-seam clipping stay in `cityjson-convert` because
they require the triangulated mesh and, for 3D Tiles geographic regions, a
CRS-aware intersection. Clipping is therefore expressed as a post-triangulation
step parameterized by an optional clip volume, not as a Tyler-level
`CityModel` transformation.

The PROJ dependency will be moved to `cityjson-lib` and implemenent `cityjson_lib::ops::reproject`.
The `cityjson_lib::ops::reproject` applies the coordinate transformation on the `CityModel.vertices`.

Logging to a file and explicit job control are v1.0 operational controls. They
are configured before the pipeline starts and are available across orchestration,
indexing, tiling, materialization, and conversion phases.

## Consequences

Good:

- Tyler gets a stable v1.0 CLI contract for the milestone's output formats.
- Multi-format runs can share input resolution, tiling, and tile model
  materialization when formats are compatible.
- Single-format GLB / 3D Tiles runs keep v0.4.1 as the performance baseline
  while the output surface expands.
- Invalid format-specific configuration fails before expensive work begins.
- The codebase gets a clearer ownership boundary: Tyler prepares tile models;
  `cityjson-convert` converts and serializes prepared models.
- Shared pure conversion logic becomes reusable by multiple writers instead of
  being hidden in the GLB path.
- CRS dependency ownership becomes more explicit and can serve both Tyler and
  library-level callers.

Trade-offs:

- Backend declarations add structure that is heavier than a single output enum.
- Validation must understand both global options and backend-specific option
  compatibility.
- Some existing writer code must move before new formats can share pure
  conversion behavior such as mesh construction and format-specific attribute
  encoding cleanly.
- Tyler must keep the tile-model preparation path explicit enough that
  converters never need to infer tiling intent from CLI configuration.
- Multi-format output makes progress reporting and error reporting more
  important because one run can contain several writer phases.

Neutral:

- The tile `CityModel` remains the boundary between Tyler and
  `cityjson-convert`.
- Individual backends may still perform format-specific optimization after
  receiving a tile model.
- Format conversion may transform the prepared model into target-format data
  structures, but it must not perform tiling-driven model preparation.
- Formats that cannot share a tiling scheme are rejected as one run rather than
  silently split into multiple independent runs.

## Rejected Alternatives

- Keep Tyler centered on a single `--format` enum and add one-off branches for
  each new output. That would be simpler for the next format only, but it would
  make shared multi-format runs awkward and would push validation into scattered
  late pipeline checks.

- Let each output format perform its own feature selection and tiling. That
  would keep backends independent, but it would repeat input scans and make
  cross-format output from one run slower and harder to reason about.

- Move tiling into `cityjson-convert`. The converter should understand how to
  turn `CityModel` data into external formats; it should not own Tyler's CLI
  semantics, tile planning, feature selection, or indexing strategy.

- Move tile model preparation into `cityjson-convert` so every writer can decide
  how to clip, filter, and shape models. That would make format writers depend
  on Tyler's tiling semantics and would duplicate orchestration decisions across
  serialization code.

- Keep mesh construction and format-specific attribute encoding private to the
  GLB writer. That preserves the current local implementation, but it prevents
  other mesh or tabular format writers from using the same pure conversion
  semantics and encourages duplicate serialization code.

- Treat logging-to-file and job control as post-v1.0 operational polish. They
  affect reproducibility and runtime behavior for large datasets, so they belong
  in the public v1.0 surface.

## Validation Plan

Documentation validation:

- confirm this file is numbered after ADR 010 and uses the ADR filename pattern
- confirm the related issue list uses only v1.0.0 milestone issue numbers
- confirm the public-surface section describes the v1.0 target behavior without
  migration wording
- confirm the internal-architecture section separates Tyler orchestration from
  `cityjson-convert` conversion responsibility
- confirm Markdown headings match the existing ADR style

Implementation validation:

- add parser and configuration-validation tests for repeated `--format` values
- add validation tests for incompatible formats and format-specific options
- add integration coverage showing compatible backends can share one tile plan
- add writer-level tests in `cityjson-convert` for CityJSON, CityJSONSeq, OBJ,
  GPKG, and GLB conversion from tile `CityModel` inputs
- keep Tyler tests for tiling-driven model preparation, including feature-level
  tile membership, filtering, LoD selection, and attribute selection
- keep `cityjson-convert` tests for pure format conversion, including mesh
  construction, post-triangulation clipping against source-CRS bounding boxes
  and 3D Tiles geographic regions, format-specific attribute encoding, semantic
  color/material serialization, CRS-sensitive serialization behavior, and GPKG
  typing
- compare representative single-format `--format 3dtiles` runs against v0.4.1
  and confirm wall time, CPU utilization, and peak memory are maintained or
  improved for equivalent inputs and options
- compare multi-format runs against separate single-format runs and confirm
  shared input resolution, tiling, and tile materialization avoid repeated work
- run `just ci-check` for the full non-mutating local validation sequence
- use `just fmt-check`, `just lint`, `just check`, `just build`, and
  `just test` for targeted validation while implementing
- if future validation needs commands not covered by these recipes, add or
  update Justfile recipes first and reference those recipes here
- CityJSON schema-valid fake data is generated with the [cityjson-fake](https://github.com/3DGI/cityjson-rs/tree/main/crates/cityjson-fake) crate. The cityjson-fake tool can generate fake CityJSON
  data in any shape and size; however, geometries are dummy values. It will be most useful for testing that `cityjson-convert` supports the full CityJSON schema.

## Notes

This ADR records the intended v1.0 target state for the grouped milestone
issues. It does not require every issue to land in one patch set. The important
architectural constraint is that new output formats should attach to the shared
tiling and tile-model pipeline instead of creating independent format-specific
pipelines inside Tyler.
