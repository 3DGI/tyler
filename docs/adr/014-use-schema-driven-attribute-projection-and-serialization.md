# Use Schema-Driven Attribute Projection and Serialization

## Status

Proposed

## Supersedes

- [ADR 011](011-define-tyler-v1-public-surface-and-internal-format-pipeline.md)
- [ADR 012](012-define-cityjson-to-geopackage-schema-mapping.md)
- [ADR 013](013-define-shared-cityjson-tabular-projection-schema.md)

## Context

Tyler and `cityjson-convert` serialize CityObject attributes into CityJSON,
CityJSONSeq, TSV, GeoPackage, and GLB, including GLBs used as 3D Tiles
content. OBJ has no attribute representation. The existing design grew from
three related but separate proposals:

- ADR 011 assigned model shaping to Tyler and physical output conversion to
  `cityjson-convert`, but used `--object-attributes name:type` for attribute
  selection and coercion.
- ADR 012 defined an input-driven GeoPackage geometry and layer mapping, but
  mixed that mapping with locally inferred attribute columns.
- ADR 013 aligned Tyler's tabular inference with CityArrow, but left the same
  value-shape rules implemented independently in `cityjson-arrow` and
  `cityjson-convert`.

The 3DBAG project has a more expressive attribute specification. Its
`attributes.json` describes an attribute's type, nullability, precision,
format, categorical values, and the output formats and locations where the
attribute belongs. It was originally designed to validate already serialized
3DBAG products, and its current JSON Schema hard-codes 3DBAG formats, layer
names, languages, shallow array items, and unqualified `int` and `float`
types.

The reusable part of this design is not 3DBAG's particular output topology. It
is the separation between an attribute's serialization contract and the
format-specific writer that realizes the contract. Tyler remains an ETL and
data-wrangling tool: the input model is transformed into the desired
CityObject, geometry, hierarchy, and semantic shape before serialization.
Attribute rules may select and transform attributes inside locations produced
from that input, but do not manufacture a 3DBAG-specific layer structure.

The inference problem is also not inherently Tyler-specific or tabular.
CityJSON, CityJSONSeq, GLB metadata, Arrow, Parquet, TSV, and GeoPackage all
need the same recursive value shapes, numeric conflict rules, schema
applicability, and mismatch behavior. A `cityjson-tabular` crate would give
non-tabular writers the wrong dependency and would prematurely introduce a
common table model for writers whose row layouts remain different.

## Decision

### Shared attribute domain

Add a lightweight `cityjson-attributes` crate to the `cityjson-rs` workspace.
The crate owns:

- the generalized attribute-specification model and semantic validation;
- recursive attribute value shapes and whole-field projection;
- numeric widening and heterogeneous `Json` fallback;
- exact `(format, location, attribute)` applicability resolution;
- schema-driven coercion and value validation; and
- structured diagnostics and immutable resolved projections.

The crate does not own layer or table creation, flattened column naming, SQL,
TSV escaping, JSON token writing, Arrow arrays, GLB buffer layout, CLI
handling, or warning presentation. Format implementations translate a resolved
logical projection into their native representation.

`cityjson-arrow` uses the shared inference and keeps its Arrow-specific
`ProjectionLayout`, schema, arrays, and package transport. `cityjson-json`
consumes a resolved projection when writing CityJSON or CityJSONSeq.
`cityjson-convert` uses the same projection in its tabular and GLB writers.
`cityjson-lib` re-exports the public attribute API.

No legacy API wrappers or aliases are retained. The repositories are under
active development, and a clean break gives the new shared domain one
vocabulary.

### Generalized specification

Generalize the 3DBAG attribute specification schema while keeping its
top-level dictionary of exact attribute names.

The type vocabulary is:

```text
int8 int16 int32 int64
uint8 uint16 uint32 uint64
float32 float64
bool string date datetime array
```

Array `items` use the same value specification recursively and are required
for `array`. Nested object schemas are not introduced. Nested input objects
remain supported through default projection and `Json` fallback, but are not
schema-addressable in this version.

Every attribute requires `type`, `nullable`, and a non-empty `appliesTo`.
`nullable` is `true`, `false`, or `null`; `null` means that the specification
does not assert nullability. Format names, location names, and translation
language keys are open strings. Descriptive 3DBAG fields become optional in
the generalized contract, while a separate 3DBAG profile requires its richer
metadata and Dutch and English translations.

`precision` is valid only for floats and is an integer from 0 through 17. It
means a maximum number of decimal places. Values are rounded to nearest and
trailing zeroes are not added.

Dates and datetimes remain strings. Their required `valueFormat` is a
Chrono-compatible strftime expression. Partial formats such as `%Y` are
validated through Chrono's lower-level parsed-field API rather than through a
constructor that requires a complete calendar date. A value that does not
match remains unchanged, retains its declared Date or Datetime field type, and
produces a warning.

The optional `values` catalogue is validated against the declared scalar type.
Unlisted values are preserved and warned about.

The existing 3DBAG specification is migrated in the same breaking change:

- `int` becomes `int32`;
- `float` becomes `float32`;
- Cesium/3D Tiles applicability becomes `glb`;
- date formats become `%Y` or `%Y-%m-%d`; and
- datetime formats become `%Y-%m-%dT%H:%M:%S%.3f`.

The old schema and vocabulary are not supported.

### Input-driven output structure

Writers derive their native locations solely from transformed input and their
established format rules. The specification cannot create, rename, map, split,
or merge locations.

The following input-driven choices from ADRs 011-013 are retained:

- Tyler owns dataset orchestration, selection, hierarchy/LoD transformations,
  parent-attribute inheritance, global preflight, and tiling.
- `cityjson-convert` owns physical conversion from a prepared model.
- GeoPackage layers remain determined by CityObject type, geometry family, and
  optional LoD splitting. Existing fixed identity and geometry constraints
  remain.
- TSV and GeoPackage retain writer-specific row layouts and deterministic
  flattened names.
- CityArrow and CityParquet retain native nested transport layouts.

The actual native location is reported to the shared resolver:

| Output | Format identifier | Location identifier |
|---|---|---|
| CityJSON | `cityjson` | CityObject or semantic type |
| CityJSONSeq | `cityjson` | CityObject or semantic type |
| GeoPackage | `gpkg` | generated table name |
| TSV | `tsv` | stable logical table name, such as `cityobjects` or `semantics` |
| GLB | `glb` | configured structural-metadata class |

CityJSONSeq deliberately reuses `cityjson` because it has the same attribute
representation and locations. Tyler's 3D Tiles output also uses `glb` because
its tile content is serialized by the same GLB writer with the same attribute
representation and locations. The tileset is an orchestration concern outside
attribute serialization.

### Optional schema-driven projection

Without an attribute specification, every writer follows its current
serialization behavior.

When `--attribute-specification <PATH>` is supplied, `appliesTo` is an
attribute whitelist for the exact output format and location:

- a declared attribute at the exact location is included and schema-driven;
- an input attribute absent from the specification or declared only for
  another location is omitted and warned about;
- a schema-expected attribute absent from a produced location is not
  synthesized and produces one aggregate warning; and
- an applicability location not produced by the input does not itself warn.

The schema applies to CityObject and semantic attributes. It does not cause a
writer to materialize semantics that the selected output would not otherwise
contain.

Tyler's `--object-attributes name:type` option is removed. The specification
provides a more precise, format/location-aware replacement for attribute
selection and transformation. `--include-parent-attributes` remains an
upstream ETL operation and runs before schema projection.

### Coercion and whole-field fallback

Schema projection resolves complete fields, not individual cells. Tyler
preflights the selected dataset before parallel output so every tile uses the
same resolution.

The existing permissive `--object-attributes` conversions are retained and
extended to the explicit-width vocabulary:

- booleans and finite numbers can become strings;
- recognized boolean strings and zero/nonzero finite numbers can become
  booleans;
- booleans and numeric strings can become numeric values;
- finite floats can become integers by truncation toward zero;
- integer and float narrowing is accepted only when the result is in range;
- date and datetime values must be strings; and
- arrays apply item conversion recursively.

If any non-null value cannot be converted, the declared attribute remains
selected but the entire field falls back to the existing inferred,
loss-preserving representation. One structured type/coercion diagnostic counts
the affected values.

The common conflict lattice remains:

- unsigned and signed integers merge to signed integer;
- integer and float values merge to float;
- compatible lists and structs merge recursively;
- null affects nullability without creating a type conflict; and
- incompatible stable shapes become `Json`.

### Nulls, missing values, and categories

Once an attribute exists in a produced location, omitted values and explicit
nulls both count as violations of `nullable: false`. An attribute absent from
every record in a produced location receives only the missing-attribute
diagnostic. `nullable: true` and `nullable: null` do not warn.

Categorical validation happens after coercion. Values absent from the optional
catalogue are preserved and counted.

### Physical encodings

Writers translate the shared decision as follows:

- CityJSON and CityJSONSeq emit native JSON values. Precision-controlled
  floats remain JSON numbers, never strings.
- TSV emits scalar text and compact JSON for arrays and `Json` fallback.
- GeoPackage maps integers to `INTEGER`, floats to `FLOAT`, booleans to
  `BOOLEAN`, strings to `TEXT`, dates to `DATE`, datetimes to `DATETIME`, and
  arrays/fallback values to compact JSON `TEXT`. Dynamic attributes stay
  physically nullable. A `uint64` field exceeding SQLite's signed integer
  range uses the lossless fallback.
- GLB maps explicit numeric widths to matching `EXT_structural_metadata`
  component types. It supports native one-dimensional arrays of simple
  scalar/string items. Unsupported nested shapes use compact JSON strings.
  Missing numeric values use a collision-free sentinel. The same encoding
  applies when Tyler uses the GLB as 3D Tiles content.
- OBJ does not invoke attribute projection. If the user supplies a valid
  specification for OBJ output, the CLI warns once and ignores it.

### Diagnostics

Specification file I/O, JSON syntax, or semantic validation errors are fatal
and occur before an output path is created or modified. Input/schema
mismatches warn and continue.

Library code returns structured diagnostics. Tyler and `cjconvert` aggregate
them by format, location, attribute, and category and print deterministic
summaries. Categories include undeclared attributes, unexpected locations,
missing expected attributes, type/coercion conflicts, non-nullable values,
invalid temporal values, and unlisted categorical values.

Default output contains one warning per aggregate key with the occurrence
count. `--verbose` prints the first 100 per-value warning details globally,
emits one cap notice, and then continues aggregate counting only. Tyler merges
reports from all workers and prints once for the complete run.

## Consequences

Positive:

- One reusable attribute domain replaces duplicated inference and conflict
  behavior.
- The same specification controls CityJSON, tabular, and structural-metadata
  outputs without imposing 3DBAG's layer topology.
- Attribute selection becomes explicit per format and native location.
- Numeric width is part of the contract, allowing compact GLB metadata without
  silently losing 64-bit values.
- Tyler remains input-driven and can evolve toward richer ETL transformations
  independently of physical writers.

Trade-offs:

- The feature spans the 3DBAG specification, `cityjson-rs`, and Tyler and must
  be delivered as coordinated breaking changes.
- Tyler needs a bounded-memory preflight pass before parallel export to ensure
  dataset-wide field decisions.
- `cityjson-attributes` adds one crate and one release artifact to the
  `cityjson-rs` workspace.
- Schema-provided attribute selection replaces `--object-attributes`; existing
  command lines must migrate.

Neutral:

- Producing the historic 3DBAG GeoPackage layer topology remains a separate
  input reshaping/layout concern.
- A future shared tabular or PostgreSQL abstraction can build on
  `cityjson-attributes`, but is not designed by this decision.
