# Define Shared CityJSON Tabular Projection Schema

## Status

Proposed

## Context

Tyler v1.0 introduces multiple tabular or partly tabular outputs for
CityJSON-compatible data. CSV and TSV are geometry-less tabular outputs.
GeoPackage stores feature geometries with relational attributes. Future
CityArrow and CityParquet outputs can preserve more nested structure while
serving the same analytical use cases.

Without a shared schema decision, each writer could infer CityJSON attributes
independently. That would make the same selected CityObjects expose different
non-geometry columns in CSV, TSV, GeoPackage, CityArrow, and CityParquet. It
would also encourage older GIS-oriented behavior where any complex attribute is
serialized as JSON text, even when it has a stable typed structure that can be
projected into queryable fields.

The `cityjson-arrow` projection semantics provide the right logical direction:
infer a stable schema from CityJSON values, preserve typed nested data where it
is compatible, widen numeric values consistently, and use an explicit JSON
fallback only for heterogeneous or otherwise incompatible paths.

## Decision

Tyler will define one logical CityJSON tabular projection schema shared by all
tabular outputs. CSV, TSV, GeoPackage, and future CityArrow or CityParquet
writers encode that logical projection in format-appropriate physical forms.

The shared projection applies to CityModel, CityObject, Geometry, and Semantic
tabular projections. All of them use the same attribute projection rules; only
their row identity, geometry handling, and physical encoding differ by output
format.

### Logical Projection Namespaces

The logical schema contains these namespaces when the source data and selected
output expose them:

- `citymodel` or `metadata` attributes for geometry-less model-level output
- `cityobject.attributes`
- `cityobject.extra`
- `geometry.extra`, if Tyler exposes geometry-level extra data
- `semantic.attributes`
- `address`, when present

Each namespace projects its values recursively. Nested objects become logical
struct fields, not JSON text by default.

### Logical Value Vocabulary

The logical value vocabulary is aligned with `cityjson-arrow`:

| Logical value | Meaning |
|---------------|---------|
| `Null` | Null-only path |
| `Boolean` | Boolean values |
| `UInt64` | Non-negative integer values that require unsigned range |
| `Int64` | Signed integer values |
| `Float64` | Floating-point values or integer/float widened numeric values |
| `Utf8` | String values |
| `GeometryRef` | Geometry references where a projection needs them |
| `List<T>` | Homogeneous lists with a projected item type |
| `Struct{...}` | Objects with recursively projected fields |
| `Json` | Explicit fallback after inference finds heterogeneous or incompatible values |

Numeric inference follows the same widening direction as `cityjson-arrow`.
Compatible integer paths stay integer where possible. Mixed integer and
floating-point paths widen to `Float64`. Heterogeneous paths that cannot be
represented by a stable typed projection become `Json`.

Lists are projected only when they have a usable homogeneous item policy for
the selected logical path. Otherwise, the path uses the explicit `Json`
fallback. The fallback is a logical category, not permission for every writer to
blindly serialize all arrays and objects.

### Physical Encodings

Physical formats encode the same logical fields differently:

- CSV and TSV use flattened columns. JSON text is used only for logical `Json`
  fields or for unsupported list/object physical cells.
- GeoPackage uses flattened SQLite columns plus its `geom` column. The
  non-geometry logical fields match geometry-less tabular output for the same
  selected CityObjects.
- CityArrow and CityParquet use nested Arrow or Parquet fields where supported,
  plus any transport or package metadata needed by those formats.

Flattened physical encodings must use deterministic names for nested paths.
For example, `cityobject.attributes.metrics.height` can become a column such as
`attributes__metrics__height` in a GeoPackage feature table. Name conflicts are
escaped or prefixed deterministically by the writer.

Addresses use the same projection rules as other attributes. A stable address
object can therefore produce typed projected fields. It only becomes compact
JSON text when the address path is inferred as logical `Json` or when a physical
format cannot represent a projected field directly.

## Consequences

Good:

- CSV, TSV, GeoPackage, and future CityArrow or CityParquet outputs expose the
  same logical non-geometry fields for the same selected CityObjects.
- Nested attributes such as `{ "metrics": { "height": 12.5 } }` become typed,
  queryable projected fields instead of JSON blobs.
- Numeric widening and nullable behavior are consistent across tabular outputs.
- GeoPackage remains QGIS-friendly through flattened physical columns while
  sharing the same logical schema as geometry-less tabular output.

Trade-offs:

- Writers need a shared projection inference layer instead of local one-off
  JSON-to-column mappings.
- Flattened formats need deterministic path naming and conflict handling.
- Some list and heterogeneous value paths still become JSON text, but only
  through the explicit `Json` fallback category.

Neutral:

- This schema is an export projection, not a CityJSON round-trip guarantee.
- Physical writers may add format-specific identity, geometry, relation, or
  package metadata columns outside the shared logical attribute projection.

## Validation Plan

Documentation validation:

- search ADRs for stale wording such as "complex attributes are compact JSON
  text" and remove or qualify it
- ensure ADR 011 links both ADR 012 and this ADR where relevant
- ensure ADR 012 no longer defines a competing standalone attribute inference
  system

Future implementation validation:

- a nested object attribute like `{ "metrics": { "height": 12.5 } }` produces a
  typed projected field, not a JSON blob
- homogeneous scalar attributes map to typed columns in CSV, TSV, and
  GeoPackage
- mixed integer and floating-point numeric paths follow the same widening rules
  as `cityjson-arrow`
- incompatible heterogeneous values map to the explicit `Json` fallback
  consistently across outputs
- GeoPackage feature tables and geometry-less tabular output expose the same
  non-geometry attribute columns for the same selected CityObjects
