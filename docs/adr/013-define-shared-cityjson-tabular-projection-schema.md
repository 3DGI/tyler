# Define Shared CityJSON Tabular Projection Schema

## Status

Proposed

## Context

Tyler v1.0 introduces multiple tabular or partly tabular outputs for
CityJSON-compatible data. CSV and TSV are geometry-less tabular outputs.
GeoPackage stores feature geometries with relational attributes.
`cityjson-rs` already includes `cityjson-arrow` and `cityjson-parquet` crates
that convert CityModels into canonical nested and tabular Arrow or Parquet
encodings for analytical use cases.

Tyler needs a concrete tabular contract so CSV, TSV, and GeoPackage do not each
invent their own attribute columns. The contract should follow the existing
`cityjson-arrow` direction: typed nested values remain typed when they have a
stable shape, numeric values widen consistently, and JSON text is an explicit
fallback for heterogeneous or otherwise incompatible paths.

## Decision

Tyler will use the logical row schemas below for tabular outputs. CSV, TSV, and
GeoPackage encode the same logical fields in format-appropriate physical forms.
The schema aligns with the existing CityArrow and CityParquet projection
semantics in `cityjson-rs`; it does not define those formats as future
Tyler-only work.

### Logical Row Schemas

Tyler's shared projection defines these logical row types:

#### `cityobjects`

One row per selected CityObject. This is the base schema for CSV and TSV.

| Logical field     | Type                    | Required | Source                                                              |
|-------------------|-------------------------|----------|---------------------------------------------------------------------|
| `cityobject_id`   | `Utf8`                  | yes      | CityJSON object key                                                 |
| `cityobject_ix`   | `UInt64`                | yes      | stable export ordinal                                               |
| `cityobject_type` | `Utf8`                  | yes      | CityObject `type`                                                   |
| `bbox`            | `FixedList<Float64, 6>` | no       | computed or source bbox, if available                               |
| `attributes`      | `Struct{...}`           | no       | CityObject `attributes`                                             |
| `extra`           | `Struct{...}`           | no       | extension/custom CityObject members not represented by fixed fields |

The fixed fields match the CityArrow `cityobjects` table conceptually. Tyler
may omit `cityobject_ix` from user-facing CSV/TSV only when no downstream join
or relation output needs it; internally it remains the stable row identity.

#### `features`

One row per exported feature geometry. This is the base schema for GeoPackage
feature layers and extends `cityobjects` with geometry identity.

| Logical field     | Type            | Required | Source                                                  |
|-------------------|-----------------|----------|---------------------------------------------------------|
| `feature_id`      | `UInt64`        | yes      | stable feature row ordinal                              |
| `cityobject_id`   | `Utf8`          | yes      | owning CityObject id                                    |
| `cityobject_ix`   | `UInt64`        | yes      | owning CityObject ordinal                               |
| `cityobject_type` | `Utf8`          | yes      | owning CityObject type                                  |
| `geometry_ix`     | `UInt64`        | yes      | geometry index within the CityObject                    |
| `geometry_type`   | `Utf8`          | yes      | CityJSON geometry type                                  |
| `lod`             | `Utf8`          | no       | CityJSON LoD value                                      |
| `geom`            | format geometry | yes      | physical geometry column, if the format stores geometry |
| `attributes`      | `Struct{...}`   | no       | owning CityObject `attributes`                          |
| `extra`           | `Struct{...}`   | no       | owning CityObject extension/custom members              |
| `geometry_extra`  | `Struct{...}`   | no       | geometry extension/custom members, if exposed           |

GeoPackage uses this row schema with the physical columns required by
[ADR 012](012-define-cityjson-to-geopackage-schema-mapping.md), including
`id INTEGER PRIMARY KEY` and `geom`.

#### `semantic_surfaces`

One row per exported semantic surface when Tyler splits semantics.

| Logical field   | Type            | Required | Source                                     |
|-----------------|-----------------|----------|--------------------------------------------|
| `semantic_id`   | `UInt64`        | yes      | stable semantic row ordinal                |
| `cityobject_id` | `Utf8`          | yes      | owning CityObject id                       |
| `cityobject_ix` | `UInt64`        | yes      | owning CityObject ordinal                  |
| `geometry_ix`   | `UInt64`        | yes      | owning geometry index                      |
| `surface_ix`    | `UInt64`        | yes      | semantic surface index within the geometry |
| `semantic_type` | `Utf8`          | yes      | semantic object `type`                     |
| `geom`          | format geometry | yes      | semantic geometry, if materialized         |
| `attributes`    | `Struct{...}`   | no       | semantic surface `attributes`              |

### Projected Field Types

Nested CityJSON values are projected recursively into this value vocabulary,
aligned with `cityjson-arrow`:

| Logical type      | Meaning                                   | Flat physical type                                                 |
|-------------------|-------------------------------------------|--------------------------------------------------------------------|
| `Null`            | null-only path                            | no standalone column unless needed for a nullable typed path       |
| `Boolean`         | boolean values                            | `BOOLEAN` or text `true`/`false`                                   |
| `UInt64`          | non-negative integers                     | unsigned integer where available, otherwise integer/text if needed |
| `Int64`           | signed integers                           | integer                                                            |
| `Float64`         | floats or widened int/float paths         | real number                                                        |
| `Utf8`            | string values                             | text                                                               |
| `FixedList<T, N>` | fixed-length homogeneous list             | native fixed list if available, otherwise compact JSON text        |
| `List<T>`         | homogeneous list values                   | native list if available, otherwise compact JSON text              |
| `Struct{...}`     | object values with projected child fields | flattened child columns                                            |
| `Json`            | heterogeneous or incompatible values      | compact JSON text                                                  |

Numeric inference follows the same widening direction as `cityjson-arrow`.
Compatible integer paths stay integer where possible. Mixed integer and
floating-point paths widen to `Float64`. Heterogeneous paths that cannot be
represented by a stable typed projection become `Json`.

Lists are projected as `List<T>` only when all present list values have a
homogeneous item type after numeric widening. Lists with incompatible item
types, mixed scalar/list/object shapes, or unstable nested list shapes become
`Json`.

### Flattened Column Names

Flat formats use the field namespace plus the nested path, joined by `__`.
Fixed identity columns are not prefixed.

| Logical path                | Flat column                   |
|-----------------------------|-------------------------------|
| `attributes.measuredHeight` | `attributes__measuredHeight`  |
| `attributes.metrics.height` | `attributes__metrics__height` |
| `extra.creationDate`        | `extra__creationDate`         |
| `geometry_extra.source.id`  | `geometry_extra__source__id`  |
| `attributes.address.street` | `attributes__address__street` |

Writers must escape path segments deterministically before joining them. The
escape rule is: replace `%` with `%25` and replace literal `__` with `%5F%5F`.
If the resulting column conflicts with a fixed column or another projected
path, append `__2`, `__3`, and so on until it is unique.

### Physical Encodings

Physical formats encode the same logical fields differently:

- CSV and TSV write one `cityobjects` row per selected CityObject. They use
  flattened columns. JSON text is used only for logical `Json` fields or for
  logical `List<T>` and `FixedList<T, N>` fields because CSV/TSV have no native
  list cell type.
- GeoPackage uses flattened SQLite columns plus its `geom` column. The
  non-geometry logical fields in each feature layer match the `cityobjects`
  projection for the owning CityObject. `List<T>` and `FixedList<T, N>` fields
  become compact JSON `TEXT` unless ADR 012 defines a row-expansion policy for
  that path.
- The existing CityArrow and CityParquet encodings in `cityjson-rs` use nested
  Arrow or Parquet fields where supported, plus any transport or package
  metadata needed by those formats.

### Examples

Given these two CityObjects:

```json
{
  "b-1": {
    "type": "Building",
    "attributes": {
      "measuredHeight": 12,
      "name": "Library",
      "metrics": {
        "roofSlope": 0.25
      },
      "tags": [ "public", "education" ],
      "mixed": [ 607, false, 28.47 ]
    }
  },
  "b-2": {
    "type": "Building",
    "attributes": {
      "measuredHeight": 14.5,
      "name": null,
      "metrics": {
        "roofSlope": 0
      },
      "tags": [ "office" ]
    }
  }
}
```

The inferred logical `attributes` struct is:

```text
attributes: Struct{
  measuredHeight: Float64,
  name: Utf8,
  metrics: Struct{
    roofSlope: Float64
  },
  tags: List<Utf8>,
  mixed: Json
}
```

CSV/TSV columns become:

| cityobject_id | cityobject_type | attributes__measuredHeight | attributes__name | attributes__metrics__roofSlope | attributes__tags         | attributes__mixed   |
|---------------|-----------------|----------------------------|------------------|--------------------------------|--------------------------|---------------------|
| `b-1`         | `Building`      | `12.0`                     | `Library`        | `0.25`                         | `["public","education"]` | `[607,false,28.47]` |
| `b-2`         | `Building`      | `14.5`                     | null             | `0.0`                          | `["office"]`             | null                |

GeoPackage feature layers use the same non-geometry columns on each exported
geometry row, plus their required feature columns:

| id  | cityobject_id | cityobject_type | geometry_type | lod   | attributes__measuredHeight | attributes__metrics__roofSlope | geom                 |
|-----|---------------|-----------------|---------------|-------|----------------------------|--------------------------------|----------------------|
| `1` | `b-1`         | `Building`      | `Solid`       | `2.2` | `12.0`                     | `0.25`                         | `MultiPolygonZ(...)` |

CityArrow/CityParquet keep the same projection nested instead of flattening it:

```text
cityobjects.attributes.measuredHeight: Float64
cityobjects.attributes.metrics.roofSlope: Float64
cityobjects.attributes.tags: List<Utf8>
cityobjects.attributes.mixed: Json
```

Addresses use the same projection rules as other attributes. A stable address
object can therefore produce typed projected fields. It only becomes compact
JSON text when the address path is inferred as logical `Json` or when a physical
format cannot represent a projected field directly.

#### 3DBAG

The 3DBAG CityJSON file looks like this: https://gist.github.com/balazsdukai/e1c8a32ec7933ded2e29cf75d47cea5f

Default output:

| cityobject\_id                   | cityobject\_type | attributes\_\_identificatie    | attributes\_\_oorspronkelijkbouwjaar | attributes\_\_status | attributes\_\_b3\_dak\_type | attributes\_\_b3\_h\_dak\_max |
|:---------------------------------|:-----------------|:-------------------------------|:-------------------------------------|:---------------------|:----------------------------|:------------------------------|
| NL.IMBAG.Pand.0935100000021359-0 | BuildingPart     | null                           | null                                 | null                 | null                        | null                          |
| NL.IMBAG.Pand.0935100000021359   | Building         | NL.IMBAG.Pand.0935100000021359 | 1975                                 | Pand in gebruik      | slanted                     | 66.181999                     |

Output with `split-semantics` enabled creates a separate semantics table:

| semantic\_id | cityobject\_id                   | cityobject\_ix | geometry\_ix | surface\_ix | semantic\_type | attributes\_\_on\_footprint\_edge | attributes\_\_b3\_h\_dak\_50p | attributes\_\_b3\_h\_dak\_70p | attributes\_\_b3\_h\_dak\_max | attributes\_\_b3\_h\_dak\_min | attributes\_\_b3\_azimut | attributes\_\_b3\_hellingshoek |
|:-------------|:---------------------------------|:---------------|:-------------|:------------|:---------------|:----------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:-------------------------|:-------------------------------|
| 0            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 0            | 0           | GroundSurface  | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 1            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 0            | 1           | WallSurface    | true                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 2            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 0            | 2           | WallSurface    | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 3            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 0            | 3           | RoofSurface    | null                              | 63.628269                     | 64.645882                     | 66.424767                     | 60.808277                     | null                     | null                           |
| 4            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 1            | 0           | GroundSurface  | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 5            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 1            | 1           | WallSurface    | true                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 6            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 1            | 2           | WallSurface    | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 7            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 1            | 3           | RoofSurface    | null                              | 63.628269                     | 64.645882                     | 66.424767                     | 60.808277                     | null                     | null                           |
| 8            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 2            | 0           | GroundSurface  | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 9            | NL.IMBAG.Pand.0935100000021359-0 | 0              | 2            | 1           | WallSurface    | true                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 10           | NL.IMBAG.Pand.0935100000021359-0 | 0              | 2            | 2           | WallSurface    | null                              | null                          | null                          | null                          | null                          | null                     | null                           |
| 11           | NL.IMBAG.Pand.0935100000021359-0 | 0              | 2            | 3           | RoofSurface    | null                              | 63.558567                     | 64.689850                     | 66.424767                     | 60.808277                     | 263.904449               | 39.630585                      |
| 12           | NL.IMBAG.Pand.0935100000021359-0 | 0              | 2            | 4           | RoofSurface    | null                              | 63.826801                     | 64.835037                     | 66.424767                     | 61.332588                     | 83.800995                | 39.800274                      |

**Notes on TSV output**

- CityObject and Semantics without values in any of the attributes are excluded from the output by default. The flag `--tsv-include-null-rows` includes all objects in the output, even those without
  attributes values.
- CityObject and Semantic hierarchy is opt-in with the `--tsv-include-hierarchy` flag.
- The `--tsv-include-cityjson-ordinal` flag adds the `cityjson_ix` field to the TSV output.
- Metadata is opt-in with `--tsv-include-metadata` flag. When set, a separate file with just the Metadata fields is written. Tyler appends the Metadata of each TSV-tile to a single Metadata file, and
  includes the WKT of the file extent for each tile.

## Consequences

Good:

- Tyler's CSV, TSV, and GeoPackage outputs expose logical non-geometry fields
  that align with the existing CityArrow and CityParquet schema direction for
  the same selected CityObjects.
- Nested attributes such as `{ "metrics": { "height": 12.5 } }` become typed,
  queryable projected fields instead of JSON blobs.
- Numeric widening and nullable behavior are consistent across tabular outputs.
- GeoPackage remains QGIS-friendly through flattened physical columns while
  sharing the same logical schema as geometry-less tabular output.

Trade-offs:

- Writers need a shared projection inference layer instead of local one-off
  JSON-to-column mappings.
- Flattened formats need deterministic path naming and conflict handling.
- Flat formats may write valid `List<T>` and `FixedList<T, N>` values as JSON
  text because they have no native list cell type.
- Heterogeneous value paths still become JSON text, but only through the
  explicit `Json` fallback category.

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
