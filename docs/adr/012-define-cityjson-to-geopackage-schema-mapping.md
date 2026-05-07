# Define CityJSON to GeoPackage Schema Mapping

## Status

Proposed

## Context

Tyler v1.0 adds GeoPackage output for CityJSON-compatible data. The output is
primarily meant for standard GIS tools such as QGIS, so the mapping prioritizes
simple feature layers, typed attributes, concrete geometry types, and queryable
relationships. The mapping is one-way: it defines how to export CityJSON to
GeoPackage, not how to reconstruct CityJSON from GeoPackage.

This mapping is informed by:

- the CityJSON 2.0.2 specification: https://www.cityjson.org/specs/2.0.2/
- the GeoPackage Encoding Standard: https://www.geopackage.org/spec/
- the 3DCityDB v5 relational schema:
  https://docs.3dcitydb.org/1.3/3dcitydb/relational-schema/

3DCityDB v5 is useful as a conceptual reference because it separates features,
geometries, relationships, metadata, appearance, and extension data. GeoPackage
and QGIS have different constraints, so this mapping adapts those concepts to a
layer-oriented Simple Features model.

## Decision

GeoPackage export uses a GIS-first mapping. Coordinates are written as
real-world coordinates in one declared CRS. CityJSON geometry templates are
resolved into explicit geometries. Appearance and metadata are omitted by
default, except for CRS metadata. CityObject hierarchy and optional semantic
hierarchy are written as separate attributes tables.

### Coordinates and CRS

CityJSON coordinates are transformed before export:

```text
real_world_coordinate = vertex * transform.scale + transform.translate
```

The GeoPackage stores real-world XYZ coordinates directly. The CityJSON
`transform` object is not written because it has already been applied.

The GeoPackage uses one CRS for all feature layers. The CRS is stored in
`gpkg_spatial_ref_sys` and referenced by `gpkg_geometry_columns`. If the input
CRS is missing or ambiguous, the converter must require an explicit output CRS.
Use the GeoPackage WKT CRS extension when needed for an accurate CRS definition.

### Layer Model

Feature layers are homogeneous by CityObject type and geometry family:

```text
<cityobject_type>_<geometry_family>
```

Examples:

```text
building_surface
building_solid
road_surface
road_line
cityfurniture_point
solitaryvegetationobject_point
genericcityobject_surface
```

Do not create one mixed-geometry layer per CityObject type by default. Concrete
geometry types are more useful in QGIS and avoid generic `GEOMETRY` layers.

CityJSON geometries map to GeoPackage feature geometries as follows:

| CityJSON geometry type | Geometry family        | GeoPackage geometry type |
|------------------------|------------------------|--------------------------|
| `MultiPoint`           | `point`                | `MultiPointZ`            |
| `MultiLineString`      | `line`                 | `MultiLineStringZ`       |
| `MultiSurface`         | `surface`              | `MultiPolygonZ`          |
| `CompositeSurface`     | `surface`              | `MultiPolygonZ`          |
| `Solid`                | `solid`                | `MultiPolygonZ`          |
| `MultiSolid`           | `solid`                | `MultiPolygonZ`          |
| `CompositeSolid`       | `solid`                | `MultiPolygonZ`          |
| `GeometryInstance`     | resolved target family | resolved target type     |

`Solid`, `MultiSolid`, and `CompositeSolid` are exported as boundary
`MultiPolygonZ` geometries.

### LoD Handling

By default, all LoDs for the same CityObject type and geometry family are stored
in the same layer with a `lod` column.

If `--gpkg-split-lod` is set, split layers by LoD:

```text
<cityobject_type>_<geometry_family>_lod<lod>
```

Examples:

```text
building_surface_lod2
building_solid_lod2
road_surface_lod1
```

LoD values are stored as text because CityJSON LoDs can be values such as
`"1"`, `"2"`, `"2.2"`, or extension-specific strings.

### Feature Rows

Each CityJSON Geometry Object becomes one GeoPackage feature row. If a
CityObject has multiple geometries, each geometry is exported as a separate row
in the appropriate layer.

Every feature layer includes these columns:

| Column            | Type                    | Description                           |
|-------------------|-------------------------|---------------------------------------|
| `id`              | `INTEGER PRIMARY KEY`   | GeoPackage row id                     |
| `cityobject_id`   | `TEXT NOT NULL`         | Source CityJSON CityObject id         |
| `cityobject_type` | `TEXT NOT NULL`         | Source CityJSON CityObject type       |
| `geometry_type`   | `TEXT NOT NULL`         | Original CityJSON geometry type       |
| `lod`             | `TEXT`                  | CityJSON LoD value                    |
| `geom`            | layer-specific geometry | Registered GeoPackage geometry column |

The registered geometry column in `gpkg_geometry_columns` is always `geom`.

### Attribute Mapping

CityObject attributes are copied onto every exported geometry row for that
CityObject. Attribute columns are inferred per output layer.

| JSON value | GeoPackage type     |
|------------|---------------------|
| string     | `TEXT`              |
| integer    | `INTEGER`           |
| number     | `REAL`              |
| boolean    | `BOOLEAN`           |
| null       | `NULL`              |
| array      | compact JSON `TEXT` |
| object     | compact JSON `TEXT` |

If an attribute is absent for a row, write `NULL`. If an attribute name conflicts
with a required column, prefix the attribute column with `attr_`.

CityJSON extension attributes follow the same rules. Scalar extension values
become normal columns. Complex extension values are compact JSON text.
Extension-defined CityObject types are exported like normal CityObject types.

CityJSON addresses are stored as compact JSON text in an `address_json` column.
No address subtable is created by default because CityJSON addresses are JSON
objects rather than the strict relational address schema used by 3DCityDB.

### CityObject Hierarchy

CityObject `parents` and `children` relationships are written to a non-spatial
attributes table:

```text
cityobject_relations
```

Columns:

| Column      | Type                  | Description          |
|-------------|-----------------------|----------------------|
| `id`        | `INTEGER PRIMARY KEY` | Row id               |
| `parent_id` | `TEXT NOT NULL`       | Parent CityObject id |
| `child_id`  | `TEXT NOT NULL`       | Child CityObject id  |

Create one row per parent-child relationship. Derive missing inverse
relationships when only `parents` or only `children` is present, and avoid
duplicates. Do not require foreign keys to feature layers because one CityObject
can produce rows in multiple layers or no geometry rows at all.

Register this table in `gpkg_contents` with `data_type = 'attributes'`.

### Semantics

Semantic objects are dropped by default.

If `--gpkg-split-semantics` is enabled, semantic primitives are exported as
separate homogeneous feature layers:

```text
semantic_<semantic_type>_<geometry_family>
```

Examples:

```text
semantic_wall_surface
semantic_roofsurface_surface
semantic_trafficarea_surface
semantic_auxiliarytrafficarea_surface
```

Each semantic feature row represents one semantic primitive or one grouped set
of primitives sharing the same semantic object.

Semantic feature layers include:

| Column          | Type                    | Description                           |
|-----------------|-------------------------|---------------------------------------|
| `id`            | `INTEGER PRIMARY KEY`   | GeoPackage row id                     |
| `semantic_id`   | `TEXT NOT NULL`         | Stable generated semantic id          |
| `semantic_type` | `TEXT NOT NULL`         | CityJSON semantic object type         |
| `cityobject_id` | `TEXT NOT NULL`         | Source CityObject id                  |
| `geom`          | layer-specific geometry | Registered GeoPackage geometry column |

Semantic scalar properties are flattened into columns using the same attribute
rules as CityObjects. CityJSON semantic properties are scalar values.

When semantic splitting is enabled, create a non-spatial attributes table:

```text
semantic_relations
```

Columns:

| Column      | Type                  | Description                                      |
|-------------|-----------------------|--------------------------------------------------|
| `id`        | `INTEGER PRIMARY KEY` | Row id                                           |
| `parent_id` | `TEXT NOT NULL`       | Parent semantic id                               |
| `child_id`  | `TEXT NOT NULL`       | Child semantic id                                |

Register this table in `gpkg_contents` with `data_type = 'attributes'`.

### Geometry Templates

Resolve CityJSON `GeometryInstance` objects before writing GeoPackage features.
The GeoPackage output does not contain reusable geometry templates or
3DCityDB-style implicit geometry tables.

### Appearance

Drop CityJSON appearance by default, including materials, textures, texture
vertices, material assignments, and texture assignments.

Optional style export may generate QGIS layer styles from CityJSON material
colors, semantic surface types, CityObject types, or CLI color options such as
`--color-*`. Style generation must not change the feature schema.

### Metadata

Drop CityJSON metadata by default, except for CRS information embedded in the
GeoPackage.

If `--gpkg-include-metadata` is enabled, store the source CityJSON `metadata`
object as compact JSON using the GeoPackage metadata extension. The metadata row
references the whole GeoPackage.

Reference: https://www.geopackage.org/spec/#metadata_example_appendix

### GeoPackage Metadata

The converter creates the standard GeoPackage metadata needed for all output
layers:

- `gpkg_spatial_ref_sys`
- `gpkg_contents`
- `gpkg_geometry_columns`

Non-spatial relation tables are registered in `gpkg_contents` as attributes
tables. Optional metadata extension tables are created only when metadata export
is explicitly enabled.

## CLI Options

This mapping defines the behavior of these GeoPackage-specific options:

| Option                    | Effect                                                            |
|---------------------------|-------------------------------------------------------------------|
| `--gpkg-split-lod`        | Split feature layers by LoD                                       |
| `--gpkg-split-semantics`  | Export semantic primitives as separate feature layers             |
| `--gpkg-include-metadata` | Store CityJSON metadata using the GeoPackage metadata extension   |
| `--color-*`               | Generate optional QGIS style definitions when style export exists |

## Consequences

Good:

- QGIS receives concrete, homogeneous feature layers.
- The output keeps useful source identifiers through `cityobject_id`,
  `geometry_type`, and `lod`.
- Solids remain visible and queryable as 3D polygonal boundary geometries.
- CityObject hierarchy and optional semantic hierarchy remain queryable without
  forcing mixed geometry or complex extension tables.

Trade-offs:

- The output is intentionally not a CityJSON preservation format.
- Appearance, metadata, geometry templates, and semantic objects are simplified,
  resolved, or omitted unless explicitly requested.
- Solid topology is not preserved as native solid topology because GeoPackage
  and QGIS workflows are surface-feature oriented.
- Splitting layers by geometry family can create more layers than a one-table
  per CityObject type mapping.

## Validation Plan

A generated GeoPackage is valid for this mapping when:

- it passes SQLite integrity checks
- GeoPackage metadata tables are internally consistent
- every feature layer has one registered `geom` geometry column
- every feature layer uses a concrete geometry type
- QGIS opens all feature layers without manual configuration
- all exported coordinates are real-world coordinates in the declared CRS
- CityObject parent-child relationships are queryable through
  `cityobject_relations`
- semantic relationships are queryable through `semantic_relations` when
  `--gpkg-split-semantics` is enabled
