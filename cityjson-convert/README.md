# cityjson-convert

`cityjson-convert` converts CityJSON data to other formats.

It supports GLB, OBJ, CityJSON and CityJSONSeq output, both through the library API
and through the `cjconvert` CLI.

By default, source builds use the `proj-system` feature, which enables
PROJ-backed clipping and reprojection through `cityjson-lib` while preferring a
system PROJ installation with native network-grid capability enabled. Use
`--no-default-features --features proj-bundled` to build with bundled PROJ
source support; bundled mode is separate from native PROJ networking.

## Library

Use the library when you want to convert parsed `CityModel` values to another
format.

### GLB

```rust
use cityjson_convert::{convert_to_glb, ExportOptions, GeometryPlacement, GeometryRepresentation};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    let options = ExportOptions {
        native_glb_color: "#FFC0CB".to_string(),
        metadata_class_name: "cityobject".to_string(),
        feature_type_colors: Default::default(),
        geometry_placement: GeometryPlacement::SourceCoordinates,
        representation: GeometryRepresentation::FullDetail,
        clip_bbox: None,
        clip_geographic_region: None,
        smooth_normals: false,
        quantize_geometry: true,
        meshopt_compression: true,
    };

    convert_to_glb(model, "output/tiles.glb", &options)
}
```

`GeometryRepresentation::Simplified { max_error_meters }` is the typed input
contract for callers that request a lower-detail mesh. The current GLB writer
validates the allowance but still writes full-detail geometry; mesh
simplification is implemented separately.

### CityJSON

```rust
use cityjson_convert::{convert_to_cityjson, JsonExportOptions};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    convert_to_cityjson(model, "output/model.city.json", &JsonExportOptions::default())
}
```

### OBJ

```rust
use cityjson_convert::{convert_to_obj, ObjExportOptions};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    convert_to_obj(model, "output/model.obj", &ObjExportOptions::default())
}
```

OBJ output is geometry-only Wavefront OBJ grouped by CityObject ID. Coordinates
are written in the source CityJSON coordinate space.

### GeoPackage

```rust
use cityjson_convert::{convert_to_gpkg, GpkgExportOptions};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    convert_to_gpkg(model, "output/model.gpkg", &GpkgExportOptions::default())
}
```

GeoPackage output requires `metadata.referenceSystem` to contain a parseable EPSG
identifier. The converter rejects missing or non-EPSG CRS metadata before it
changes the output path; it intentionally has no output-CRS override.

### CityJSONSeq

```rust
use cityjson_convert::{convert_to_cityjsonseq, CityJsonSeqExportOptions};
use cityjson_lib::CityModel;

fn export(base_root: &CityModel, feature_models: &[CityModel]) -> anyhow::Result<()> {
    convert_to_cityjsonseq(
        base_root,
        feature_models,
        "output/features.city.jsonl",
        &CityJsonSeqExportOptions::default(),
    )
}
```

`convert_to_cityjsonseq` expects feature models and writes a `CityJSONSeq`
stream with a `CityJSON` header followed by `CityJSONFeature` items.

## CLI

Convert a CityJSON file to a GLB file from the command line:

```shell
cjconvert input.city.json --output output/model.glb
```

The default format is `glb`. You can set it explicitly with `--format glb`:

```shell
cjconvert input.city.json --output output/model.glb --format glb
```

Convert a CityJSON file to a CityJSON file:

```shell
cjconvert input.city.json --output output/model.city.json --format cityjson
```

Convert a CityJSON file to an OBJ file:

```shell
cjconvert input.city.json --output output/model.obj --format obj
```

Convert a CityJSONSeq or CityJSONFeature stream to CityJSONSeq:

```shell
cjconvert input.city.jsonl --output output/features.city.jsonl --format cityjsonseq
```

`--format cityjsonseq` accepts CityJSONSeq/CityJSONFeature-stream input. A
single CityJSON document is rejected because decomposing a merged document into
feature models is not implemented.

The CLI also accepts a dataset directory. This includes the legacy split
CityJSONFeature layout with a `metadata.json` base document and one or more
`*.city.jsonl` feature files:

```text
dataset/
├── metadata.json
├── feature1.city.jsonl
└── feature2.city.jsonl
```

```shell
cjconvert dataset --output output/model.glb --format glb
```

Directory inputs are discovered and read through `cityjson-index`. The input 
directory is recursively scanned for CityJSONFeature files. The CLI creates or 
refreshes the `.cityjson-index.sqlite` sidecar in the dataset
directory. GLB, OBJ and CityJSON output merge all indexed features into one
model in memory; CityJSONSeq preserves them as separate features.

You can also customize the output:

```shell
cjconvert input.city.json \
  --output output/model.glb \
  --native-glb-color "#FFC0CB" \
  --3dtiles-metadata-class cityobject \
  --smooth-normals
```

### GeoPackage options

Write GeoPackage output with `--format gpkg`:

```shell
cjconvert input.city.json --output output/model.gpkg --format gpkg
```

The source `CityJSON` must provide a parseable EPSG identifier in
`metadata.referenceSystem`. All GeoPackage-specific options are disabled by
default.

| Option | Effect |
|--------|--------|
| `--gpkg-split-lod` | Writes separate feature layers for each CityObject type, geometry family, and LoD. Without this option, LoDs share a layer and are stored in its `lod` column. |
| `--gpkg-include-semantics` | Adds a `semantics` feature layer containing semantic primitive rows. |
| `--gpkg-include-hierarchy` | Adds the non-spatial `cityobject_hierarchy` and `semantic_hierarchy` tables. |
| `--gpkg-include-address` | Adds an `addresses` feature layer from `CityObject.extra.address` values. |
| `--gpkg-include-metadata` | Exports source CityJSON metadata through the GeoPackage metadata extension. |

For example, to export every optional GeoPackage product and split the main
feature layers by LoD:

```shell
cjconvert input.city.json --output output/model.gpkg --format gpkg \
  --gpkg-split-lod \
  --gpkg-include-semantics \
  --gpkg-include-hierarchy \
  --gpkg-include-address \
  --gpkg-include-metadata
```
