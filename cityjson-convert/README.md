# cityjson-convert

`cityjson-convert` converts CityJSON data to other formats.

It supports GLB, CityJSON and CityJSONSeq output, both through the library API
and through the `cjconvert` CLI.

By default, source builds use the `proj-system` feature, which enables
PROJ-backed clipping and reprojection through `cityjson-lib` while preferring a
system PROJ installation. Use `--no-default-features --features proj-bundled`
to build with bundled PROJ source support.

## Library

Use the library when you want to convert parsed `CityModel` values to another
format.

### GLB

```rust
use cityjson_convert::{convert_to_glb, ExportOptions, GeometryPlacement};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    let options = ExportOptions {
        native_glb_color: "#FFC0CB".to_string(),
        metadata_class_name: "cityobject".to_string(),
        feature_type_colors: Default::default(),
        geometry_placement: GeometryPlacement::SourceCoordinates,
        clip_bbox: None,
        clip_geographic_region: None,
        smooth_normals: false,
        quantize_geometry: true,
        meshopt_compression: true,
    };

    convert_to_glb(model, "output/tiles.glb", &options)
}
```

### CityJSON

```rust
use cityjson_convert::{convert_to_cityjson, JsonExportOptions};
use cityjson_lib::CityModel;

fn export(model: &CityModel) -> anyhow::Result<()> {
    convert_to_cityjson(model, "output/model.city.json", &JsonExportOptions::default())
}
```

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

Convert a CityJSONSeq or CityJSONFeature stream to CityJSONSeq:

```shell
cjconvert input.city.jsonl --output output/features.city.jsonl --format cityjsonseq
```

`--format cityjsonseq` accepts CityJSONSeq/CityJSONFeature-stream input. A
single CityJSON document is rejected because decomposing a merged document into
feature models is not implemented.

You can also customize the output:

```shell
cjconvert input.city.json \
  --output output/model.glb \
  --native-glb-color "#FFC0CB" \
  --3dtiles-metadata-class cityobject \
  --smooth-normals
```
