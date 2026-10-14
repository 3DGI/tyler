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
        voxel_remesh: None,
    };

    let report = convert_to_glb(model, "output/tiles.glb", &options)?;
    println!("GLB: {} bytes", report.output_bytes);
    Ok(())
}
```

### Experimental voxel remeshing

Set `ExportOptions.voxel_remesh` to `Some(VoxelRemeshOptions { geometric_error: 10.0 })`
to merge all supplied geometry into one voxel reconstruction, then simplify it.
Input and error units must be Cartesian metres. Select source geometry/LoDs before
calling the library; it does not select them automatically. Both clipping options
are incompatible with remeshing.

The optional `--lod` flag selects an exact source LoD for every output format:

```shell
cjconvert buildings.city.json --output proxy.glb --lod 2.2 --remesh-error 10
```

Without `--lod`, all source LoDs are exported, including when remeshing.
Selection preserves object attributes and relationships, including geometry-less
parents. A LoD that matches no geometry in the input returns an error.
`--remesh-error` is available only for GLB output and does not require `--lod`.
The CLI reports rectangular grid dimensions,
voxel size, error components, triangle counts and GLB bytes to stderr. The voxel
edge length is `geometric_error / (2 * sqrt(3))`, rounded toward finer quality.
Each grid axis is sized independently, with empty padding; no resolution clamping
is applied. One grid merges all supplied geometry, including different object
types. Memory grows with grid volume, so finer voxels over large areas require more
memory. Unrepresentable grids return an error instead of changing the voxel size.

Half the requested budget controls voxel reconstruction (using one voxel diagonal
as its allowance), and half controls absolute simplification error. The returned
`GlbExportReport.remeshing` includes `estimated_geometric_error`: voxel diagonal +
simplifier-reported error + quantization/float-rounding allowance. This is a
selection estimate, not a certified surface-distance bound; remeshing can merge
gaps, erase thin wires, and fill interiors. Simplification can also remove isolated
small buildings, making source-to-proxy distance much larger than this estimate.
Quantization and float rounding can increase the estimate beyond the requested
budget.

After simplification, each triangle receives the ID of the source triangle nearest
to its centroid. All three corners share that ID. Original metadata rows are
retained, including rows with no surviving geometry; ownership around merged
buildings is approximate. Normals are regenerated and remeshed materials use
backface culling. Quantization, compression, colours, and smooth normals remain
available. Remeshing can increase triangle count for already simple source meshes;
the report shows this explicitly.

For future 3D Tiles parent content, validate a surface-distance error before using
the proxy for the parent's `geometricError`; the current estimate alone is
insufficient when disconnected objects disappear. Keep parent errors at least as
large as child errors, use
`REPLACE` refinement, and generate each proxy from the selected source geometry.
Tyler's explicit content planner enables voxel remeshing for simplified internal
tiles using each level's planned allowance. Full-detail leaves and implicit
content use `voxel_remesh: None`. These planned tile errors remain selection
estimates; they are not replaced by certified surface-distance measurements.

See [ADR 015](../docs/adr/015-use-rust-voxel-remeshing-for-optional-glb-proxies.md)
for the implementation decision, geometric-error limitations, and alternatives.
See the [Rust remesher evaluation](../docs/performance/voxel-remeshing/rust-voxel-remeshing-evaluation.md) for
current memory/runtime comparisons. Surface-distance observations are recorded
in the [native rectangular grid evaluation](../docs/performance/voxel-remeshing/rectangular-voxel-grid-evaluation.md)
and the [original evaluation](../docs/performance/voxel-remeshing/voxel-remeshing-evaluation.md).

The reconstruction code lives in `cityjson-convert` as a Rust port of
meshoptimizer v1.3's remesher. Its MIT attribution and licence are retained in
[`src/voxel_remesh/LICENSE.md`](src/voxel_remesh/LICENSE.md). Native simplification,
remapping, optimization, and compression use upstream `meshopt` master, pinned
to `3c8c9126cd98e70f34c6d9c54d02f41af6749406`, with meshoptimizer v1.2. No local
`meshopt` fork is required. v1.2 has position-updating simplification and vertex
locks, but lacks v1.3's `PreserveFolds` flag, so the current pipeline uses
`ErrorAbsolute` alone. Updated measurements are in the
[upstream v1.2 evaluation](../docs/performance/voxel-remeshing/upstream-meshopt-1.2-evaluation.md).

The public `cityjson_convert::voxel_remesh::remesh` API takes a `VoxelGrid` (origin and voxel size) and a
`VoxelExtent` (global processing cell range), and returns shared positions,
triangle indices, grid dimensions, and queryable cell identities. The converter
uses a centred grid for the current whole-object workflow and welds coincident
positions before simplification. Packed occupancy and indexed output reduce
memory; reconstruction builds the grid once and samples the source twice.
The checked native `simplify_positions_with_update` wrapper in the same module
retains optional per-vertex locks. Reconstruction modules forbid unsafe code.

These interfaces prepare for shared-grid regional reconstruction. Shell mode can
process a cropped domain; solid mode requires complete input coverage with padding.
Clipped GLB remeshing remains unsupported. Future tiled integration needs buffered
source geometry, valid solid boundary classification, locked clip boundaries, and
consistent encoding; neighbouring tiles at different resolutions also need a
transition policy. Shared grid coordinates alone do not guarantee matching GLBs.

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
