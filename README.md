# tyler

<p align="center">
  <img src="https://github.com/3DGI/tyler/blob/master/tyler.png" />
</p>

*tyler* creates tiles from 3D city objects.

As input, *tyler* takes a single directory, containing regular `CityJSON`, `CityJSONSeq` or the legacy `CityJSONFeature`-files.
Tyler no longer walks raw `metadata.city.json` plus [CityJSON Features](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature) trees on its own, but it
uses [cityjson-index](https://github.com/3DGI/cityjson/tree/main/crates/cityjson-index) for input resolution.

As output, *tyler* can create:

- [3D Tiles v1.1](https://docs.ogc.org/cs/22-025r4/22-025r4.html)
- tiled `CityJSON`
- tiled `CityJSONSeq`

Details of the 3D Tiles output:

- The tileset content is binary glTF (.glb).
- The glTF assets contain feature metadata (per CityObject), using the [EXT_mesh_features](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_mesh_features)
  and [EXT_structural_metadata](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata) extensions.
- The features are colored to default values, and the colors can by set per CityObject type.
- The glTF files are compressed, using the [KHR_mesh_quantization](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Khronos/KHR_mesh_quantization)
  and [EXT_meshopt_compression](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Vendor/EXT_meshopt_compression) extensions.
- Implicit tiling is supported (optional).

Additional information about the internals of *tyler* you will find in the [design document](https://github.com/3DGI/tyler/blob/master/docs/design_document.md).

## Installation

Pull the docker image with `docker pull 3dgi/tyler:<version>`, e.g. `docker pull 3dgi/tyler:0.4.1-beta1`.

### Compiling from source

*tyler* is written in Rust and you need the [Rust toolchain](https://www.rust-lang.org/learn/get-started) to compile it.

After downloading the source code from GitHub, navigate into the tyler directory and you can install *tyler* with *cargo*.

```shell
cargo install .
```

Plain `cargo build` and `cargo install .` use Tyler's default `proj-system`
feature. This enables PROJ-backed APIs and asks `proj-sys` to prefer a system
PROJ installation with native network-grid capability enabled.

`proj-sys 0.27.0` requires PROJ 9.6.2 or newer for system discovery. If no
acceptable system PROJ is found, `proj-sys` may fall back to building PROJ from
source. Build explicitly in bundled mode when you want the source build path:

```shell
cargo build --no-default-features --features proj-bundled
```

For Linux development, the recommended repeatable environment is the
devcontainer in this repository. It installs the native PROJ, SQLite, TIFF,
`pkg-config`, C/C++ build tooling and Rust tooling needed by both system and
bundled PROJ modes. This is especially useful when your host package manager
cannot provide the PROJ version required by `proj-sys`.

#### On Windows

Use [MSYS2](https://www.msys2.org/) with `UCRT64` environment.

Required libraries (prefix: `mingw-w64-ucrt-x86_64-`):

* clang
* cmake
* libtiff
* make
* rust
* sqlite3

Local Linux validation assumes the devcontainer and uses separate PROJ modes:

```shell
just ci-check-system
just ci-check-bundled
```

## Usage

*tyler* is a command line application.

Use `--help` to see the help menu.

```shell
tyler --help
```

Execution logs are outputted to the console.
You can control the loging level (`debug`, `info`, `error`) by setting the `RUST_LOG` environment variable.
For instance turn on the debug messages.

```shell
RUST_LOG=debug tyler ...
```

Tyler uses the [proj](https://proj.org/) library for reprojecting the input to the required CRS.
Set [PROJ_DATA](https://proj.org/usage/environmentvars.html#envvar-PROJ_DATA) if your environment requires a custom PROJ data directory.

The packaged Nix build uses the reproducible/offline path: it vendors the Dutch
`nl_nsgi_nlgeo2018.tif` and `nl_nsgi_rdtrans2018.tif` grids into `PROJ_DATA`, so
required transforms do not depend on runtime downloads.

Native PROJ networking is a separate, opt-in path for non-packaged builds. It
requires both a system PROJ build with curl support and explicit runtime
enablement in the caller. When enabled, PROJ fetches the needed remote grid
chunks on demand and stores them in its user-writable cache database; it does
not materialize full `.tif` files into `PROJ_DATA`. Reproducible or offline
deployments should continue to provide required grids locally.

### Performance

You can control the level of parallelism by setting the `RAYON_NUM_THREADS` environment variables.
By default *tyler* (rayon) will use the same number of threads as the number of CPUs available.
Note that on systems with hyperthreading enabled this equals the number of logical cores and not the physical ones.

### Calculating the extent and counting features

The input dataset is passed as the positional `input` argument, and the
included feature types (`CityObject` types) can be restricted with the
`--object-type` argument. See above for the details.

Firstly, *tyler* calculates the complete extent of the input from the bounding
box of each feature (of the specified type) that it can find in the input
dataset.
The result of this operation is reported in the logs.
The example below shows that *tyler* found `436` features of type `Building`
and `BuildingPart` in the input dataset.
The extent of the input data were calculated from these `436` features.
The computed extent is a 3D bounding box in the CRS of the input data, and it is also reported in the logs.
In the example below, the coordinates are in *RD New (EPSG: 7415)*.

```commandline
[2023-07-05T08:52:06Z INFO  tyler::parser] Found 436 features of type Some([Building, BuildingPart])
[2023-07-05T08:52:06Z INFO  tyler::parser] Ignored 0 features of type []
[2023-07-05T08:52:06Z DEBUG tyler::parser] extent_qc: BboxQc([-86804720, -26383186, -5333, -86155251, -25703867, 52882])
[2023-07-05T08:52:06Z DEBUG tyler::parser] Computed extent from features in real-world coordinates: [84995.28, 446316.814, -5.333, 85644.749, 446996.133, 52.882]
```

For large `cjindex` datasets, Tyler processes the extent in chunks and reuses
the dataset index while it works. The practical requirement is still the same:
give Tyler a dataset root that `cjindex` can resolve.

### Output formats

Tyler writes 3D Tiles by default. Select another output format with `--format`:

- `3dtiles`: writes a 3D Tiles tileset and `.glb` tile content.
- `cityjson`: writes one merged `CityJSON` document per non-empty tile.
- `cityjsonseq`: writes one `CityJSONSeq` stream per non-empty tile, with a
  `CityJSON` header and one `CityJSONFeature` item per selected source feature.

The tiled `CityJSON` and `CityJSONSeq` outputs use the same grid, quadtree,
feature filtering, LoD selection, attribute filtering and parent-attribute
inheritance path as the 3D Tiles output. They do not write `tileset.json`,
`subtrees/`, or 3D Tiles metadata.

### Exporting 3D Tiles

An example command for generating 3D Tiles.
The argument details are explained in the text below.

```shell
tyler \
    /data \
    --output /3dtiles \
    --3dtiles-implicit \
    --object-type LandUse \
    --object-type PlantCover \
    --object-type WaterBody \
    --object-type Road \
    --object-type GenericCityObject \
    --object-type Bridge \
    --object-attributes objectid:int,bronhouder:string \
    --3dtiles-metadata-class terrain \
    --grid-minz=-5 \
    --grid-maxz=300
```

### Exporting CityJSON

Use `--format cityjson` to write tiled `CityJSON` documents.

```shell
tyler \
    /data \
    --output /cityjson-tiles \
    --format cityjson \
    --object-type Building \
    --object-type BuildingPart
```

The output tiles are written under `t/<level>/<x>/<y>.city.json`.
Each tile is a merged `CityJSON` document containing the selected features for
that tile.

### Exporting CityJSONSeq

Use `--format cityjsonseq` to write tiled `CityJSONSeq` streams.

```shell
tyler \
    /data \
    --output /cityjsonseq-tiles \
    --format cityjsonseq \
    --object-type Building \
    --object-type BuildingPart
```

The output tiles are written under `t/<level>/<x>/<y>.city.jsonl`.
Each tile stream starts with a `CityJSON` header, followed by one
`CityJSONFeature` item per selected source feature. This is the same stream
shape used by `--debug-dump-data`.

#### Input data

`tyler` accepts a single `input` dataset directory with regular CityJSON files, CityJSONSeq files, or the legacy `feature-files` layout.

For example:

`tyler /some/dataset-root --output /some/output`

In case you have the legacy `feature-files` layout and the current CityJSON version (v2.0), you need to rename the CityJSON file(s) that contain the `transform` property to `metadata.json`, so the
`cityjson-index` can discover them.
For example:

```
metadata.json # incldues the transform property
feature1.city.jsonl
feature2.city.jsonl
```

In case you have the legacy `feature-files` layout and old CityJSON version, you can upgrade like this:

```
./script/concat_city_jsonl.sh metadata.city.jsonl FolderWithFeatures combined.city.jsonl
cat combined.city.jsonl | cjseq collect > combined.city.json
cjio combined.city.json upgrade save combined/combined_upg.city.json
```

#### Output

`--output`

The output is written to the directory set in `--output`.
For 3D Tiles output, it will contain a `tileset.json` file and `t/` directory with the glTF files.
In case of implicit tiling, also a `subtrees/` directory is written with the subtrees.
For `--format cityjson`, it will contain a `t/` directory with `.city.json` tile files.
For `--format cityjsonseq`, it will contain a `t/` directory with `.city.jsonl` tile files.

For `cjindex` datasets, Tyler writes a derived metadata file under `metadata/`. Per-tile CityJSONFeature streams are kept in memory by default. To inspect them, pass `--debug-dump-data`; Tyler then
writes `debug/inputs/<tile>.city.jsonl`.

#### CityObject type

CityJSON data can contain different types of CityObjects, like Building, PlantCover or Road.
It is possible to only include the selected CityObject types in the tiled output.
The CityObject types are selected with the `--object-type` argument.
This argument can be specified multiple times to select multiple object types.

For example:

`tyler … --object-type Building --object-type BuildingPart`

#### 3D Tiles metadata class

The 3D Tiles metadata specification uses the concept of classes to categorize features.
With the `--3dtiles-metadata-class` argument it is possible to set the metadata class for the features in the 3D Tiles output.
It defaults to `citymodel`.
The metadata class works in conjunction with selecting the CityObject types. Such that one declares a metadata class for a set of CityObject types.

For example:

`tyler … --3dtiles-metadata-class building --object-type Building --object-type BuildingPart`

#### Level of Detail (LoD)

CityJSON can store city objects with multiple levels of detail.
For each CityObject type, its LoD can be specified explicitly.
This is the LoD defined in the input data.
The LoD value for each CityObject type is set with the `--lod-<cityobject type>` arguments. The `<cityobject type>` is the CityJSON CityObject type, such as BuildingPart or LandUse.
The arguments are lower-case, thus “BuildingPart” becomes “building-part” and “LandUse” becomes “land-use”.
If no `--lod-<cityobject type>` is provided, Tyler selects the highest available LoD for each CityObject.

For example:

`tyler … --lod-land-use 1 --lod-building-part 1.3`

#### Attributes

Attributes on the glTF features are set with the `--object-attributes` argument.
The argument takes the attribute name and attribute value type as its value.
The attribute name and type are separated by a colon “:” and concatenated into a single string, such as “name:type”.
The possible value types are “string”, “int”, “float”, “bool”.
Multiple mappings are passed as a comma-separated list.

For example:

`tyler … --object-attributes bouwjaar:int,objectid:int,bagpandid:string,bgt_type:string`

Some datasets, such as 3DBAG or roofer exports, store the attributes on a parent `Building`
while the geometry lives on the child `BuildingPart`. In that case, enable
`--include-parent-attributes` so Tyler copies inherited parent attributes onto the
geometry-bearing child before GLB conversion.

For example:

`tyler … --object-type Building --object-type BuildingPart --include-parent-attributes`

#### Mesh normals

Use `--smooth-normals` when you want smooth vertex normals in the GLB output.

#### Colors

Colors on the glTF features are set with the `--color-<cityobject type>` arguments.
The `<cityobject type>` is the CityJSON CityObject type, such as BuildingPart or LandUse.
The arguments are lower-case, thus “BuildingPart” becomes “building-part” and “LandUse” becomes “land-use”.
The argument value is the hexadecimal rgb color value. For instance “#FF0000” is red.

For example:

`tyler … --color-building-part #FF0000`

#### Bounding volumes

*tyler* represents 3D Tiles tile bounding volumes as
[`region`](https://docs.ogc.org/cs/22-025r4/22-025r4.html#core-region)
values in EPSG:4979.

For explicit tilesets, tiles are still generated from the source CRS grid and
each tile writes its own transformed geographic region.

For implicit tilesets, content tile IDs follow the implicit geographic
subdivision of the root region: `x` subdivides longitude and `y` subdivides
latitude. Implicit tilesets use additive refinement so parent and child content
can coexist when source-leaf content appears at mixed implicit levels. The GLB
content remains encoded in the shared root ENU frame.

#### Geometric error

*tyler* currently writes full-detail content only on leaf tiles. Internal tiles
are traversal nodes, not simplified renderable representations. Their
`geometricError` is therefore computed from spatial tile size:

```text
geometricError = tile_width * geometric_error_factor
```

Leaf tiles use `geometricError: 0.0`. Tune the factor with
`--3dtiles-geometric-error-factor`; higher values make detailed content refine earlier.

For explicit tilesets, it is possible to add a tightly-fitted bounding volume to the [tile's content](https://docs.ogc.org/cs/22-025r4/22-025r4.html#core-content-bounding-volume).
You can enable this with the `--3dtiles-content-add-bv` option.

If you do want a content bounding volume, but you want it to follow the tile bounding volume exactly, you can force this with the option `--3dtiles-content-bv-from-tile`.
Usually, this happens for content that is clipped to the tile boundaries, such as terrain.
If you clip the geometry before writing the GLB with `--3dtiles-content-clip-to-tile-bounds`, this is the flag to pair with it.

## Debugging

Run *tyler* in debug mode, by setting the logging level to `debug` in the `RUST_LOG` environment variable.

```shell
RUST_LOG=debug tyler ...
```

In debug mode, or when `--debug-dump-data` is passed, *tyler* will write the `world`, `quadtree` and `tiles_failed` instances as [bincode](https://crates.io/crates/bincode) under `debug/`.
In case of a large area and lots of features (eg. an entire country and multiple millions of features), the `world.bincode` file can become a couple GB in size.
When `--debug-dump-data` is enabled, Tyler also writes intermediary per-tile CityJSONFeature streams under `debug/inputs/`.

The bincode files can be loaded by passing the directory with the bincode files to the  `--debug-load-data` parameter. When *tyler* load the instance data from the file, it will skip the instance
creation and use the loaded data instead.

The order in which *tyler* creates the instances:

1. world
2. quadtree
3. tileset
4. (implicit tileset)
5. tiles_failed
6. pruned tileset

In addition to the instance data, *tyler* can export the grid (part of the `world`), quadtree and tileset data to Tab-separated values (`.tsv`), which you can load into a GIS.
You can enable the `.tsv` export with the `--debug-dump-grid` flag.
With the `--debug-dump-grid-features` flag, also the feature centroids and their grid cell assignment will be exported.
Use `--debug-load-grid` to reload a previously exported `grid.tsv` and matching `features.tsv` before quadtree recomputation.
Only use this for a small number of features.

In debug mode, *tyler* will write the unpruned tileset too, together with the tileset that was pruned after the glTF conversion.

It is possible to only generate the tileset, without running the glTF conversion.
This can be helpful for debugging the tileset itself.
You can enable this with the `--debug-3dtiles-tileset-only` option.

## Contributing

### Profiling

Use `just profile -- --help` to see the repo-local profiling wrapper options.
Run a profile with `just profile -- --profile <case-spec>`.
The wrapper builds `tyler` in release mode with debuginfo and frame pointers,
then resolves the input and the case-specific `tyler` arguments from Geodepot.
It reads the Geodepot executable path from the repo-local `.env` file via
`GEODEPOT_BIN`, and expects each case to carry a `profile-tyler.json` data item
at the case root, for example `bvz-dh-coast-5/profile-tyler.json` for the
`bvz-dh-coast-5/bvz_dh` case spec, with the matching command-line flags.
When `--profile` is used, `--input` and `--label` are rejected because the case
spec already determines both.
Use `--runner perf`, `--runner massif`, or the default `--runner all` to choose
which profiler(s) run. This is useful when Massif is too slow for a routine
iteration.

The profiler then captures the selected `perf stat` and/or Valgrind Massif
summaries into a new directory under `docs/performance/runs/`.
Tyler's own output is stored only in a temporary staging area during the run
and is removed before the final profiling directory is published.
If the script exits with an error, the staging directory is retained as a
`<run-id>.failed/` directory under `docs/performance/runs/` so you can inspect
partial results without rerunning long profiles.

Example invocations:

```bash
just profile -- --profile bvz-dh-coast-5/bvz_dh
just profile -- --profile bvz-dh-coast-5/bvz_dh --runner perf
just profile -- --profile bvz-dh-coast-5/bvz_dh --runner massif
just profile -- --profile bvz-dh-coast-5/bvz_dh --output-root docs/performance/runs
just profile -- --input /data/cases/demo --label manual-smoke
```

Profiling Dependencies:

- [geodepot](https://github.com/3DBAG/geodepot)
- [jq](https://jqlang.github.io/jq/)
- [massif](https://valgrind.org/docs/manual/ms-manual.html)
- [perf](https://perfwiki.github.io/main/)

## Roadmap

- [x] Parallel extent computation
- [x] Parallel grid indexing
- [x] Integrate the glTF converter
- [x] Integrate cjlib
- [x] Read regular CityJSON files, not only CityJSONFeatures
- [x] Additional export formats:
    - [x] CityJSON
    - [x] CityJSONSeq
- [ ] Additional export formats:
    - [ ] Wavefront OBJ
    - [ ] GeoPackage

## Funding

- Version 0.3 (3D Tiles) was funded by the [Dutch Kadaster](https://www.kadaster.nl/).
- Version 0.4.1 was funded by the [Dutch Kadaster](https://www.kadaster.nl/).
