# tyler

<p align="center">
  <img src="https://github.com/3DGI/tyler/blob/master/tyler.png" />
</p>

*tyler* creates tiles from 3D city objects.

As input, *tyler* takes a dataset root. That root can either be a legacy
directory containing `metadata.city.json` plus
[CityJSON Features](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature),
or a `cjindex` dataset root.

As output, *tyler* can create:

- [3D Tiles v1.1](https://docs.ogc.org/cs/22-025r4/22-025r4.html)

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

Pull the docker image with `docker pull 3dgi/tyler:<version>`, e.g. `docker pull 3dgi/tyler:0.3.12`.

### Using the pre-compiled binaries on windows

1. Download the latest Tyler binary package for windows from the [Tyler release page](https://github.com/3DGI/tyler/releases)
1. Unzip the Tyler binary package to a folder, for example `C:\software\tyler`
1. You can now run Tyler using the `run_tyler_example.bat` file inside this directory by double clicking on it. You can also copy and open this file in a text editor to change the parameters (eg.
   input and output data directories) used for running.

For testing purposes you download [this sample data](https://data.3dgi.xyz/3dtiles-test-data/download/3D-basisvoorziening-2021-30dz1_01.zip). Create a `data` folder in same the folder as the `.bat`
file mentioned above and unzip the contents there.

Contents of the `run_tyler_example.bat` file :

```
set RUST_LOG=debug
set PROJ_DATA=%~dp0\share\proj

%~dp0\bin\tyler.exe ^
%~dp0\data ^
--output %~dp0\data-out\3dtiles-terrain ^
--3dtiles-implicit ^
--object-type LandUse ^
--object-type PlantCover ^
--object-type WaterBody ^
--object-type Road ^
--object-type GenericCityObject ^
--object-type Bridge ^
--object-attribute objectid:int ^
--object-attribute bronhouder:string ^
--object-attribute bgt_fysiekvoorkomen:string ^
--object-attribute bgt_type:string ^
--3dtiles-metadata-class terrain ^
--grid-minz=-15 ^
--grid-maxz=400 >> log.txt 2>&1
```

### Compiling from source

*tyler* is written in Rust and you need the [Rust toolchain](https://www.rust-lang.org/learn/get-started) to compile it.

After downloading the source code from GitHub, navigate into the tyler directory and you can install *tyler* with *cargo*.

```shell
cargo install .
```

#### On Windows

Use [MSYS2](https://www.msys2.org/) with `UCRT64` environment.

Required libraries (prefix: `mingw-w64-ucrt-x86_64-`):

* clang
* cmake
* libtiff
* make
* rust
* sqlite3

CI also installs `libtiff` on Linux because `proj-sys` falls back to building bundled PROJ with TIFF support when a suitable system `libproj` is not available.

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

### Performance

For large input, like multiple millions of features you need have an SSD.
Running on a HDD is not feasible for large areas.

There are three resource intensive steps, 1) computing the extent of the input, 2) indexing the input with the grid, 3) converting the tiles.
Each of the three steps are executed concurrently, with the help of the [rayon library](https://crates.io/crates/rayon).

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

For legacy feature-file datasets, the extent calculation is done in parallel
for each input subdirectory, if there are any. The contents of each
subdirectory are processed sequentially. Individual files directly under the
dataset root are processed sequentially after the subdirectories. Therefore,
for optimal performance, organize legacy feature files into subdirectories.

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
    --object-attribute objectid:int \
    --object-attribute bronhouder:string \
    --3dtiles-metadata-class terrain \
    --grid-minz=-5 \
    --grid-maxz=300
```

#### Input data

`tyler` accepts a single `input` dataset root:

- Legacy feature-file dataset:
  the directory must contain `metadata.city.json` and one or more
  `.city.jsonl` feature files under that root.
- `cjindex` dataset:
  the directory must resolve as a `cjindex` dataset root in `ndjson`,
  `cityjson`, or `feature-files` layout.

For example:

`tyler /some/dataset-root --output /some/output`

#### Output

`--output`

The output is written to the directory set in `--output`.
For 3D Tiles output, it will contain a `tileset.json` file and `t/` directory with the glTF files.
In case of implicit tiling, also a `subtrees/` directory is written with the subtrees.

For `cjindex` datasets, Tyler writes a derived metadata file under `metadata/`. Per-tile CityJSONFeature streams are kept in memory by default. To inspect them, pass `--debug-tile-inputs`; Tyler then
writes `inputs/<tile>.city.jsonl` and keeps the `inputs/` directory.

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
The metadata class works in conjunction with selecting the CityObject types. Such that one declares a metadata class for a set of CityObject types.

For example:

`tyler … --3dtiles-metadata-class building --object-type Building --object-type BuildingPart`

#### Level of Detail (LoD)

CityJSON can store city objects with multiple levels of detail.
For each CityObject type, its LoD needs to be specified as well.
This is the LoD defined in the input data.
The LoD value for each CityObject type is set with the `--lod-<cityobject type>` arguments. The `<cityobject type>` is the CityJSON CityObject type, such as BuildingPart or LandUse.
The arguments are lower-case, thus “BuildingPart” becomes “building-part” and “LandUse” becomes “land-use”.
If the value of `--lod-<cityobject type>` is an empty string (this is the default), then Tyler will select the highest available LoD for the city object.

For example:

`tyler … --lod-land-use 1 --lod-building-part 1.3`

#### Attributes

Attributes on the glTF features are set with the `--object-attribute` argument.
The argument takes the attribute name and attribute value type as its value.
The attribute name and type are separated by a colon “:” and concatenated into a single string, such as “name:type”.
The possible value types are “string”, “int”, “float”, “bool”.
The `--object-attribute` argument can be specified multiple times to include multiple attributes.

For example:

`tyler … --object-attribute bouwjaar:int --object-attribute objectid:int --object-attribute bagpandid:string --object-attribute bgt_type:string`

Some datasets, such as 3DBAG or roofer exports, store the attributes on a parent `Building`
while the geometry lives on the child `BuildingPart`. In that case, enable
`--include-parent-attributes` so Tyler copies inherited parent attributes onto the
geometry-bearing child before GLB conversion.

For example:

`tyler … --object-type Building --object-type BuildingPart --include-parent-attributes`

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
`--geometric-error-factor`; higher values make detailed content refine earlier.

For explicit tilesets, it is possible to add a tightly-fitted bounding volume to the [tile's content](https://docs.ogc.org/cs/22-025r4/22-025r4.html#core-content-bounding-volume).
You can enable this with the `--3dtiles-content-add-bv` option.

If you do want a content bounding volume, but you want it to follow the tile bounding volume exactly, you can force this with the option `--3dtiles-content-bv-from-tile`.
Usually, this happens for content that is clipped to the tile boundaries, such as terrain.

## Debugging

Run *tyler* in debug mode, by setting the logging level to `debug` in the `RUST_LOG` environment variable.

```shell
RUST_LOG=debug tyler ...
```

In debug mode, *tyler* will write the `world`, `quadtree` and `tiles_failed` instances to [bincode](https://crates.io/crates/bincode) to the working directory.
In case of a large area and lots of features (eg. an entire country and multiple millions of features), the `world.bincode` file can become a couple GB in size.

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
You can enable the `.tsv` export with the `--grid-export` flag.
With the `--grid-export-features` flag, also the feature feature centorids and their grid cell assignment will be exported.
Only use this for small amount of features.

In debug mode, *tyler* will write the unpruned tileset too, together with the tileset that was pruned after the glTF conversion.

It is possible to only generate the tileset, without running the glTF conversion.
This can be helpful for debugging the tileset itself.
You can enable this with the `--3dtiles-tileset-only` option.

## Roadmap

- [x] Parallel extent computation
- [x] Parallel grid indexing
- [x] Integrate the glTF converter
- [x] Integrate cjlib
- [x] Read regular CityJSON files, not only CityJSONFeatures
- [ ] Additional export formats:
    - [ ] CityJSON
    - [ ] Wavefront OBJ
    - [ ] GeoPackage

## Funding

Version 0.3 (3D Tiles) was funded by the [Dutch Kadaster](https://www.kadaster.nl/).
