# tyler

<p align="center">
  <img src="https://github.com/3DGI/tyler/blob/master/tyler.png" />
</p>

*tyler* creates tiles from 3D city objects.

As input, *tyler* takes [CityJSON Features](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature), where each feature is stored in a separate file.

As output, *tyler* can create:

- [3D Tiles v1.1](https://docs.ogc.org/cs/22-025r4/22-025r4.html)

Details of the 3D Tiles output:

- The tileset content if binary glTF (.glb).
- The glTF assets contain feature metadata (per CityObject), using the [EXT_mesh_features](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_mesh_features) and [EXT_structural_metadata](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata) extensions.
- The features are colored to default values, and the colors can by set per CityObject type.
- The glTF files are compressed, using the [KHR_mesh_quantization](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Khronos/KHR_mesh_quantization) and [EXT_meshopt_compression](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Vendor/EXT_meshopt_compression) extensions.
- Implicit tiling is supported (optional).

Additional information about the internals of *tyler* you will find in the [design document](https://github.com/3DGI/tyler/blob/master/docs/design_document.md).

## Installation

For the time being, *tyler* depends on the [geoflow-bundle](https://github.com/geoflow3d/geoflow-bundle) for converting CityJSONFeatures to glTF.
Unless you want to install the *geoflow-bundle* yourself, we strongly recommend to use the provided docker image for running *tyler*, because it contains the *geoflow-bundle*.

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

### Exporting 3D Tiles

An example command for generating 3D Tiles. 
The agrument details are explained in the text below.

```shell
tyler \
    --metadata metadata.city.json \
    --features features/ \
    --output /3dtiles \
    --3dtiles-implicit \
    --object-type LandUse \
    --object-type PlantCover \
    --object-type WaterBody \
    --object-type Road \
    --object-type GenericCityObject \
    --object-type Bridge \
    --object-attribute objectid:int,bronhouder:string \
    --3dtiles-metadata-class terrain \
    --grid-minz=-5 \
    --grid-maxz=300
```

#### Input data

1. A main `.city.json` file, containing at least the [CRS](https://www.cityjson.org/specs/1.1.3/#referencesystem-crs) and [transform](https://www.cityjson.org/specs/1.1.3/#transform-object) objects.
2. A directory (or directory tree) of `.city.jsonl` files, each containing one [CityJSON Feature](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature), including all its children City Objects.

`--metadata`

A main `.city.json` file, containing at least the `CRS` and `transform` objects, set by the argument.

`--features`

A directory (or directory tree) of `.city.jsonl` files, each containing one CityJSON Feature, including all its children City Objects.

For example:

`tyler --metadata metadata.city.json --features /some/directory/`

#### Output

`--output`

The output is written to the directory set in `--output`. 
For 3D Tiles output, it will contain a `tileset.json` file and `tiles/` directory with the glTF files. 
In case of implicit tiling, also a `subtrees/` directory is written with the subtrees.

During the operation of Tyler, also an `input/` directory is created with text files, but this directory is removed with all its content after Tyler finished processing the tiles (except when debug mode is enabled).

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

#### Colors

Colors on the glTF features are set with the `--color-<cityobject type>` arguments. 
The `<cityobject type>` is the CityJSON CityObject type, such as BuildingPart or LandUse. 
The arguments are lower-case, thus “BuildingPart” becomes “building-part” and “LandUse” becomes “land-use”.
The argument value is the hexadecimal rgb color value. For instance “#FF0000” is red.

For example:

`tyler … --color-building-part #FF0000`
