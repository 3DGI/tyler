# Changelog

## tyler 0.3.7 (unreleased)

### Fixed
- Geometric error calculation for the parents of the deepest leaves, where the leaf has an edge length of the grid cellsize.
- Replace maxz in leaf content if minz is greater than maxz.
- Tile bounding volume when tile is empty.
- Fixed the hardcoded timeout.

### Added
- Set content bounding volume from tile bounding volume with the `--3dtiles-content-bv-from-tile` option. Required when the tile content is clipped to the tile's extent, for example for terrain.
- Write quadtree content bounding box to .tsv with `--grid-export`.
- Log all arguments in debug.

### Changed
- Debug data, incl. `--grid-export`, is written to the `debug` directory within the `--output` directory.
- Remove logging from geof. Speeds up the conversion and fixes the extreme memory consumption when geof emits a large amount of messages.
- Use BufWriter for writing the input paths.
- Implement parallel computation for the extent, where the direct subdirectories of `--features` are visited in parallel (but their contents are processed sequentially).
- Floats in the 3D Tiles tileset.json are written with 6 decimals

## tyler 0.3.6 (2023-07-17)

### Changed
- The `--grid-export` switch does not export the feature centroids anymore. Use the `--grid-export-features` if you want to export the feature centroids together with the grid cells.
- Write both pruned and unpruned tilesets. Unpruned tiles are only written in debug mode.
- Reduced the logging in debug mode.

### Fixed
- Invalid subtree for implicit tiling, in case of very large areas (eg. the Netherlands).

### Added
- Write the `world`, `quadtree` and `tiles_failed` instances to bincode when running in debug mode. The instances can be loaded for debugging with the `--debug-load-data`, in which case tyler will load the instance data instead of generating it.
- Describe how to generate debug data.

## tyler 0.3.5 (2023-06-28)

### Added
- Option to only generate and write the 3D Tiles tileset, without running the glTF export. Enable with `--3dtiles-tileset-only`.
- Timeout `--timeout` in seconds for the converter subprocesses. If speficied, tyler will kill the subprocess after the provided seconds, otherwise it will wait for the process to finish.

### Fixed
- geoflow version reporting

## tyler 0.3.4 (2023-06-22)

### Fixed
- geoflow-bundle version in the docker image

## tyler 0.3.3 (2023-06-19)

### Fixed
- Infinite loop in the glb conversion in geoflow-bundle, in a rare case during the mesh simplification.
- Infinite loop in case there is not a single CityJSONFeature file in the directory tree.

### Changed
- tyler and geoflow versions are reported as info

## tyler 0.3.2 (2023-05-05)

### Fixed
- *proj* related fixes for running Tyler under Windows

## tyler 0.3.1 (2023-04-06)

### Changed

- Add the *proj* crate as a submodule, because the proj-sys build script need to be changed so that proj-sys can be built with MYSYS2 on Windows (see https://github.com/georust/proj/pull/156).
- Warning instead of error when the gltf export in the subprocess fails (fixes #36).
- The geoflow flowchart directory can be set at runtime with the `TYLER_RESOURCES_DIR` environment variable. By default, the directory is set to the `CARGO_MANIFEST_DIR` at compile time.
- Features are assigned to tiles based on their bounding box, instead of only their vertices (fixes #28).
- Update the geoflow docker image.
- Improve documentation.

## tyler 0.3.0 (2023-03-17)

First public release of Tyler.
With this version Tyler can generate [3D Tiles v1.1](https://docs.ogc.org/cs/22-025r4/22-025r4.html) from [CityJSONFeatures](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature) that are stored individually in files.

Details of the 3D Tiles output:

- The tileset content if binary glTF (.glb).
- The glTF assets contain feature metadata (per CityObject), using the [EXT_mesh_features](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_mesh_features) and [EXT_structural_metadata](https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata) extensions.
- The features are colored to default values, and the colors can by set per CityObject type.
- The glTF files are compressed, using the [KHR_mesh_quantization](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Khronos/KHR_mesh_quantization) and [EXT_meshopt_compression](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Vendor/EXT_meshopt_compression) extensions.
- Implicit tiling is supported (optional).

The current version depends on the [geoflow-bundle](https://github.com/geoflow3d/geoflow-bundle) for converting CityJSONFeatures to glTF.
Therefore, we strongly recommend to use the provided docker image for running Tyler.
