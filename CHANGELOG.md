# Changelog

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
