//! Output formats for the tiles.
//! All format

pub mod cesium3dtiles {
    //! Cesium [3D Tiles](https://github.com/CesiumGS/3d-tiles).
    //! Supported version: 1.1.
    //! Not supported: `extras`.
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::Path;

    use log::debug;
    use serde::Serialize;

    use crate::proj::Proj;

    /// [Tileset](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tileset).
    ///
    /// Not supported: `extras`.
    #[derive(Serialize, Default, Debug)]
    #[serde(rename_all = "camelCase")]
    pub struct Tileset {
        asset: Asset,
        geometric_error: GeometricError,
        root: Tile,
        #[serde(skip_serializing_if = "Option::is_none")]
        properties: Option<Properties>,
        #[serde(skip_serializing_if = "Option::is_none")]
        extensions_used: Option<Vec<ExtensionName>>,
        #[serde(skip_serializing_if = "Option::is_none")]
        extensions_required: Option<Vec<ExtensionName>>,
        #[serde(skip_serializing_if = "Option::is_none")]
        extensions: Option<Extensions>,
    }

    impl From<&crate::spatial_structs::SquareGrid> for Tileset {
        fn from(grid: &crate::spatial_structs::SquareGrid) -> Self {
            let crs_from = format!("EPSG:{}", grid.epsg);
            // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
            let crs_to = "EPSG:4978";
            let transformer = Proj::new_known_crs(&crs_from, &crs_to, None).unwrap();

            let mut root_children: Vec<Tile> = Vec::with_capacity(grid.length ^ 2);
            for (cellid, _cell) in grid {
                let cell_bbox = grid.cell_bbox(&cellid);
                debug!("{}-{} bbox: {:?}", cellid[0], cellid[1], &cell_bbox);
                let bounding_volume =
                    BoundingVolume::from_bbox_reproject(&cell_bbox, &transformer).unwrap();
                debug!(
                    "{}-{} boundingVolume: {:?}",
                    cellid[0], cellid[1], &bounding_volume
                );
                // The geometric error of a tile is its 'size'.
                // Since we have square tiles, we compute its size as the length of
                // its side on the x-axis.
                // let geometric_error = cell_bbox[3] - cell_bbox[0];
                // but because now this is a leaf, the geometric error has to be close or equal to 0
                let geometric_error = 0.0;
                let content = Content {
                    bounding_volume: None,
                    uri: format!("{}-{}.glb", cellid[0], cellid[1]),
                };

                // TODO: For each LoD add a child

                root_children.push(Tile {
                    bounding_volume,
                    geometric_error,
                    viewer_request_volume: None,
                    refine: None,
                    transform: None,
                    content: Some(content),
                    children: None,
                })
            }

            let root_volume =
                BoundingVolume::from_bbox_reproject(&grid.bbox, &transformer).unwrap();
            let root_geometric_error = grid.bbox[3] - grid.bbox[0];
            let root = Tile {
                bounding_volume: root_volume,
                geometric_error: root_geometric_error,
                viewer_request_volume: None,
                refine: Some(Refinement::Replace),
                transform: None,
                content: None,
                children: Some(root_children),
            };

            // Using gltf tile content
            let mut extensions: Extensions = HashMap::new();
            let e1 = Extension::ContentGtlf {
                extensions_used: None,
                extensions_required: None,
            };
            extensions.insert(ExtensionName::ContentGltf, e1);

            Self {
                asset: Default::default(),
                geometric_error: root_geometric_error * 1.5,
                root,
                properties: None,
                extensions_used: Some(vec![ExtensionName::ContentGltf]),
                extensions_required: Some(vec![ExtensionName::ContentGltf]),
                extensions: Some(extensions),
            }
        }
    }

    impl Tileset {
        /// Write the tileset to a `tileset.json` file
        pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
            let file_out = File::create(path.as_ref())?;
            serde_json::to_writer(&file_out, self)?;
            Ok(())
        }
    }

    /// [Asset](https://github.com/CesiumGS/3d-tiles/tree/main/specification#asset).
    ///
    /// Not supported: `extensions, extras`.
    #[derive(Serialize, Debug)]
    #[serde(rename_all = "camelCase")]
    struct Asset {
        version: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        tileset_version: Option<String>,
    }

    impl Default for Asset {
        fn default() -> Self {
            Self {
                version: String::from("1.1"),
                tileset_version: None,
            }
        }
    }

    /// [geometricError](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tilesetgeometricerror-white_check_mark).
    /// Must be `>=0`.
    type GeometricError = f64;

    /// [Properties](https://github.com/CesiumGS/3d-tiles/tree/main/specification#properties).
    ///
    /// Not supported: `extras`.
    #[derive(Serialize, Default, Debug)]
    #[serde(rename_all = "camelCase")]
    struct Properties {
        maximum: f64,
        minimum: f64,
    }

    type Extensions = HashMap<ExtensionName, Extension>;

    #[derive(Serialize, Default, Debug)]
    #[serde(untagged)]
    enum Extension {
        #[default]
        None,
        /// [3DTILES_content_gltf](https://github.com/CesiumGS/3d-tiles/tree/main/extensions/3DTILES_content_gltf).
        /// Allows to use glTF tile content instead of Batched 3DModel (.b3dm).
        #[serde(rename = "3DTILES_content_gltf")]
        ContentGtlf {
            #[serde(rename = "extensionsUsed", skip_serializing_if = "Option::is_none")]
            extensions_used: Option<Vec<ExtensionName>>,
            #[serde(rename = "extensionsRequired", skip_serializing_if = "Option::is_none")]
            extensions_required: Option<Vec<ExtensionName>>,
        },
    }

    #[derive(Serialize, Default, Debug, PartialEq, Eq, Hash)]
    enum ExtensionName {
        #[default]
        None,
        #[serde(rename = "3DTILES_content_gltf")]
        ContentGltf,
        #[serde(rename = "EXT_mesh_features")]
        MeshFeatures,
        #[serde(rename = "EXT_structural_metadata")]
        StructuralMetadata,
        #[serde(rename = "3DTILES_implicit_tiling")]
        ImplicitTiling,
    }

    /// [Tile](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tile).
    #[derive(Serialize, Default, Debug)]
    #[serde(rename_all = "camelCase")]
    struct Tile {
        bounding_volume: BoundingVolume,
        geometric_error: GeometricError,
        #[serde(skip_serializing_if = "Option::is_none")]
        viewer_request_volume: Option<BoundingVolume>,
        #[serde(skip_serializing_if = "Option::is_none")]
        refine: Option<Refinement>,
        #[serde(skip_serializing_if = "Option::is_none")]
        transform: Option<Transform>,
        #[serde(skip_serializing_if = "Option::is_none")]
        content: Option<Content>,
        #[serde(skip_serializing_if = "Option::is_none")]
        children: Option<Vec<Tile>>,
    }

    /// [boundingVolume](https://github.com/CesiumGS/3d-tiles/tree/main/specification#bounding-volume).
    #[derive(Serialize, Debug, Copy, Clone)]
    #[serde(rename_all = "lowercase")]
    enum BoundingVolume {
        Box([f64; 12]),
        Region([f64; 6]),
        Sphere([f64; 4]),
    }

    impl Default for BoundingVolume {
        fn default() -> Self {
            Self::Region([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        }
    }

    /// Compute the oriented-bounding box for 3D Tiles, from a 'regular' bounding box.
    ///
    /// The output is an array of 12 numbers that define an oriented bounding box in a
    /// right-handed 3-axis (x, y, z) Cartesian coordinate system where the z-axis is up.
    /// The first three elements define the x, y, and z values for the center of the box.
    /// The next three elements (with indices 3, 4, and 5) define the x-axis direction
    /// and half-length.
    /// The next three elements (indices 6, 7, and 8) define the y-axis direction and
    /// half-length.
    /// The last three elements (indices 9, 10, and 11) define the z-axis direction and
    /// half-length.
    impl From<&crate::Bbox> for BoundingVolume {
        fn from(bbox: &crate::Bbox) -> Self {
            let dx = bbox[0] - bbox[3];
            let dy = bbox[1] - bbox[4];
            let dz = bbox[2] - bbox[5];
            let center: [f64; 3] = [bbox[0] + dx * 0.5, bbox[1] + dy * 0.5, bbox[2] + dz * 0.5];
            // The x-direction and half-length
            let x_axis_dir: [f64; 3] = [dx * 0.5, 0.0, 0.0];
            // The y-direction and half-length
            let y_axis_dir: [f64; 3] = [0.0, dy * 0.5, 0.0];
            // The z-direction and half-length
            let z_axis_dir: [f64; 3] = [0.0, 0.0, dz * 0.5];
            // an array cannot be built directly from an iterator, so we need to
            // "try_(to convert the vector)_into" an array
            let bounding_volume_array: [f64; 12] = [center, x_axis_dir, y_axis_dir, z_axis_dir]
                .into_iter()
                .flatten()
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap();
            Self::Box(bounding_volume_array)
        }
    }

    impl BoundingVolume {
        fn from_bbox_reproject(
            bbox: &crate::Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            let min_coord = transformer.convert((bbox[0], bbox[1], bbox[2]))?;
            let max_coord = transformer.convert((bbox[3], bbox[4], bbox[5]))?;
            Ok(BoundingVolume::from(&[
                min_coord.0,
                min_coord.1,
                min_coord.2,
                max_coord.0,
                max_coord.1,
                max_coord.2,
            ]))
        }
    }

    /// [Tile.refine](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tilerefine).
    #[derive(Serialize, Debug)]
    #[serde(rename_all = "UPPERCASE")]
    enum Refinement {
        Add,
        Replace,
    }

    /// [Tile.transform](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tiletransform)
    #[derive(Serialize, Debug)]
    struct Transform([f64; 16]);

    impl Default for Transform {
        fn default() -> Self {
            Self([
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ])
        }
    }

    /// [Tile.content](https://github.com/CesiumGS/3d-tiles/tree/main/specification#content).
    #[derive(Serialize, Default, Debug)]
    #[serde(rename_all = "camelCase")]
    struct Content {
        #[serde(skip_serializing_if = "Option::is_none")]
        bounding_volume: Option<BoundingVolume>,
        uri: String,
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use serde_json::to_string_pretty;

        #[test]
        fn test_refinement() {
            let r = Refinement::Replace;
            let j = serde_json::to_string(&r).unwrap();
            assert_eq!(j, r#""REPLACE""#.to_string());
        }

        #[test]
        fn test_boundingvolume_from_bbox() {
            let bbox: crate::Bbox = [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
            let bounding_volume = BoundingVolume::from(&bbox);
            println!("{:?}", bounding_volume);
        }

        /// Verify that we can serialize the 3DTILES_content_gltf as an empty object when it is just
        /// declared as a 'null' extension in the tileset.
        /// Should print:
        /// ```json
        ///   "extensions": {
        ///     "3DTILES_content_gltf": {}
        ///   },
        /// ```
        #[test]
        fn test_3d_tiles_empty_extension() {
            let mut extensions: Extensions = HashMap::new();
            let e1 = Extension::ContentGtlf {
                extensions_used: None,
                extensions_required: None,
            };
            extensions.insert(ExtensionName::ContentGltf, e1);

            let t = Tileset {
                asset: Default::default(),
                geometric_error: 0.0,
                properties: None,
                extensions_used: Some(vec![ExtensionName::ContentGltf]),
                extensions_required: Some(vec![ExtensionName::ContentGltf]),
                extensions: Some(extensions),
                root: Default::default(),
            };
            println!("{}", to_string_pretty(&t).unwrap());
        }

        /// Verify that we can serialize the 3DTILES_content_gltf extension that uses other extensions.
        #[test]
        fn test_3d_tiles_with_extension() {
            let mut extensions: Extensions = HashMap::new();
            let e1 = Extension::ContentGtlf {
                extensions_used: Some(vec![
                    ExtensionName::MeshFeatures,
                    ExtensionName::StructuralMetadata,
                ]),
                extensions_required: None,
            };
            extensions.insert(ExtensionName::ContentGltf, e1);

            let t = Tileset {
                asset: Default::default(),
                geometric_error: 0.0,
                properties: None,
                extensions_used: Some(vec![ExtensionName::ContentGltf]),
                extensions_required: Some(vec![ExtensionName::ContentGltf]),
                extensions: Some(extensions),
                root: Default::default(),
            };
            println!("{}", to_string_pretty(&t).unwrap());
        }
    }
}
