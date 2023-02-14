//! Output formats for the tiles.
//! All format

pub mod cesium3dtiles {
    //! Cesium [3D Tiles](https://github.com/CesiumGS/3d-tiles).
    //! Supported version: 1.1.
    //! Not supported: `extras`.
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::Path;

    use log::{debug, error};
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

    impl Tileset {
        /// Write the tileset to a `tileset.json` file
        pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
            let file_out = File::create(path.as_ref())?;
            serde_json::to_writer(&file_out, self)?;
            Ok(())
        }

        pub fn from_quadtree(
            quadtree: &crate::spatial_structs::QuadTree,
            world: &crate::parser::World,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
        ) -> Self {
            let crs_from = format!("EPSG:{}", world.crs.to_epsg().unwrap());
            // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
            let crs_to = "EPSG:4979";
            let transformer = Proj::new_known_crs(&crs_from, crs_to, None).unwrap();
            // y-up to z-up transform needed because we are using gltf assets, which is y-up
            // https://github.com/CesiumGS/3d-tiles/tree/main/specification#y-up-to-z-up
            // let y_up_to_z_up = Transform([
            //     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            // ]);

            let mut root = Self::generate_tiles(quadtree, world, &transformer, arg_minz, arg_maxz);
            // root.transform = Some(y_up_to_z_up);

            // Using gltf tile content
            let mut extensions: Extensions = HashMap::new();
            let e1 = Extension::ContentGtlf {
                extensions_used: None,
                extensions_required: None,
            };
            extensions.insert(ExtensionName::ContentGltf, e1);

            Self {
                asset: Default::default(),
                geometric_error: root.geometric_error * 1.5,
                root,
                properties: None,
                extensions_used: None,
                extensions_required: None,
                extensions: None,
            }
        }

        fn generate_tiles(
            quadtree: &crate::spatial_structs::QuadTree,
            world: &crate::parser::World,
            transformer: &Proj,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
        ) -> Tile {
            if !quadtree.children.is_empty() {
                if quadtree.children.len() != 4 {
                    error!("Quadtree does not have 4 children {:?}", &quadtree);
                }
                // Tile bounding volume
                let mut tile_bbox = quadtree.bbox(&world.grid);
                // Set the bounding volume height from the content height
                tile_bbox[2] = world.grid.bbox[2];
                tile_bbox[5] = world.grid.bbox[5];
                let mut bounding_volume =
                    BoundingVolume::region_from_bbox(&tile_bbox, transformer).unwrap();
                match bounding_volume {
                    BoundingVolume::Box(_) => {}
                    BoundingVolume::Region(ref mut region) => {
                        if region[5] < region[4] {
                            // This happens with the 3D Basisvoorziening data
                            debug!(
                                "Parent tile {:?} (in input CRS) bounding volume region maxz {} is less than minz {}. Replacing maxz with minz + minz * 0.01.",
                                &tile_bbox, region[5], region[4]
                            );
                            region[5] = region[4] + region[4] * 0.01;
                        }
                    }
                    BoundingVolume::Sphere(_) => {}
                }

                // The geometric error of a tile is its 'size'.
                // Since we have square tiles, we compute its size as the length of
                // its side on the x-axis.
                let dz = tile_bbox[5] - tile_bbox[2];
                if dz < 0.0 {
                    debug!("dz is negative in parent");
                }
                let mut tile_children: Vec<Tile> = Vec::new();
                for child in quadtree.children.iter() {
                    tile_children.push(Self::generate_tiles(
                        child,
                        world,
                        transformer,
                        arg_minz,
                        arg_maxz,
                    ));
                }
                Tile {
                    bounding_volume,
                    geometric_error: dz,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: None,
                    children: Some(tile_children),
                }
            } else {
                // Compute the tile content bounding box <-- the bbox of all the cells in a tile
                let mut tile_content_bbox_qc = crate::spatial_structs::BboxQc([0, 0, 0, 0, 0, 0]);
                for cellid in &quadtree.cells {
                    let cell = world.grid.cell(cellid);
                    if !cell.feature_ids.is_empty() {
                        tile_content_bbox_qc = world.features[cell.feature_ids[0]].bbox_qc.clone();
                        break;
                    }
                }
                for cellid in &quadtree.cells {
                    let cell = world.grid.cell(cellid);
                    for fi in cell.feature_ids.iter() {
                        tile_content_bbox_qc.update_with(&world.features[*fi].bbox_qc);
                    }
                }
                // If the limit-minz/maxz arguments are set, also limit the z of the
                // bounding volume. We could also just use the grid.bbox values to limit the z,
                // however at this point we don't know if that was computed from the data or set by
                // the argument. Setting the argument signals intent, so only then do we override
                // the values.
                let tile_content_bbox_rw =
                    tile_content_bbox_qc.to_bbox(&world.transform, arg_minz, arg_maxz);

                // Tile bounding volume
                let mut tile_bbox = quadtree.bbox(&world.grid);
                // Set the bounding volume height from the content height
                tile_bbox[2] = tile_content_bbox_rw[2];
                tile_bbox[5] = tile_content_bbox_rw[5];
                let mut bounding_volume =
                    BoundingVolume::region_from_bbox(&tile_bbox, transformer).unwrap();
                match bounding_volume {
                    BoundingVolume::Box(_) => {}
                    BoundingVolume::Region(ref mut region) => {
                        if region[5] < region[4] {
                            // This happens with the 3D Basisvoorziening data
                            debug!(
                                "Child tile {:?} (in input CRS) bounding volume region maxz {} is less than minz {}. Replacing maxz with minz + minz * 0.01.",
                                &tile_bbox, region[5], region[4]
                            );
                            region[5] = region[4] + region[4] * 0.01;
                        }
                    }
                    BoundingVolume::Sphere(_) => {}
                }

                // The geometric error of a tile is its 'size'.
                // Since we have square tiles, we compute its size as the length of
                // its side on the x-axis.
                let dz = tile_bbox[5] - tile_bbox[2];
                if dz < 0.0 {
                    debug!("dz is negative in child");
                }
                let content_bounding_voume =
                    BoundingVolume::region_from_bbox(&tile_content_bbox_rw, transformer).unwrap();
                match content_bounding_voume {
                    BoundingVolume::Box(_) => {}
                    BoundingVolume::Region(ref region) => {
                        if region[5] < region[4] {
                            debug!(
                                "content bounding volume region maxz {} is less than minz {}",
                                region[5], region[4]
                            )
                        }
                    }
                    BoundingVolume::Sphere(_) => {}
                }

                // FIXME: this is a hack to replace the tile bounding volume with the content bounding volume if the content is larger than the tile
                match bounding_volume {
                    BoundingVolume::Box(_) => {}
                    BoundingVolume::Region(ref mut region) => match content_bounding_voume {
                        BoundingVolume::Box(_) => {}
                        BoundingVolume::Region(ref content_region) => {
                            if content_region[0] < region[0] {
                                region[0] = content_region[0];
                            }
                            if content_region[1] < region[1] {
                                region[1] = content_region[1];
                            }
                            if content_region[4] < region[4] {
                                region[4] = content_region[4];
                            }
                            if content_region[2] > region[2] {
                                region[2] = content_region[2];
                            }
                            if content_region[3] > region[3] {
                                region[3] = content_region[3];
                            }
                            if content_region[5] > region[5] {
                                region[5] = content_region[5];
                            }
                            debug!("Updated child tile {:?} (in input CRS) bounding region from content region, because the content was larger", &tile_bbox);
                        }
                        BoundingVolume::Sphere(_) => {}
                    },
                    BoundingVolume::Sphere(_) => {}
                }

                Tile {
                    bounding_volume,
                    geometric_error: 0.0,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: Some(Content {
                        bounding_volume: Some(content_bounding_voume),
                        uri: format!("tiles/{}.glb", quadtree.id()),
                    }),
                    children: None,
                }
            }
        }

        #[allow(dead_code)]
        pub fn from_grid(
            grid: &crate::spatial_structs::SquareGrid,
            citymodel: &crate::parser::CityJSONMetadata,
            feature_set: &crate::parser::FeatureSet,
        ) -> Self {
            let crs_from = format!(
                "EPSG:{}",
                citymodel.metadata.reference_system.to_epsg().unwrap()
            );
            // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
            let crs_to = "EPSG:4979";
            let transformer = Proj::new_known_crs(&crs_from, crs_to, None).unwrap();

            let mut root_children: Vec<Tile> = Vec::with_capacity(grid.length * grid.length);
            for (cellid, cell) in grid {
                if cell.feature_ids.is_empty() {
                    // Empty cell, don't create tiles for it
                    debug!("cell {} is empty", cellid);
                    continue;
                }

                let mut content_bbox_qc = feature_set[cell.feature_ids[0]].bbox_qc.clone();
                for fi in cell.feature_ids.iter() {
                    content_bbox_qc.update_with(&feature_set[*fi].bbox_qc);
                }
                let content_bbox_rw = content_bbox_qc.to_bbox(&citymodel.transform, None, None);
                let content_bounding_voume =
                    BoundingVolume::region_from_bbox(&content_bbox_rw, &transformer).unwrap();

                let mut cell_bbox = grid.cell_bbox(&cellid);
                // Set the bounding volume height from the content height
                cell_bbox[2] = content_bbox_rw[2];
                cell_bbox[5] = content_bbox_rw[5];
                let bounding_volume =
                    BoundingVolume::region_from_bbox(&cell_bbox, &transformer).unwrap();

                // We are adding a child for each LoD.
                // TODO: but we are cheating here now, because we know that the input data has 3 LoDs...

                // The geometric error of a tile is its 'size'.
                // Since we have square tiles, we compute its size as the length of
                // its side on the x-axis.
                let dz = cell_bbox[5] - cell_bbox[2];

                // this is a leaf node, so the geometric_error is 0
                // LoD2.2
                let tile_lod22 = Tile {
                    bounding_volume,
                    geometric_error: 0.0,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: Some(Content {
                        bounding_volume: Some(content_bounding_voume),
                        uri: format!("tiles/{}-0-0.glb", cellid),
                    }),
                    children: None,
                };

                // LoD 1.3
                let tile_lod13 = Tile {
                    bounding_volume,
                    geometric_error: dz * 0.1,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: Some(Content {
                        bounding_volume: Some(content_bounding_voume),
                        uri: format!("tiles/{}-0.glb", cellid),
                    }),
                    children: Some(vec![tile_lod22]),
                };

                // LoD 1.2
                root_children.push(Tile {
                    bounding_volume,
                    geometric_error: dz * 0.3,
                    // geometric_error: 10.0,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: Some(Content {
                        bounding_volume: Some(content_bounding_voume),
                        uri: format!("tiles/{}.glb", cellid),
                    }),
                    children: Some(vec![tile_lod13]),
                });
            }

            let root_volume = BoundingVolume::region_from_bbox(&grid.bbox, &transformer).unwrap();
            debug!("root bbox: {:?}", &grid.bbox);
            debug!("root boundingVolume: {:?}", &root_volume);
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
                extensions_used: None,
                extensions_required: None,
                extensions: None,
            }
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

    #[allow(dead_code)]
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
    #[allow(dead_code)]
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

    /// Compute the boundingVolume.box from a 'regular' bounding box.
    ///
    /// This function does not reproject the bounding box coordinates.
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
    impl From<&crate::spatial_structs::Bbox> for BoundingVolume {
        fn from(bbox: &crate::spatial_structs::Bbox) -> Self {
            let dx = bbox[3] - bbox[0];
            let dy = bbox[4] - bbox[1];
            let dz = bbox[5] - bbox[2];
            debug!("bounding volume box dx {} dy {} dz {}", dx, dy, dz);
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
        /// Compute the boundingVolume.box from a 'regular' bounding box.
        ///
        /// This function does reproject the bounding box coordinates.
        ///
        /// The CRS transformation `transformer` must have `EPSG:4978` as target CRS in
        /// order to get a correct `boundingVolume.box`. The `transformer` is not initialized in
        /// this function, because it is expected that this function is called in a loop, and thus
        /// it is more optimal to init the transformation outside of the loop.
        ///
        /// ## Note on CRS transformation
        /// We assume that the transformation is from a projected, Cartesian system to
        /// [ECEF](https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system).
        /// The coordinates of a polygon on the northern-hemisphere in ECEF look like in the image
        /// below. Note that the minimum point is the top-left corner of the polygon.
        ///
        /// ![ecef polygon](../docs/img/ecef.jpg)
        ///
        /// In a projected Cartesian CRS we have the minimum point of a polygon in the lower-left
        /// corner (at least in the Netherlands...).
        ///
        /// ```shell
        /// (0,10) (10,10)   (0,0)  (0,10)
        ///      +--+           +------+
        ///      |  |     --->  | ECEF |
        ///      +--+           +------+
        /// (0,0)  (10,0)    (10,0) (10,10)
        /// ```
        ///
        /// Then the input bbox (in 3D) of `[0,0,0,10,10,10]` becomes `[10,0,0,0,10,10]` after the
        /// coordinate transformation.
        /// Therefore, we have to swap the values at the indices `0` and `3` in order to get a
        /// correct `[minX, minY, minZ, maxX, maxY, maxZ]` bbox after the transformation.
        ///
        #[allow(dead_code)]
        fn box_from_bbox(
            bbox: &crate::spatial_structs::Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            let min_coord = transformer.convert((bbox[0], bbox[1], bbox[2]))?;
            let max_coord = transformer.convert((bbox[3], bbox[4], bbox[5]))?;
            debug!(
                "bounding volume reprojected box dz {}",
                max_coord.2 - min_coord.2
            );
            Ok(BoundingVolume::from(&[
                max_coord.0,
                min_coord.1,
                min_coord.2,
                min_coord.0,
                max_coord.1,
                max_coord.2,
            ]))
        }

        fn region_from_bbox(
            bbox: &crate::spatial_structs::Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            let (west, south, minh) = transformer.convert((bbox[0], bbox[1], bbox[2]))?;
            let (east, north, maxh) = transformer.convert((bbox[3], bbox[4], bbox[5]))?;
            Ok(BoundingVolume::Region([
                west.to_radians(),
                south.to_radians(),
                east.to_radians(),
                north.to_radians(),
                minh,
                maxh,
            ]))
        }
    }

    /// [Tile.refine](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tilerefine).
    #[allow(dead_code)]
    #[derive(Serialize, Debug)]
    #[serde(rename_all = "UPPERCASE")]
    enum Refinement {
        Add,
        Replace,
    }

    /// [Tile.transform](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tiletransform)
    #[derive(Serialize, Debug, Copy, Clone)]
    struct Transform([f64; 16]);

    impl Default for Transform {
        #[rustfmt::skip]
        fn default() -> Self {
            Self([
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0,
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
            let bbox: crate::spatial_structs::Bbox =
                [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
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
