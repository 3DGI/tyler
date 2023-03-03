//! Output formats for the tiles.
//! All format

pub mod cesium3dtiles {
    //! Cesium [3D Tiles](https://github.com/CesiumGS/3d-tiles).
    //! Supported version: 1.1.
    //! Not supported: `extras`.
    use std::collections::HashMap;
    use std::collections::VecDeque;
    use std::fmt::{Display, Formatter};
    use std::fs::File;
    use std::io::Write;
    use std::path::Path;

    use bitvec::prelude as bv;
    use log::{debug, error};
    use serde::{Serialize, Serializer};
    use serde_repr::Serialize_repr;

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

            let root = Self::generate_tiles(quadtree, world, &transformer, arg_minz, arg_maxz);
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
                let d = tile_bbox[3] - tile_bbox[0];
                if d < 0.0 {
                    debug!("d is negative in parent");
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
                    id: TileId::from(&quadtree.id),
                    bounding_volume,
                    geometric_error: d,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: None,
                    children: Some(tile_children),
                    implicit_tiling: None,
                }
            } else {
                // Compute the tile content bounding box <-- the bbox of all the cells in a tile
                let mut tile_content_bbox_qc = crate::spatial_structs::BboxQc([0, 0, 0, 0, 0, 0]);
                for cellid in quadtree.cells() {
                    let cell = world.grid.cell(cellid);
                    if !cell.feature_ids.is_empty() {
                        tile_content_bbox_qc = world.features[cell.feature_ids[0]].bbox_qc.clone();
                        break;
                    }
                }
                for cellid in quadtree.cells() {
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
                let d = tile_bbox[3] - tile_bbox[0];
                if d < 0.0 {
                    debug!("d is negative in child");
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
                    id: TileId::from(&quadtree.id),
                    bounding_volume,
                    geometric_error: 0.0,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content: Some(Content {
                        bounding_volume: Some(content_bounding_voume),
                        uri: format!("tiles/{}.glb", quadtree.id),
                    }),
                    children: None,
                    implicit_tiling: None,
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
                    id: TileId::new(cellid.column, cellid.row, 3),
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
                    implicit_tiling: None,
                };

                // LoD 1.3
                let tile_lod13 = Tile {
                    id: TileId::new(cellid.column, cellid.row, 2),
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
                    implicit_tiling: None,
                };

                // LoD 1.2
                root_children.push(Tile {
                    id: TileId::new(cellid.column, cellid.row, 1),
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
                    implicit_tiling: None,
                });
            }

            let root_volume = BoundingVolume::region_from_bbox(&grid.bbox, &transformer).unwrap();
            debug!("root bbox: {:?}", &grid.bbox);
            debug!("root boundingVolume: {:?}", &root_volume);
            let root_geometric_error = grid.bbox[3] - grid.bbox[0];

            let root = Tile {
                id: TileId::new(0, 0, 0),
                bounding_volume: root_volume,
                geometric_error: root_geometric_error,
                viewer_request_volume: None,
                refine: Some(Refinement::Replace),
                transform: None,
                content: None,
                children: Some(root_children),
                implicit_tiling: None,
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

        /// Flatten the tile hierarchy, visiting each tile in the quadtree.
        /// If 'levels_up' is provided, the tiles will be flattened only
        /// 'n levels upwards from the leaves', outputting only the flattened tiles
        /// (instead of the whole tree).
        pub fn flatten(&self, levels_up: Option<u16>) -> Vec<&Tile> {
            self.root.flatten(levels_up)
        }

        pub fn add_content(&mut self, levels_up: Option<u16>) {
            self.root.add_content_from_level(levels_up);
        }

        /// The number of levels in the quadtree.
        pub fn available_levels(&self) -> u16 {
            self.root.max_level()
        }

        /// Convert to implicit tiling.
        /// It modifies the tileset and deletes the explicit tiles.
        /// Expects that explicit tiling is already created.
        pub fn make_implicit(
            &mut self,
            grid: &crate::spatial_structs::SquareGrid,
            qtree: &crate::spatial_structs::QuadTree,
        ) {
            let subtree_sections: usize = 2;
            let subtree_levels =
                (self.available_levels() as f32 / subtree_sections as f32).ceil() as u16;
            let subtrees = Subtrees::default();
            let implicittiling = ImplicitTiling {
                subdivision_scheme: SubdivisionScheme::Quadtree,
                subtree_levels,
                available_levels: self.available_levels(),
                subtrees,
            };
            self.root.implicit_tiling = Some(implicittiling);

            // We want to have a sparse implicit tileset, which stores the availability
            //  of each tile. We do not have a full quadtree, because of the node
            //  capacity that is set.
            // In order to create this sparse implicit tileset, we need to have a full
            //  theoretical quadtree of the same height as the explicit quadtree
            //  (tileset), and for each node in the theoretical tree, we need to mark
            //  whether that node (tile) is indeed available and whether it has a
            //  content in the explicit tileset.
            let mut tiles_queue = VecDeque::new();
            tiles_queue.push_back(&self.root);
            let extent_side_length = grid.bbox[3] - grid.bbox[0];
            let grid_buffer = None;
            let mut level_cntr: u32 = 0;
            for section in 0..subtree_sections {
                let mut buffer_bytearray: Vec<u8> = Vec::new();
                // let mut tile_availability_bitstream = bv::bitvec![u8, bv::Msb0;];
                // let mut content_availability_bitstream = bv::bitvec![u8, bv::Msb0;];
                for level in 0..subtree_levels as usize {
                    let nr_tiles_current_level = 4_usize.pow(level_cntr);
                    let nr_tiles_current_level_side = 2_usize.pow(level_cntr);
                    level_cntr += 1;
                    let tile_side_length = extent_side_length / nr_tiles_current_level_side as f64;
                    let tile_side_length_u16 = tile_side_length as u16;

                    debug!("creating grid_for_level {}-{}", &section, &level);
                    let grid_for_level = crate::spatial_structs::SquareGrid::new(
                        &grid.bbox,
                        tile_side_length_u16,
                        grid.epsg,
                        grid_buffer,
                    );

                    // DEBUG
                    let mut file_grid =
                        File::create(format!("grid_for_level-{}-{}.tsv", &section, &level))
                            .unwrap();
                    for (cellid, _) in &grid_for_level {
                        let wkt = grid_for_level.cell_to_wkt(&cellid);
                        writeln!(file_grid, "{}\t{}", &cellid, wkt);
                    }

                    let mut mortoncodes: Vec<u64> = grid_for_level
                        .into_iter()
                        .map(|(cellid, _)| {
                            let [minx, miny, _minz, _maxx, _maxy, _maxz] =
                                grid_for_level.cell_bbox(&cellid);
                            crate::spatial_structs::interleave(&(minx as u64), &(miny as u64))
                        })
                        .collect();
                    mortoncodes.sort();
                    // let mut grid_cells_in_morton: HashMap<u64, crate::spatial_structs::CellId> =
                    //     HashMap::new();
                    // for (cellid, _) in &grid_for_level {
                    //     let [minx, miny, _minz, _maxx, _maxy, _maxz] =
                    //         grid_for_level.cell_bbox(&cellid);
                    //     let mortoncode = crate::spatial_structs::interleave(
                    //                     &(minx as u64),
                    //                     &(miny as u64),
                    //                 );
                    //     grid_cells_in_morton.insert(mortoncode, cellid);
                    // }
                    let mut grid_cells_in_morton: Vec<(
                        u64,
                        crate::spatial_structs::CellId,
                        f64,
                        f64,
                    )> = grid_for_level
                        .into_iter()
                        .map(|(cellid, _)| {
                            let [minx, miny, _minz, _maxx, _maxy, _maxz] =
                                grid_for_level.cell_bbox(&cellid);
                            (
                                // crate::spatial_structs::interleave(
                                //     &(minx as u64),
                                //     &(miny as u64),
                                // ),
                                crate::spatial_structs::interleave(
                                    &(cellid.column as u64),
                                    &(cellid.row as u64),
                                ),
                                cellid,
                                minx,
                                miny,
                            )
                        })
                        .collect();
                    grid_cells_in_morton.sort_by_key(|k| k.0);

                    for (cellid, _) in &grid_for_level {
                        let [minx, miny, _minz, _maxx, _maxy, _maxz] =
                            grid_for_level.cell_bbox(&cellid);
                        println!("[{},{}]", minx, miny);
                    }
                    for (cellid, _) in &grid_for_level {
                        println!("[{},{}]", cellid.column, cellid.row);
                    }

                    for (mc, cellid, minx, miny) in &grid_cells_in_morton {
                        let coords = crate::spatial_structs::deinterleave(mc);
                        println!("{}\t{:?}\t[{},{}]\t{}", mc, coords, minx, miny, cellid);
                    }
                    for (mc, cellid, minx, miny) in &grid_cells_in_morton {
                        let coords = crate::spatial_structs::deinterleave(mc);
                        println!("{:?}", coords);
                    }
                    for (mc, cellid, minx, miny) in &grid_cells_in_morton {
                        let coords = crate::spatial_structs::deinterleave(mc);
                        println!("[{},{}]", minx, miny);
                    }
                    for (mc, cellid, minx, miny) in &grid_cells_in_morton {
                        let coords = crate::spatial_structs::deinterleave(mc);
                        println!("{}", cellid);
                    }
                    debug!("grid_cells_in_morton {:?}", grid_cells_in_morton);
                    // for (i, mc) in mortoncodes.iter().enumerate() {
                    //     assert_eq!(*mc, grid_cells_in_morton[i].0)
                    // }

                    let nr_bytes = (nr_tiles_current_level as f32 / 8_f32).ceil() as u32;
                    let mut tile_availability_for_level = bv::bitvec![u8, bv::Msb0;];
                    tile_availability_for_level.resize((nr_bytes * 8) as usize, false);
                    let mut content_availability_for_level = bv::bitvec![u8, bv::Msb0;];
                    content_availability_for_level.resize((nr_bytes * 8) as usize, false);
                    debug!(
                        "tile_availability_for_level {}-{}: {:?}",
                        &section, &level, &tile_availability_for_level
                    );

                    // Check if the tile exists (== present in tileset)
                    let mut children_current_level: Vec<&Tile> = Vec::new();
                    for t in tiles_queue.iter() {
                        if let Some(ref ch) = t.children {
                            for c in ch {
                                let tileid = &c.id;
                                let qtree_nodeid: crate::spatial_structs::QuadTreeNodeId =
                                    tileid.into();
                                let cell = qtree.node(&qtree_nodeid).unwrap();
                                if cell.nr_items > 0 {
                                    children_current_level.push(c)
                                } else {
                                    debug!("found empty tile {}", tileid);
                                }
                            }
                        }
                    }

                    while let Some(tile) = tiles_queue.pop_front() {
                        let tileid = &tile.id;
                        let qtree_nodeid: crate::spatial_structs::QuadTreeNodeId = tileid.into();
                        let tile_bbox = qtree.node(&qtree_nodeid).unwrap().bbox(grid);
                        let tile_minx = tile_bbox[0];
                        let tile_miny = tile_bbox[1];
                        let tile_mortoncode = crate::spatial_structs::interleave(
                            &(tile_minx as u64),
                            &(tile_miny as u64),
                        );
                        // Set the tile and content available
                        if let Ok(tile_idx) = mortoncodes.binary_search(&tile_mortoncode) {
                            debug!("index in mortoncodes {}", tile_idx);
                            let (_, cellid, minx, miny) = grid_cells_in_morton[tile_idx];
                            debug!("tile {} matched grid cell {}", tileid, cellid);
                            let mut tile_bit =
                                tile_availability_for_level.get_mut(tile_idx).unwrap();
                            *tile_bit = true;
                            if tile.content.is_some() {
                                let mut content_bit =
                                    content_availability_for_level.get_mut(tile_idx).unwrap();
                                *content_bit = true;
                            }
                        } else {
                            error!(
                                "tile {} morton-code should match one of the generated cells",
                                tile.id
                            );
                        }
                    }
                    tiles_queue.extend(children_current_level);

                    // DEBUG
                    debug!(
                        "tile_availability_for_level {}-{}: {:?}",
                        &section, &level, &tile_availability_for_level
                    );
                    let mut file_implicit_tileset_at_level =
                        File::create(format!("implicit-level-{}-{}.tsv", section, level)).unwrap();
                    for (i, cellref) in grid_for_level.into_iter().enumerate() {
                        let wkt = grid_for_level.cell_to_wkt(&cellref.0);
                        let tile_available = tile_availability_for_level.get(i).unwrap();
                        let content_available = content_availability_for_level.get(i).unwrap();
                        writeln!(
                            file_implicit_tileset_at_level,
                            "{}\t{}\t{}\t{}",
                            &cellref.0,
                            tile_available.as_ref(),
                            content_available.as_ref(),
                            wkt
                        )
                        .unwrap();
                    }
                }
                // write subtree file
                let buffer = Buffer {
                    name: None,
                    byte_length: buffer_bytearray.len(),
                };
                let subtree = Subtree {
                    buffers: Some(vec![buffer]),
                    buffer_views: None,
                    tile_availability: Default::default(),
                    content_availability: None,
                    child_subtree_availability: Default::default(),
                };
            }

            self.root.children = None;
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
    pub struct Tile {
        #[serde(skip)]
        pub id: TileId,
        bounding_volume: BoundingVolume,
        pub geometric_error: GeometricError,
        #[serde(skip_serializing_if = "Option::is_none")]
        viewer_request_volume: Option<BoundingVolume>,
        #[serde(skip_serializing_if = "Option::is_none")]
        refine: Option<Refinement>,
        #[serde(skip_serializing_if = "Option::is_none")]
        transform: Option<Transform>,
        #[serde(skip_serializing_if = "Option::is_none")]
        content: Option<Content>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub children: Option<Vec<Tile>>,
        #[serde(skip_serializing_if = "Option::is_none")]
        implicit_tiling: Option<ImplicitTiling>,
    }

    impl Tile {
        fn flatten_recurse<'collect>(
            &'collect self,
            nodes: &mut Vec<&'collect Tile>,
            limit_upwards: &u16,
        ) {
            if self.id.level < *limit_upwards {
                if let Some(ref children) = self.children {
                    debug!("nr of children {}", children.len());
                    for child in children {
                        child.flatten_recurse(nodes, limit_upwards);
                    }
                }
            } else {
                nodes.push(self);
                if let Some(ref children) = self.children {
                    debug!("nr of children {}", children.len());
                    for child in children {
                        child.flatten_recurse(nodes, limit_upwards);
                    }
                }
            }
        }

        /// Flatten the tile hierarchy, visiting each tile in the quadtree.
        /// If 'levels_up' is provided, the tiles will be flattened only
        /// 'n levels upwards from the leaves', outputting only the flattened tiles
        /// (instead of the whole tree).
        pub fn flatten(&self, levels_up: Option<u16>) -> Vec<&Tile> {
            let max_level = self.max_level();
            let mut limit_upwards: u16 = 0;
            if let Some(limit) = levels_up {
                if limit < max_level {
                    limit_upwards = max_level - limit;
                }
            }
            let mut flat_tiles: Vec<&Tile> = Vec::new();
            self.flatten_recurse(&mut flat_tiles, &limit_upwards);
            flat_tiles
        }

        fn add_content_from_level(&mut self, levels_up: Option<u16>) {
            let max_level = self.max_level();
            let mut lower_limit: u16 = 0;
            if let Some(limit) = levels_up {
                if limit < max_level {
                    lower_limit = max_level - limit;
                }
            }
            let mut q = VecDeque::new();
            q.push_back(self);
            while let Some(node) = q.pop_front() {
                if node.id.level >= lower_limit {
                    node.add_content();
                }
                if let Some(ref mut children) = node.children {
                    for child in children.iter_mut() {
                        q.push_back(child);
                    }
                }
            }
        }

        fn max_level(&self) -> u16 {
            let mut max_level: u16 = 0;
            self.max_level_recurse(&mut max_level);
            max_level
        }

        fn max_level_recurse(&self, current_level: &mut u16) {
            if let Some(ref children) = self.children {
                for child in children {
                    child.max_level_recurse(current_level);
                }
            } else if self.id.level > *current_level {
                *current_level = self.id.level;
            }
        }

        // Adds `Content` to the Tile, by generating a content bounding volume from the
        // tile's bounding volume and a filepath from the tile.id.
        pub fn add_content(&mut self) {
            self.content = Some(Content {
                bounding_volume: Some(self.bounding_volume),
                uri: format!("tiles/{}.glb", self.id),
            })
        }
    }

    #[derive(Clone, Debug, Default, Eq, PartialEq)]
    pub struct TileId {
        x: usize,
        y: usize,
        level: u16,
    }

    impl TileId {
        pub fn new(x: usize, y: usize, level: u16) -> Self {
            Self { x, y, level }
        }
    }

    impl Display for TileId {
        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}/{}/{}", self.level, self.x, self.y)
        }
    }

    impl From<&crate::spatial_structs::QuadTreeNodeId> for TileId {
        fn from(value: &crate::spatial_structs::QuadTreeNodeId) -> Self {
            Self {
                x: value.x,
                y: value.y,
                level: value.level,
            }
        }
    }

    impl Into<crate::spatial_structs::QuadTreeNodeId> for &TileId {
        fn into(self) -> crate::spatial_structs::QuadTreeNodeId {
            crate::spatial_structs::QuadTreeNodeId::new(self.x, self.y, self.level)
        }
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

    /// Implicit tiling object.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#implicit-root-tile
    /// https://github.com/CesiumGS/3d-tiles/blob/1.1/specification/schema/tile.implicitTiling.schema.json
    #[derive(Serialize, Default, Debug)]
    #[serde(rename_all = "camelCase")]
    struct ImplicitTiling {
        subdivision_scheme: SubdivisionScheme,
        subtree_levels: u16,
        available_levels: u16,
        subtrees: Subtrees,
    }

    /// Implicit tiling subtree subdivision scheme.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#subdivision-scheme
    #[allow(dead_code)]
    #[derive(Serialize, Debug, Default)]
    #[serde(rename_all = "UPPERCASE")]
    enum SubdivisionScheme {
        #[default]
        Quadtree,
        Octree,
    }

    /// Implicit tiling subtrees.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#subtrees
    #[derive(Serialize, Debug)]
    struct Subtrees {
        uri: String,
    }

    impl Default for Subtrees {
        fn default() -> Self {
            Self {
                uri: String::from("subtrees/{level}/{x}/{y}.subtree"),
            }
        }
    }

    /// Implicit tiling subtree object.
    /// Metadata is not supported.
    /// https://github.com/CesiumGS/3d-tiles/blob/1.1/specification/schema/Subtree/subtree.schema.json
    #[derive(Serialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct Subtree {
        #[serde(skip_serializing_if = "Option::is_none")]
        buffers: Option<Vec<Buffer>>,
        #[serde(skip_serializing_if = "Option::is_none")]
        buffer_views: Option<Vec<BufferView>>,
        tile_availability: Availability,
        #[serde(skip_serializing_if = "Option::is_none")]
        content_availability: Option<Vec<Availability>>,
        child_subtree_availability: Availability,
    }

    #[derive(Serialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct Buffer {
        #[serde(skip_serializing_if = "Option::is_none")]
        name: Option<String>,
        byte_length: usize,
    }

    #[derive(Serialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct BufferView {
        buffer: u8,
        byte_offset: usize,
        byte_length: usize,
        #[serde(skip_serializing_if = "Option::is_none")]
        name: Option<String>,
    }

    /// Implicit tiling subtree availability.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#availability-1
    /// https://github.com/CesiumGS/3d-tiles/blob/1.1/specification/schema/Subtree/availability.schema.json
    #[derive(Serialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct Availability {
        #[serde(skip_serializing_if = "Option::is_none")]
        bitstream: Option<u8>,
        #[serde(skip_serializing_if = "Option::is_none")]
        available_count: Option<u16>,
        #[serde(skip_serializing_if = "Option::is_none")]
        constant: Option<AvailabilityConstant>,
    }

    /// Integer indicating whether all of the elements are available (1) or all are unavailable (0).
    #[allow(dead_code)]
    #[derive(Debug, Default, Serialize_repr)]
    #[repr(u8)]
    enum AvailabilityConstant {
        #[default]
        Unavailable = 0,
        Available = 1,
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use serde_json::to_string_pretty;
        use std::path::PathBuf;

        fn test_data_dir() -> PathBuf {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("resources")
                .join("data")
        }

        #[test]
        fn test_implicittiling() {
            // 85162.9 447106.8 85562.9 447706.8
            // let bbox: crate::spatial_structs::Bbox =
            //     [85162.9, 447106.8, -10.7, 85962.9, 447906.8, 320.5];
            // let grid = crate::spatial_structs::SquareGrid::new(&bbox, 200, 7415, Some(10.0));

            let mut world = crate::parser::World::new(
                test_data_dir()
                    .join("features_3dbag_5909")
                    .join("metadata.city.json"),
                test_data_dir()
                    .join("features_3dbag_5909")
                    .join("3dbag_v21031_7425c21b_5909_subset"),
                200,
                Some(vec![
                    crate::parser::CityObjectType::Building,
                    crate::parser::CityObjectType::BuildingPart,
                ]),
                None,
                None,
            )
            .unwrap();
            world.index_with_grid();

            world.export_grid();

            let quadtree = crate::spatial_structs::QuadTree::from_world(
                &world,
                crate::spatial_structs::QuadTreeCapacity::Vertices(15000),
            );
            quadtree.export(&world.grid).unwrap();

            let mut tileset = Tileset::from_quadtree(&quadtree, &world, None, None);

            tileset.make_implicit(&world.grid, &quadtree);

            let mut i = ImplicitTiling::default();
            println!("{}", serde_json::to_string(&i).unwrap());
        }

        #[test]
        fn test_availability() {
            let a = AvailabilityConstant::Available;
            assert_eq!("1", serde_json::to_string(&a).unwrap());
        }

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
