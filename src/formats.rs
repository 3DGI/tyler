//! Output formats for the tiles.
//! All format
// Copyright 2023 Balázs Dukai, Ravi Peters
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
    use log::{debug, error, info, warn};
    use morton_encoding::morton_encode;
    use serde::{Deserialize, Serialize};
    use serde_repr::{Deserialize_repr, Serialize_repr};

    use crate::proj::Proj;
    use crate::spatial_structs::{Bbox, CellId, QuadTree, QuadTreeNodeId, SquareGrid};

    /// [Tileset](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tileset).
    ///
    /// Not supported: `extras`.
    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
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
            let mut ser =
                serde_json::ser::Serializer::with_formatter(file_out, BoundingVolumeFormatter);
            self.serialize(&mut ser)?;
            Ok(())
        }

        pub fn export_bincode(
            &self,
            name: Option<&str>,
            output_dir: Option<&Path>,
        ) -> bincode::Result<()> {
            let file_name: &str = name.unwrap_or("tileset");
            let file = match output_dir {
                None => File::create(format!("{file_name}.bincode"))?,
                Some(outdir) => File::create(outdir.join(format!("{file_name}.bincode")))?,
            };
            bincode::serialize_into(file, self)
        }

        /// Write the tile boundingVolume and content boundingVolume for each level of the tileset
        /// to .tsv to the working directory.
        pub fn export(&self, output_dir: Option<&Path>) -> std::io::Result<()> {
            let mut q = VecDeque::new();
            q.push_back(&self.root);
            let mut tileset_level: u16 = self.root.id.level;
            let [mut outdir_tileset, mut outdir_tileset_content] = match output_dir {
                None => [Path::new(""), Path::new("")],
                Some(outdir) => [outdir, outdir],
            };
            let mut file_tileset =
                File::create(outdir_tileset.join(format!("tileset_level-{tileset_level}.tsv")))?;
            let mut file_tileset_content = File::create(
                outdir_tileset_content.join(format!("tileset_content_level-{tileset_level}.tsv")),
            )?;

            file_tileset
                .write_all("id\tlevel\thas_content\twkt\n".as_bytes())
                .expect("cannot write tileset TSV header");
            file_tileset_content
                .write_all("id\tlevel\twkt\n".as_bytes())
                .expect("cannot write tileset content TSV header");

            while let Some(tile) = q.pop_front() {
                if tile.id.level != tileset_level {
                    tileset_level = tile.id.level;
                    file_tileset = File::create(
                        outdir_tileset.join(format!("tileset_level-{tileset_level}.tsv")),
                    )?;
                    file_tileset_content = File::create(
                        outdir_tileset_content
                            .join(format!("tileset_content_level-{tileset_level}.tsv")),
                    )?;
                    file_tileset
                        .write_all("id\tlevel\thas_content\twkt\n".as_bytes())
                        .expect("cannot write tileset TSV header");
                    file_tileset_content
                        .write_all("id\tlevel\twkt\n".as_bytes())
                        .expect("cannot write tileset content TSV header");
                }
                let wkt = tile.bounding_volume.to_wkt();
                file_tileset
                    .write_all(
                        format!(
                            "{}\t{}\t{}\t{}\n",
                            tile.id,
                            tile.id.level,
                            tile.content.is_some(),
                            wkt
                        )
                        .as_bytes(),
                    )
                    .expect("cannot write tileset tile");
                if let Some(ref content) = tile.content {
                    if let Some(ref bv) = content.bounding_volume {
                        let wkt_content_bbox = bv.to_wkt();
                        file_tileset_content
                            .write_all(
                                format!("{}\t{}\t{}\n", tile.id, tile.id.level, wkt_content_bbox)
                                    .as_bytes(),
                            )
                            .expect("cannot write tileset tile content");
                    }
                }
                if let Some(ref children) = tile.children {
                    for child in children {
                        q.push_back(child);
                    }
                }
            }
            Ok(())
        }

        pub fn from_quadtree(
            quadtree: &QuadTree,
            world: &crate::parser::World,
            geometric_error_above_leaf: f64,
            arg_cellsize: u32,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
            content_bv_from_tile: bool,
        ) -> Self {
            let crs_from = format!("EPSG:{}", world.crs.to_epsg().unwrap());
            // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
            let crs_to = "EPSG:4978";
            let transformer = Proj::new_known_crs(&crs_from, crs_to, None).unwrap();
            // y-up to z-up transform needed because we are using gltf assets, which is y-up
            // https://github.com/CesiumGS/3d-tiles/tree/main/specification#y-up-to-z-up
            // let y_up_to_z_up = Transform([
            //     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            // ]);

            let root = Self::generate_tiles(
                quadtree,
                world,
                &transformer,
                geometric_error_above_leaf,
                arg_cellsize,
                arg_minz,
                arg_maxz,
                content_bv_from_tile,
            );
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
            quadtree: &QuadTree,
            world: &crate::parser::World,
            transformer: &Proj,
            geometric_error_above_leaf: f64,
            arg_cellsize: u32,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
            content_bv_from_tile: bool,
        ) -> Tile {
            if !quadtree.children.is_empty() {
                let tile_id = TileId::from(&quadtree.id);

                if quadtree.children.len() != 4 {
                    warn!("Quadtree does not have 4 children {:?}", &quadtree);
                }
                // Tile bounding volume
                // Set the bounding volume height from the grid height, which can be set with
                // an argument, or else calculated from the data (content).
                let mut tile_bbox = quadtree.bbox(&world.grid);
                // But it can happen with faulty data, eg. 3D Basisvoorziening,
                // that maxz is less than minz.
                if tile_bbox[5] < tile_bbox[2] {
                    debug!("Internal tile {tile_id} {:?} (in input CRS) bbox maxz {} is less than minz {}. Replacing maxz with minz + minz * 0.01.", &tile_bbox, tile_bbox[5], tile_bbox[2]);
                    tile_bbox[5] = tile_bbox[2] + tile_bbox[2] * 0.01;
                }
                let bounding_volume =
                    BoundingVolume::box_from_bbox(&tile_bbox, transformer).unwrap();

                let _s = serde_json::to_string(&bounding_volume).unwrap();
                debug!("{}", _s);

                // The geometric error of a tile is computed based on the specified error
                // for the nodes have leafs as children (assuming all leaf nodes are at the same level)
                let level_multiplier = (tile_bbox[3] - tile_bbox[0]) / (arg_cellsize as f64) - 2.0;
                let mut d = geometric_error_above_leaf * level_multiplier;
                let d_string = format!("{d:.2}");
                if d < 0.0 {
                    warn!("d is negative in internal tile {tile_id}");
                } else if d_string == *"0.00" {
                    // Because, for instance we have a —grid-cellsize 250, then a parent of the deepest level will have an edge length of 2 * 250.
                    // So for the 'level_multiplier' formula we get:
                    // 500 / 250 - 2.0 = 0
                    // Which then results in a 'd' of 0.
                    d = geometric_error_above_leaf;
                }
                let mut tile_children: Vec<Tile> = Vec::new();
                for child in quadtree.children.iter() {
                    tile_children.push(Self::generate_tiles(
                        child,
                        world,
                        transformer,
                        geometric_error_above_leaf,
                        arg_cellsize,
                        arg_minz,
                        arg_maxz,
                        content_bv_from_tile,
                    ));
                }
                Tile {
                    id: tile_id,
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
                let tile_id = TileId::from(&quadtree.id);

                // Tile bounding volume
                let mut tile_bbox = quadtree.bbox(&world.grid);
                if tile_bbox[5] < tile_bbox[2] {
                    // See explanation above
                    debug!("Leaf tile {tile_id} {:?} (in input CRS) bbox maxz {} is less than minz {}. Replacing maxz with minz + minz * 0.01.", &tile_bbox, tile_bbox[5], tile_bbox[2]);
                    tile_bbox[5] = tile_bbox[2] + tile_bbox[2] * 0.01;
                }
                let mut bounding_volume =
                    BoundingVolume::box_from_bbox(&tile_bbox, transformer).unwrap();
                let mut content: Option<Content> = None;

                if quadtree.nr_items > 0 {
                    let mut tile_content_bbox_rw =
                        quadtree.node_content_bbox(&world, arg_minz, arg_maxz);
                    if content_bv_from_tile {
                        tile_content_bbox_rw = tile_bbox;
                    } else {
                        // Stretch the tile bbox so that it covers the content bbox
                        tile_content_bbox_rw[0..3]
                            .iter()
                            .enumerate()
                            .for_each(|(i, min_c)| {
                                if min_c < &tile_bbox[i] {
                                    tile_bbox[i] = min_c - min_c * 0.01;
                                }
                            });
                        tile_content_bbox_rw[3..6]
                            .iter()
                            .enumerate()
                            .for_each(|(i0, max_c)| {
                                let i = 3 + i0;
                                if max_c > &tile_bbox[i] {
                                    tile_bbox[i] = max_c + max_c * 0.01;
                                }
                            });
                    }

                    if tile_content_bbox_rw[5] < tile_content_bbox_rw[2] {
                        // See explanation above
                        debug!("Leaf tile content {tile_id} {:?} (in input CRS) bbox maxz {} is less than minz {}. Replacing maxz with minz + minz * 0.01.", &tile_content_bbox_rw, tile_content_bbox_rw[5], tile_content_bbox_rw[2]);
                        tile_content_bbox_rw[5] =
                            tile_content_bbox_rw[2] + tile_content_bbox_rw[2] * 0.01;
                    }
                    let content_bounding_volume =
                        BoundingVolume::box_from_bbox(&tile_content_bbox_rw, transformer).unwrap();

                    content = Some(Content {
                        bounding_volume: Some(content_bounding_volume),
                        uri: format!("tiles/{}.glb", quadtree.id),
                    });
                }

                Tile {
                    id: tile_id,
                    bounding_volume,
                    geometric_error: 0.0,
                    viewer_request_volume: None,
                    refine: Some(Refinement::Replace),
                    transform: None,
                    content,
                    children: None,
                    implicit_tiling: None,
                }
            }
        }

        #[allow(dead_code)]
        pub fn from_grid(
            grid: &SquareGrid,
            citymodel: &crate::parser::CityJSONMetadata,
            feature_set: &crate::parser::FeatureSet,
        ) -> Self {
            let crs_from = format!(
                "EPSG:{}",
                citymodel.metadata.reference_system.to_epsg().unwrap()
            );
            // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
            let crs_to = "EPSG:4978";
            let transformer = Proj::new_known_crs(&crs_from, crs_to, None).unwrap();

            let mut root_children: Vec<Tile> = Vec::with_capacity(grid.length * grid.length);
            for (cellid, cell) in grid {
                if cell.feature_ids.is_empty() {
                    // Empty cell, don't create tiles for it
                    continue;
                }

                let mut content_bbox_qc = feature_set[cell.feature_ids[0]].bbox_qc.clone();
                for fi in cell.feature_ids.iter() {
                    content_bbox_qc.update_with(&feature_set[*fi].bbox_qc);
                }
                let content_bbox_rw = content_bbox_qc.to_bbox(&citymodel.transform, None, None);
                let content_bounding_voume =
                    BoundingVolume::box_from_bbox(&content_bbox_rw, &transformer).unwrap();

                let mut cell_bbox = grid.cell_bbox(&cellid);
                // Set the bounding volume height from the content height
                cell_bbox[2] = content_bbox_rw[2];
                cell_bbox[5] = content_bbox_rw[5];
                let bounding_volume =
                    BoundingVolume::box_from_bbox(&cell_bbox, &transformer).unwrap();

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

            let root_volume = BoundingVolume::box_from_bbox(&grid.bbox, &transformer).unwrap();
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

        pub fn collect_leaves(&self) -> Vec<&Tile> {
            self.root.collect_leaves()
        }

        pub fn add_content(&mut self, levels_up: Option<u16>) {
            self.root.add_content_from_level(levels_up);
        }

        /// The number of levels in the quadtree, which is `max_level + 1`.
        pub fn available_levels(&self) -> u16 {
            self.root.max_level() + 1
        }

        /// Convert to implicit tiling.
        /// It modifies the tileset and deletes the explicit tiles.
        /// Expects that explicit tiling is already created.
        pub fn make_implicit(
            &mut self,
            grid: &SquareGrid,
            qtree: &QuadTree,
            grid_export: bool,
            subtrees_dir: Option<&str>,
            output_dir_debug: Option<&Path>,
        ) -> (Vec<(Tile, TileId)>, Vec<(TileId, Vec<u8>)>) {
            let mut subtrees_vec: Vec<(TileId, Vec<u8>)> = Vec::new();
            let mut flat_tiles_with_content: Vec<(Tile, TileId)> = Vec::new();
            // https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc37
            let subtree_sections: usize = 1;
            let subtree_levels =
                (self.available_levels() as f32 / subtree_sections as f32).ceil() as u16;
            let subtrees = match subtrees_dir {
                None => Subtrees::default(),
                Some(dirname) => Subtrees::new(dirname),
            };
            let implicittiling = ImplicitTiling {
                subdivision_scheme: SubdivisionScheme::Quadtree,
                subtree_levels,
                available_levels: self.available_levels(),
                subtrees,
            };
            debug!("{:?}", &implicittiling);
            self.root.implicit_tiling = Some(implicittiling);

            // We want to have a sparse implicit tileset, which stores the availability
            //  of each tile. We do not have a full quadtree, because of the node
            //  capacity that is set.
            // In order to create this sparse implicit tileset, we need to have a full
            //  theoretical quadtree of the same height as the explicit quadtree
            //  (tileset), and for each node in the theoretical tree, we need to mark
            //  whether that node (tile) is indeed available and whether it has a
            //  content in the explicit tileset.

            let grid_epsg = grid.epsg;
            let level_subtree_root: u32 = 0; // subtree root level
            let mut subtree_queue = VecDeque::new();
            let rootid = &self.root.id;
            let cellid: CellId = rootid.into();
            subtree_queue.push_back((level_subtree_root, cellid, &self.root));

            while let Some((level_subtree_root, cellid, tile)) = subtree_queue.pop_front() {
                let subtree_id = TileId::new(cellid.column, cellid.row, level_subtree_root as u16);
                let mut buffer_vec: Vec<u8> = Vec::new();
                let mut tile_availability_bitstream: bv::BitVec<u8, bv::Lsb0> = bv::BitVec::new();
                let mut content_availability_bitstream: bv::BitVec<u8, bv::Lsb0> =
                    bv::BitVec::new();

                let tileid = &tile.id;
                let qtree_nodeid: QuadTreeNodeId = tileid.into();
                let tile_bbox = qtree.node(&qtree_nodeid).unwrap().bbox(grid);
                let extent_width = tile_bbox[3] - tile_bbox[0];

                let mut tiles_queue = VecDeque::new();
                tiles_queue.push_back(tile);
                let mut level_quadtree: u32 = 0;
                for level_subtree in 0..subtree_levels as u32 {
                    level_quadtree = level_subtree_root + level_subtree;
                    // The number of tiles on the current level of the full quadtree
                    let _nr_tiles = 4_usize.pow(level_quadtree);
                    // The number of tiles on the current level within the subtree. Each subtree
                    //  with a single root tile, on level 0. Regardless where the subtree is in the
                    //  full quadtree hierarchy.
                    let nr_tiles_subtree = 4_usize.pow(level_subtree);
                    let grid_coordinate_map = Self::grid_coordinate_map(
                        level_subtree,
                        extent_width,
                        &tile_bbox,
                        grid_epsg,
                        grid_export,
                        output_dir_debug,
                    );
                    let grid_coordinate_map_global = Self::grid_coordinate_map(
                        level_quadtree,
                        extent_width,
                        &tile_bbox,
                        grid_epsg,
                        false,
                        None,
                    );

                    let mut children_current_level: Vec<&Tile> = Vec::new();
                    for t in tiles_queue.iter() {
                        if let Some(ref ch) = t.children {
                            for c in ch {
                                let tileid = &c.id;
                                let qtree_nodeid: QuadTreeNodeId = tileid.into();
                                let cell = qtree.node(&qtree_nodeid).unwrap();
                                if cell.nr_items > 0 {
                                    children_current_level.push(c)
                                }
                            }
                        }
                    }

                    let mut tile_availability_for_level: bv::BitVec<u8, bv::Lsb0> =
                        bv::BitVec::new();
                    tile_availability_for_level.resize(nr_tiles_subtree, false);
                    let mut content_availability_for_level: bv::BitVec<u8, bv::Lsb0> =
                        bv::BitVec::new();
                    content_availability_for_level.resize(nr_tiles_subtree, false);

                    let mut tile_availability_for_level_vec: Vec<bool> = Vec::new();
                    tile_availability_for_level_vec.resize(nr_tiles_subtree, false);
                    let mut content_availability_for_level_vec: Vec<bool> = Vec::new();
                    content_availability_for_level_vec.resize(nr_tiles_subtree, false);

                    // FIXME DEBUG
                    let mut tileids_contiguous_vec: Vec<String> = Vec::new();
                    tileids_contiguous_vec.resize(nr_tiles_subtree, String::default());

                    while let Some(tile) = tiles_queue.pop_front() {
                        let tile_corner_coord = Self::tile_corner_coordinate(grid, qtree, tile);
                        // Set the tile and content available
                        if let Some((_cellid_grid_level, i_z_curve)) =
                            grid_coordinate_map.get(&tile_corner_coord)
                        {
                            if let Some((cellid_grid_global, ..)) =
                                grid_coordinate_map_global.get(&tile_corner_coord)
                            {
                                tile_availability_for_level.set(*i_z_curve, true);
                                tile_availability_for_level_vec[*i_z_curve] = true;
                                if tile.content.is_some() {
                                    content_availability_for_level.set(*i_z_curve, true);
                                    content_availability_for_level_vec[*i_z_curve] = true;
                                    let tileid_continuous = TileId::new(
                                        cellid_grid_global.column,
                                        cellid_grid_global.row,
                                        level_quadtree as u16,
                                    );
                                    tileids_contiguous_vec[*i_z_curve] =
                                        tileid_continuous.to_string();
                                    flat_tiles_with_content.push((tile.clone(), tileid_continuous));
                                }
                            } else {
                                debug!(
                                    "could not locate tile {} in grid_coordinate_map_global",
                                    tile.id
                                );
                            }
                        } else {
                            debug!("could not locate tile {} in grid_coordinate_map", tile.id);
                        }
                    }

                    if grid_export {
                        let nr_tiles = 4_usize.pow(level_subtree);
                        // Grid for the current level
                        let tile_width = (extent_width / (nr_tiles as f64).sqrt()) as u32;
                        let grid_for_level =
                            SquareGrid::new(&tile_bbox, tile_width, grid_epsg, None);
                        let outdir = output_dir_debug.unwrap_or(Path::new(""));
                        let filename = outdir.join(format!(
                            "implicit-level-{}-{}-{}.tsv",
                            &level_quadtree, &tile.id.x, &tile.id.y
                        ));
                        debug!(
                            "Exporting the subtree {}/{}/{} to TSV file to {}",
                            &level_quadtree,
                            &tile.id.x,
                            &tile.id.y,
                            filename.display()
                        );
                        let mut file_implicit_tileset_at_level = File::create(&filename).unwrap();
                        writeln!(
                            file_implicit_tileset_at_level,
                            "cell_id\ttile_id_subtree\ttile_available\tcontent_available\twkt",
                        )
                        .unwrap();
                        for (cellid_grid_level, i_z_curve) in grid_coordinate_map.values() {
                            let wkt = grid_for_level.cell_to_wkt(cellid_grid_level);
                            let va = tile_availability_for_level.get(*i_z_curve);
                            let vc = content_availability_for_level.get(*i_z_curve);
                            let tile_id_subtree = &tileids_contiguous_vec[*i_z_curve];
                            if va.is_none() {
                                error!("tileAvailability bitstream is inconsistent, there is no value at index {i_z_curve}");
                            };
                            if vc.is_none() {
                                error!("contentAvailability bitstream is inconsistent, there is no value at index {i_z_curve}");
                            }
                            let tile_available = va.unwrap();
                            let content_available = vc.unwrap();
                            writeln!(
                                file_implicit_tileset_at_level,
                                "{}\t{}\t{}\t{}\t{}",
                                cellid_grid_level,
                                tile_id_subtree,
                                tile_available.as_ref(),
                                content_available.as_ref(),
                                wkt
                            )
                            .unwrap();
                        }
                    }

                    tile_availability_for_level.set_uninitialized(false);
                    content_availability_for_level.set_uninitialized(false);
                    tile_availability_bitstream.extend_from_bitslice(&tile_availability_for_level);
                    content_availability_bitstream
                        .extend_from_bitslice(&content_availability_for_level);
                    tiles_queue.extend(children_current_level);
                }

                let level_child_subtree = level_quadtree + 1;
                let nr_tiles_child_level = 4_usize.pow(level_child_subtree);
                let nr_tiles_total_subtree = (4_usize.pow(subtree_levels as u32) - 1) / 3;
                assert_eq!(
                    tile_availability_bitstream.len(),
                    nr_tiles_total_subtree,
                    "tileAvailability bitstream must have {} elements, but it has {}",
                    &nr_tiles_total_subtree,
                    &tile_availability_bitstream.len()
                );
                let grid_coordinate_map = Self::grid_coordinate_map(
                    level_child_subtree,
                    extent_width,
                    &tile_bbox,
                    grid_epsg,
                    false,
                    None,
                );

                let mut child_subtree_availability_bitstream: bv::BitVec<u8, bv::Lsb0> =
                    bv::BitVec::new();
                child_subtree_availability_bitstream.resize(nr_tiles_child_level, false);
                for child in tiles_queue.iter() {
                    let tile_corner_coord = Self::tile_corner_coordinate(grid, qtree, child);
                    if let Some((cellid_grid_level, i_z_curve)) =
                        grid_coordinate_map.get(&tile_corner_coord)
                    {
                        child_subtree_availability_bitstream.set(*i_z_curve, true);
                        subtree_queue.push_back((level_child_subtree, *cellid_grid_level, child));
                    } else {
                        debug!(
                            "could not locate tile {} in grid_for_child_subtree",
                            child.id
                        );
                    }
                }

                // Bufferviews
                let mut bufferviews: Vec<BufferView> = Vec::with_capacity(3);
                let mut bufferview_idx: usize = 0;

                let tile_availability =
                    Self::create_availability(bufferview_idx, &mut tile_availability_bitstream);
                if tile_availability.constant.is_none() {
                    // add padding in buffer_vec to next 8 byte boundary
                    Self::add_padding(&mut buffer_vec, 8);
                    Self::add_bitstream(
                        &mut buffer_vec,
                        &mut bufferviews,
                        tile_availability_bitstream,
                    );
                    bufferview_idx += 1;
                }

                let content_availability =
                    Self::create_availability(bufferview_idx, &mut content_availability_bitstream);
                if content_availability.constant.is_none() {
                    // add padding in buffer_vec to next 8 byte boundary
                    Self::add_padding(&mut buffer_vec, 8);
                    Self::add_bitstream(
                        &mut buffer_vec,
                        &mut bufferviews,
                        content_availability_bitstream,
                    );
                    bufferview_idx += 1;
                }

                let child_subtree_availability = Self::create_availability(
                    bufferview_idx,
                    &mut child_subtree_availability_bitstream,
                );
                if child_subtree_availability.constant.is_none() {
                    // add padding in buffer_vec to next 8 byte boundary
                    Self::add_padding(&mut buffer_vec, 8);
                    Self::add_bitstream(
                        &mut buffer_vec,
                        &mut bufferviews,
                        child_subtree_availability_bitstream,
                    );
                }

                // pad our buffer_vec to have a length that is a multiple of 8 bytes
                Self::add_padding(&mut buffer_vec, 8);

                debug!("Writing subtree {}", &subtree_id);
                let buffer = Buffer {
                    name: None,
                    byte_length: buffer_vec.len(),
                };
                let subtree = Subtree {
                    buffers: Some(vec![buffer]),
                    buffer_views: Some(bufferviews),
                    tile_availability,
                    content_availability: Some(vec![content_availability]),
                    child_subtree_availability,
                };
                let mut subtree_json = serde_json::to_string(&subtree)
                    .expect("failed to serialize the subtree to json");
                let remainder = subtree_json.as_bytes().len() % 8;
                let mut padding = 0;
                if remainder > 0 {
                    padding = 8 - remainder;
                }
                for _i in 0..padding {
                    subtree_json += " ";
                }
                let subtree_json_bytes = subtree_json.as_bytes();

                let mut subtree_bytes: Vec<u8> = Vec::new();
                // Header
                let mut header: Vec<u8> = Vec::new();
                let magic = 0x74627573u32.to_le_bytes(); // subt
                let version = 1_u32.to_le_bytes();
                let json_byte_length = (subtree_json_bytes.len() as u64).to_le_bytes();
                let buffer_byte_length = ((buffer_vec.len()) as u64).to_le_bytes();
                header.extend_from_slice(&magic);
                header.extend_from_slice(&version);
                header.extend_from_slice(&json_byte_length);
                header.extend_from_slice(&buffer_byte_length);
                if header.len() != 24 {
                    warn!(
                        "Subtree {} binary header must be 24 bytes long, it is {}",
                        &subtree_id,
                        header.len()
                    );
                }
                subtree_bytes.extend(header);

                // Content
                subtree_bytes.extend(subtree_json_bytes);
                let buffer_bytes: Vec<u8> = buffer_vec
                    .iter()
                    .map(|i| i.to_le_bytes())
                    .flat_map(|bytearray| bytearray.to_vec())
                    .collect();
                subtree_bytes.extend(buffer_bytes);
                subtrees_vec.push((subtree_id, subtree_bytes));
            }

            self.root.content = Some(Content {
                bounding_volume: None,
                uri: "tiles/{level}/{x}/{y}.glb".to_string(),
            });
            self.root.children = None;
            (flat_tiles_with_content, subtrees_vec)
        }

        fn add_padding(buffer_vec: &mut Vec<u8>, align_by: usize) {
            let padding = (align_by - (buffer_vec.len() % align_by)) % align_by;
            for _i in 0..padding {
                buffer_vec.push(0);
            }
        }

        fn add_bitstream(
            buffer_vec: &mut Vec<u8>,
            bufferviews: &mut Vec<BufferView>,
            availability_bitstream: bv::BitVec<u8, bv::Lsb0>,
        ) {
            let availability_vec = availability_bitstream.into_vec();
            bufferviews.push(BufferView {
                buffer: 0,
                byte_offset: buffer_vec.len(),
                byte_length: availability_vec.len(),
                name: None,
            });
            buffer_vec.extend(availability_vec);
        }

        fn create_availability(
            bf_availability: usize,
            availability_bitstream: &mut bv::BitVec<u8, bv::Lsb0>,
        ) -> Availability {
            availability_bitstream.set_uninitialized(false);
            if availability_bitstream.not_any() {
                Availability {
                    bitstream: None,
                    available_count: None,
                    constant: Some(AvailabilityConstant::Unavailable),
                }
            } else if availability_bitstream.all() {
                Availability {
                    bitstream: None,
                    available_count: None,
                    constant: Some(AvailabilityConstant::Available),
                }
            } else {
                Availability {
                    bitstream: Some(bf_availability),
                    available_count: Some(availability_bitstream.count_ones()),
                    constant: None,
                }
            }
        }

        fn tile_corner_coordinate(grid: &SquareGrid, qtree: &QuadTree, tile: &Tile) -> String {
            let tileid = &tile.id;
            let qtree_nodeid: QuadTreeNodeId = tileid.into();
            let tile_bbox = qtree.node(&qtree_nodeid).unwrap().bbox(grid);
            let [minx, miny, ..] = tile_bbox;
            format!("{:.0},{:.0}", minx, miny)
        }

        /// Build a map of grid-cell-corner-coorinates and cell ID-s.
        /// It is used for matching the "theoretical grid" to the quadtree nodes, based
        /// on the coordinates.
        fn grid_coordinate_map(
            level_current: u32,
            extent_width: f64,
            bbox: &Bbox,
            epsg: u16,
            do_export: bool,
            output_dir_debug: Option<&Path>,
        ) -> HashMap<String, (CellId, usize)> {
            let nr_tiles = 4_usize.pow(level_current);

            // Grid for the current level
            let tile_width = (extent_width / (nr_tiles as f64).sqrt()) as u32;
            let grid_for_level = SquareGrid::new(bbox, tile_width, epsg, None);

            // Map of:
            //  - x,y coordinate of the min coordinate of the lower-left cell
            //  - (cell ID, index in the z-order array)
            let mut grid_for_level_corner_coords: HashMap<String, (CellId, usize)> = HashMap::new();

            let mut mortoncodes: Vec<(u128, CellId)> = grid_for_level
                .into_iter()
                .map(|(cellid, _)| {
                    (
                        morton_encode([cellid.row as u64, cellid.column as u64]),
                        cellid,
                    )
                })
                .collect();
            mortoncodes.sort_by_key(|k| k.0);

            for (i, (_mc, cellid)) in mortoncodes.iter().enumerate() {
                let [minx, miny, ..] = grid_for_level.cell_bbox(cellid);
                // Since the input for grid_cellsize is u16 and expected to be in the range
                //  of several (hundreds) of meters, we don't care about decimal precision.
                let corner_coord_string = format!("{:.0},{:.0}", minx, miny);
                grid_for_level_corner_coords.insert(corner_coord_string, (*cellid, i));
            }

            if do_export {
                // DEBUG
                let first_cell = grid_for_level.into_iter().next().unwrap();
                let outdir = output_dir_debug.unwrap_or(Path::new(""));
                let mut file_grid = File::create(outdir.join(format!(
                    "grid_for_level-{}-{}.tsv",
                    &grid_for_level.length, &first_cell.0
                )))
                .unwrap();
                writeln!(file_grid, "cell_id\twkt").unwrap();
                for (cellid, _) in &grid_for_level {
                    let wkt = grid_for_level.cell_to_wkt(&cellid);
                    writeln!(file_grid, "{}\t{}", &cellid, wkt).unwrap();
                }

                let mut file_grid_morton = File::create(outdir.join(format!(
                    "grid_for_level_morton-{}-{}.tsv",
                    &grid_for_level.length, &first_cell.0
                )))
                .unwrap();
                writeln!(file_grid_morton, "idx\tcell_id\twkt").unwrap();
                for (i, (_mc, cellid)) in mortoncodes.iter().enumerate() {
                    let [minx, miny, ..] = grid_for_level.cell_bbox(cellid);
                    writeln!(
                        file_grid_morton,
                        "{}\t{}\tPOINT({} {})",
                        i, &cellid, minx, miny
                    )
                    .unwrap();
                }
            }

            grid_for_level_corner_coords
        }

        /// Prune the tileset by removing the tiles in `tiles_to_remove`.
        /// In addition, it also removes that with `nr_items == 0`.
        pub fn prune(&mut self, tiles_to_remove: &Vec<Tile>, qtree: &QuadTree) {
            self.root.prune(tiles_to_remove, qtree);
        }

        pub fn debug_bounding_volume_string(&self) -> String {
            serde_json::to_string(&self.root.bounding_volume).unwrap()
        }
    }

    /// [Asset](https://github.com/CesiumGS/3d-tiles/tree/main/specification#asset).
    ///
    /// Not supported: `extensions, extras`.
    #[derive(Serialize, Deserialize, Debug, Clone)]
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
    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
    #[serde(rename_all = "camelCase")]
    struct Properties {
        maximum: f64,
        minimum: f64,
    }

    type Extensions = HashMap<ExtensionName, Extension>;

    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
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
    #[derive(Serialize, Deserialize, Default, Debug, Clone, PartialEq, Eq, Hash)]
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
    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
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

    /// Tile equality is evaluated on the tile ID.
    impl PartialEq for Tile {
        fn eq(&self, other: &Self) -> bool {
            self.id == other.id
        }
    }

    /// Tile equality is evaluated on the tile ID.
    impl Eq for Tile {}

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

        fn collect_leaves_recurse<'collect>(&'collect self, leaves: &mut Vec<&'collect Tile>) {
            if let Some(ref children) = self.children {
                for child in children {
                    child.collect_leaves_recurse(leaves);
                }
            } else {
                leaves.push(self);
            }
        }

        pub fn collect_leaves(&self) -> Vec<&Self> {
            let mut leaves: Vec<&Tile> = Vec::new();
            self.collect_leaves_recurse(&mut leaves);
            leaves
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

        fn prune(&mut self, tiles_to_remove: &Vec<Tile>, qtree: &QuadTree) {
            if let Some(mut children) = self.children.take() {
                let mut children_new: Vec<Tile> = Vec::with_capacity(4);
                for child in children.iter_mut() {
                    if !tiles_to_remove.contains(&*child) {
                        let tileid: &TileId = &child.id;
                        let qtree_nodeid: QuadTreeNodeId = tileid.into();
                        if let Some(qtree_node) = qtree.node(&qtree_nodeid) {
                            if qtree_node.nr_items > 0 {
                                child.prune(tiles_to_remove, qtree);
                                children_new.push(child.clone());
                            }
                        } else {
                            error!("Did not find matching QuadTree node for TileId {}", tileid);
                        }
                    }
                }
                self.children = Some(children_new);
            }
        }
    }

    #[derive(Clone, Debug, Default, Eq, PartialEq)]
    pub struct TileId {
        pub(crate) x: usize,
        pub(crate) y: usize,
        pub(crate) level: u16,
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

    impl From<&QuadTreeNodeId> for TileId {
        fn from(value: &QuadTreeNodeId) -> Self {
            Self {
                x: value.x,
                y: value.y,
                level: value.level,
            }
        }
    }

    impl From<&TileId> for QuadTreeNodeId {
        fn from(val: &TileId) -> Self {
            QuadTreeNodeId::new(val.x, val.y, val.level)
        }
    }

    impl From<&TileId> for CellId {
        fn from(val: &TileId) -> Self {
            CellId {
                column: val.x,
                row: val.y,
            }
        }
    }

    /// Format the BoundingVolume coordinates to 6 decimal places in the JSON output.
    /// 6 decimal places, because that gives 0.11112m precision.
    /// See https://wiki.openstreetmap.org/wiki/Precision_of_coordinates
    struct BoundingVolumeFormatter;

    impl serde_json::ser::Formatter for BoundingVolumeFormatter {
        fn write_f64<W>(&mut self, writer: &mut W, value: f64) -> std::io::Result<()>
        where
            W: ?Sized + Write,
        {
            write!(writer, "{:.6}", value)
        }
    }

    /// [boundingVolume](https://github.com/CesiumGS/3d-tiles/tree/main/specification#bounding-volumes).
    #[allow(dead_code)]
    #[derive(Serialize, Deserialize, Debug, Copy, Clone)]
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
    impl From<&Bbox> for BoundingVolume {
        fn from(bbox: &Bbox) -> Self {
            unimplemented!()
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
        #[allow(dead_code)]
        fn box_from_bbox(
            bbox: &Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            // An empirical value from playing around in Cesium. Without this shift, the
            // Box center is placed way too high with implicit tiling.
            let magic_z_shift = 0.0;
            // Input CRS
            let dx = bbox[3] - bbox[0];
            let dy = bbox[4] - bbox[1];
            let dz = bbox[5] - bbox[2];
            let center: [f64; 3] = [
                bbox[0] + dx * 0.5,
                bbox[1] + dy * 0.5,
                (bbox[2] + dz * 0.5) + magic_z_shift,
            ];
            // ????
            // The center points on the max. surface of the bbox
            let x_max_pt: [f64; 3] = [bbox[3], center[1], center[2]];
            let y_max_pt: [f64; 3] = [center[0], bbox[4], center[2]];
            let z_max_pt: [f64; 3] = [center[0], center[1], bbox[5] + magic_z_shift];

            // Target CRS (ECEF)
            let center_ecef = transformer.convert((center[0], center[1], center[2]))?;
            let x_pt_ecef = transformer.convert((x_max_pt[0], x_max_pt[1], x_max_pt[2]))?;
            let y_pt_ecef = transformer.convert((y_max_pt[0], y_max_pt[1], y_max_pt[2]))?;
            let z_pt_ecef = transformer.convert((z_max_pt[0], z_max_pt[1], z_max_pt[2]))?;

            // Vectors that define the size and orientation of the OBB
            let vx = (
                x_pt_ecef.0 - center_ecef.0,
                x_pt_ecef.1 - center_ecef.1,
                x_pt_ecef.2 - center_ecef.2,
            );
            let vy = (
                y_pt_ecef.0 - center_ecef.0,
                y_pt_ecef.1 - center_ecef.1,
                y_pt_ecef.2 - center_ecef.2,
            );
            let vz = (
                z_pt_ecef.0 - center_ecef.0,
                z_pt_ecef.1 - center_ecef.1,
                z_pt_ecef.2 - center_ecef.2,
            );

            Ok(Self::Box([
                center_ecef.0,
                center_ecef.1,
                center_ecef.2,
                vx.0,
                vx.1,
                vx.2,
                vy.0,
                vy.1,
                vy.2,
                vz.0,
                vz.1,
                vz.2,
            ]))
        }

        fn region_from_bbox(
            bbox: &Bbox,
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

        /// Cast to 2D WKT
        fn to_wkt(&self) -> String {
            let [minx, miny, _minz, maxx, maxy, _maxz] = match self {
                BoundingVolume::Box(bbox) => {
                    let center = &bbox[0..3];
                    let vx = &bbox[3..6];
                    let vy = &bbox[6..9];
                    // 3┌───────▲───────┐2
                    //  │    -X │       │
                    //  │       │       │
                    //  │       │     Y │
                    //  ◄───────┼───────►
                    //  │ -Y    │       │
                    //  │       │       │
                    //  │       │  X    │
                    // 4└───────▼───────┘1
                    // Corner vectors
                    let corner_1_v = [vx[0] + vy[0], vx[1] + vy[1], vx[2] + vy[2]];
                    let corner_2_v = [vy[0] + -vx[0], vy[1] + -vx[1], vy[2] + -vx[2]];
                    let corner_3_v = [-corner_1_v[0], -corner_1_v[1], -corner_1_v[2]];
                    let corner_4_v = [-corner_2_v[0], -corner_2_v[1], -corner_2_v[2]];
                    // Lot of unnecessary iterations and array allocations here, but we only use
                    // WKT for debugging and rather have things here explicit here, for clarity.
                    let corner_points: Vec<[f64; 3]> =
                        [corner_1_v, corner_2_v, corner_3_v, corner_4_v]
                            .iter()
                            .map(|corner| {
                                [
                                    center[0] + corner[0],
                                    center[1] + corner[1],
                                    center[2] + corner[2],
                                ]
                            })
                            .collect();
                    let minx = corner_points.iter().map(|a| a[0]).reduce(f64::min).unwrap();
                    let miny = corner_points.iter().map(|a| a[1]).reduce(f64::min).unwrap();
                    let minz = corner_points.iter().map(|a| a[2]).reduce(f64::min).unwrap();
                    let maxx = corner_points.iter().map(|a| a[0]).reduce(f64::max).unwrap();
                    let maxy = corner_points.iter().map(|a| a[1]).reduce(f64::max).unwrap();
                    let maxz = corner_points.iter().map(|a| a[2]).reduce(f64::max).unwrap();
                    [minx, miny, minz, maxx, maxy, maxz]
                }
                BoundingVolume::Region(bbox) => [
                    bbox[0].to_degrees(),
                    bbox[1].to_degrees(),
                    bbox[4],
                    bbox[2].to_degrees(),
                    bbox[3].to_degrees(),
                    bbox[5],
                ],
                BoundingVolume::Sphere(_) => {
                    unimplemented!()
                }
            };
            format!(
                "POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
                minx = minx,
                miny = miny,
                maxx = maxx,
                maxy = maxy
            )
        }
    }

    /// [Tile.refine](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tilerefine).
    #[allow(dead_code)]
    #[derive(Serialize, Deserialize, Debug, Clone)]
    #[serde(rename_all = "UPPERCASE")]
    enum Refinement {
        Add,
        Replace,
    }

    /// [Tile.transform](https://github.com/CesiumGS/3d-tiles/tree/main/specification#tiletransform)
    #[derive(Serialize, Deserialize, Debug, Copy, Clone)]
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
    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
    #[serde(rename_all = "camelCase")]
    struct Content {
        #[serde(skip_serializing_if = "Option::is_none")]
        bounding_volume: Option<BoundingVolume>,
        uri: String,
    }

    /// Implicit tiling object.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#implicit-root-tile
    /// https://github.com/CesiumGS/3d-tiles/blob/1.1/specification/schema/tile.implicitTiling.schema.json
    #[derive(Serialize, Deserialize, Default, Debug, Clone)]
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
    #[derive(Serialize, Deserialize, Debug, Default, Clone)]
    #[serde(rename_all = "UPPERCASE")]
    enum SubdivisionScheme {
        #[default]
        Quadtree,
        Octree,
    }

    /// Implicit tiling subtrees.
    /// https://github.com/CesiumGS/3d-tiles/tree/1.1/specification/ImplicitTiling#subtrees
    #[derive(Serialize, Deserialize, Debug, Clone)]
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

    impl Subtrees {
        fn new(subtrees_dir: &str) -> Self {
            Self {
                uri: format!("{}/{{level}}/{{x}}/{{y}}.subtree", subtrees_dir),
            }
        }
    }

    /// Implicit tiling subtree object.
    /// Metadata is not supported.
    /// https://github.com/CesiumGS/3d-tiles/blob/1.1/specification/schema/Subtree/subtree.schema.json
    #[derive(Serialize, Deserialize, Debug, Default)]
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

    #[derive(Serialize, Deserialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct Buffer {
        #[serde(skip_serializing_if = "Option::is_none")]
        name: Option<String>,
        byte_length: usize,
    }

    #[derive(Serialize, Deserialize, Debug, Default)]
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
    #[derive(Serialize, Deserialize, Debug, Default)]
    #[serde(rename_all = "camelCase")]
    struct Availability {
        /// An integer index that identifies the buffer view containing the availability bitstream.
        #[serde(skip_serializing_if = "Option::is_none")]
        bitstream: Option<usize>,
        #[serde(skip_serializing_if = "Option::is_none")]
        available_count: Option<usize>,
        #[serde(skip_serializing_if = "Option::is_none")]
        constant: Option<AvailabilityConstant>,
    }

    /// Integer indicating whether all of the elements are Available (1) or all are Unavailable (0).
    #[allow(dead_code)]
    #[derive(Debug, Default, Serialize_repr, Deserialize_repr)]
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

            world.export_grid(false, None).unwrap();

            let quadtree = QuadTree::from_world(
                &world,
                crate::spatial_structs::QuadTreeCapacity::Vertices(15000),
            );
            quadtree.export(&world, None).unwrap();

            let tileset = Tileset::from_quadtree(&quadtree, &world, 16_f64, 200, None, None, true);

            // tileset.make_implicit(&world.grid, &quadtree, );

            let i = ImplicitTiling::default();
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
            let crs_to = "EPSG:4978";
            let transformer = Proj::new_known_crs("EPSG:7415", crs_to, None).unwrap();
            let bbox: Bbox = [171790.0, 472690.0, -15.0, 274190.0, 575090.0, 400.0];
            let bbox: Bbox = [
                84362.90299999999,
                446306.814,
                -20.66,
                212362.903,
                574306.814,
                62.882,
            ];
            let bounding_volume = BoundingVolume::box_from_bbox(&bbox, &transformer).unwrap();
            println!("{:?}", serde_json::to_string(&bounding_volume));
        }

        #[test]
        fn test_boundingvolume_region_from_bbox() {
            let crs_to = "EPSG:4979";
            let transformer = Proj::new_known_crs("EPSG:7415", crs_to, None).unwrap();
            let bbox: Bbox = [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
            let bounding_volume = BoundingVolume::region_from_bbox(&bbox, &transformer);
            println!("{:?}", bounding_volume);
        }

        #[test]
        fn test_boundingvolume_json_precision() {
            let bounding_volume = BoundingVolume::Region([
                0.0762316935076296,
                0.9075853268461559,
                0.07639433996161832,
                0.9076933035069118,
                38.18193225618456,
                96.3895906072444,
            ]);
            println!("{}", serde_json::to_string(&bounding_volume).unwrap());
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
