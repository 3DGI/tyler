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
    use std::collections::HashSet;
    use std::collections::VecDeque;
    use std::fmt::{Display, Formatter};
    use std::fs::File;
    use std::io::Write;
    use std::path::Path;

    use bitvec::prelude as bv;
    use log::{debug, error, warn};
    use morton_encoding::morton_encode;
    use serde::{Deserialize, Serialize};
    use serde_repr::{Deserialize_repr, Serialize_repr};

    use crate::coordinates::RootEnuFrame;
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

        #[allow(dead_code)]
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
            let [outdir_tileset, outdir_tileset_content] = match output_dir {
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
                let wkt = tile.bounding_volume.as_wkt();
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
                        let wkt_content_bbox = bv.as_wkt();
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
            geometric_error_factor: f64,
            _arg_cellsize: u32,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
            content_bv_from_tile: bool,
            content_add_bv: bool,
            root_enu_frame: &RootEnuFrame,
        ) -> Self {
            let crs_from = format!("EPSG:{}", world.crs.to_epsg().unwrap());
            let transformer = Proj::new_known_crs(&crs_from, "EPSG:4979", None).unwrap();

            let root = Self::generate_tiles(
                quadtree,
                world,
                &transformer,
                geometric_error_factor,
                arg_minz,
                arg_maxz,
                content_bv_from_tile,
                content_add_bv,
            );
            log::info!(
                "Root ENU frame - input CRS origin: [{:.2}, {:.2}, {:.2}], geodetic: [{:.8}, {:.8}, {:.2}], ECEF: [{:.2}, {:.2}, {:.2}]",
                root_enu_frame.source_origin[0],
                root_enu_frame.source_origin[1],
                root_enu_frame.source_origin[2],
                root_enu_frame.geodetic_origin[0],
                root_enu_frame.geodetic_origin[1],
                root_enu_frame.geodetic_origin[2],
                root_enu_frame.ecef_origin[0],
                root_enu_frame.ecef_origin[1],
                root_enu_frame.ecef_origin[2],
            );
            let root_transform = Transform(root_enu_frame.transform());

            // Using gltf tile content
            let mut extensions: Extensions = HashMap::new();
            let e1 = Extension::ContentGtlf {
                extensions_used: None,
                extensions_required: None,
            };
            extensions.insert(ExtensionName::ContentGltf, e1);

            // Apply root transform to root tile
            let mut root_with_transform = root;
            root_with_transform.transform = Some(root_transform);

            Self {
                asset: Default::default(),
                geometric_error: root_with_transform.geometric_error,
                root: root_with_transform,
                properties: None,
                extensions_used: Some(vec![ExtensionName::ContentGltf]),
                extensions_required: Some(vec![ExtensionName::ContentGltf]),
                extensions: Some(extensions),
            }
        }

        // TODO: The function outputs a Tile even if the quadtree node has 0 items,
        //  because the function is recursive and it must output a Tile. It would be
        //  more elegant to output Option<Tile>, but that needs refactoring downstream
        //  (eg. serialization).
        fn generate_tiles(
            quadtree: &QuadTree,
            world: &crate::parser::World,
            transformer: &Proj,
            geometric_error_factor: f64,
            arg_minz: Option<i32>,
            arg_maxz: Option<i32>,
            content_bv_from_tile: bool,
            content_add_bv: bool,
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
                    BoundingVolume::region_from_bbox(&tile_bbox, transformer).unwrap();

                // Tyler currently writes full-detail content only at leaves. Internal
                // tiles are traversal nodes, so their geometric error is a spatial
                // refinement trigger proportional to the tile footprint.
                let tile_width = (tile_bbox[3] - tile_bbox[0]).max(tile_bbox[4] - tile_bbox[1]);
                let geometric_error = tile_width * geometric_error_factor;
                let mut tile_children: Vec<Tile> = Vec::new();
                for child in quadtree.children.iter() {
                    tile_children.push(Self::generate_tiles(
                        child,
                        world,
                        transformer,
                        geometric_error_factor,
                        arg_minz,
                        arg_maxz,
                        content_bv_from_tile,
                        content_add_bv,
                    ));
                }
                Tile {
                    id: tile_id,
                    bounding_volume,
                    geometric_error,
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
                let bounding_volume =
                    BoundingVolume::region_from_bbox(&tile_bbox, transformer).unwrap();
                let mut content: Option<Content> = None;

                if quadtree.nr_items > 0 {
                    let mut tile_content_bbox_rw =
                        quadtree.node_content_bbox(world, arg_minz, arg_maxz);
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
                        BoundingVolume::region_from_bbox(&tile_content_bbox_rw, transformer)
                            .unwrap();

                    content = Some(Content {
                        bounding_volume: if content_add_bv {
                            Some(content_bounding_volume)
                        } else {
                            None
                        },
                        uri: format!("t/{}.glb", quadtree.id),
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

        /// Flatten the tile hierarchy, visiting each tile in the quadtree.
        /// If 'levels_up' is provided, the tiles will be flattened only
        /// 'n levels upwards from the leaves', outputting only the flattened tiles
        /// (instead of the whole tree).
        #[allow(dead_code)]
        pub fn flatten(&self, levels_up: Option<u16>) -> Vec<&Tile> {
            self.root.flatten(levels_up)
        }

        pub fn collect_leaves(&self) -> Vec<&Tile> {
            self.root.collect_leaves()
        }

        #[allow(dead_code)]
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
        #[allow(dead_code)]
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
                        let grid_for_level = SquareGrid::new(&tile_bbox, tile_width, grid_epsg);
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
                uri: "t/{level}/{x}/{y}.glb".to_string(),
            });
            self.root.refine = Some(Refinement::Add);
            self.root.children = None;
            (flat_tiles_with_content, subtrees_vec)
        }

        /// Convert the root tile to an implicit tileset using content tile IDs
        /// that already follow the implicit subdivision coordinate system.
        pub fn make_implicit_from_content_tile_ids(
            &mut self,
            content_tile_ids: &[TileId],
            subtrees_dir: Option<&str>,
        ) -> Vec<(TileId, Vec<u8>)> {
            let available_levels = self.available_levels();
            let subtree_levels = available_levels;
            let subtrees = match subtrees_dir {
                None => Subtrees::default(),
                Some(dirname) => Subtrees::new(dirname),
            };
            self.root.implicit_tiling = Some(ImplicitTiling {
                subdivision_scheme: SubdivisionScheme::Quadtree,
                subtree_levels,
                available_levels,
                subtrees,
            });
            self.root.content = Some(Content {
                bounding_volume: None,
                uri: "t/{level}/{x}/{y}.glb".to_string(),
            });
            self.root.refine = Some(Refinement::Add);
            self.root.children = None;

            vec![Self::implicit_root_subtree_from_content_tile_ids(
                available_levels,
                content_tile_ids,
            )]
        }

        fn implicit_root_subtree_from_content_tile_ids(
            available_levels: u16,
            content_tile_ids: &[TileId],
        ) -> (TileId, Vec<u8>) {
            let mut available_tiles = HashSet::new();
            let mut content_tiles = HashSet::new();

            for tile_id in content_tile_ids {
                if tile_id.level >= available_levels {
                    warn!(
                        "Skipping implicit content tile {} because availableLevels is {}",
                        tile_id, available_levels
                    );
                    continue;
                }

                content_tiles.insert(tile_id.clone());
                for ancestor_level in 0..=tile_id.level {
                    let shift = tile_id.level - ancestor_level;
                    available_tiles.insert(TileId::new(
                        tile_id.x >> shift,
                        tile_id.y >> shift,
                        ancestor_level,
                    ));
                }
            }

            let nr_tiles_total_subtree = (4_usize.pow(available_levels as u32) - 1) / 3;
            let mut tile_availability_bitstream: bv::BitVec<u8, bv::Lsb0> = bv::BitVec::new();
            tile_availability_bitstream.resize(nr_tiles_total_subtree, false);
            let mut content_availability_bitstream: bv::BitVec<u8, bv::Lsb0> = bv::BitVec::new();
            content_availability_bitstream.resize(nr_tiles_total_subtree, false);

            for tile_id in available_tiles {
                let index = Self::implicit_tile_bit_index(&tile_id);
                if index < nr_tiles_total_subtree {
                    tile_availability_bitstream.set(index, true);
                }
            }
            for tile_id in content_tiles {
                let index = Self::implicit_tile_bit_index(&tile_id);
                if index < nr_tiles_total_subtree {
                    content_availability_bitstream.set(index, true);
                }
            }

            let nr_tiles_child_level = 4_usize.pow(available_levels as u32);
            let mut child_subtree_availability_bitstream: bv::BitVec<u8, bv::Lsb0> =
                bv::BitVec::new();
            child_subtree_availability_bitstream.resize(nr_tiles_child_level, false);

            let subtree_id = TileId::new(0, 0, 0);
            (
                subtree_id,
                Self::subtree_bytes_from_availability(
                    &tile_availability_bitstream,
                    &content_availability_bitstream,
                    &child_subtree_availability_bitstream,
                ),
            )
        }

        fn implicit_tile_bit_index(tile_id: &TileId) -> usize {
            let level_offset = (4_usize.pow(tile_id.level as u32) - 1) / 3;
            let morton_index = morton_encode([tile_id.y as u64, tile_id.x as u64]) as usize;
            level_offset + morton_index
        }

        fn subtree_bytes_from_availability(
            tile_availability_bitstream: &bv::BitVec<u8, bv::Lsb0>,
            content_availability_bitstream: &bv::BitVec<u8, bv::Lsb0>,
            child_subtree_availability_bitstream: &bv::BitVec<u8, bv::Lsb0>,
        ) -> Vec<u8> {
            let mut buffer_vec: Vec<u8> = Vec::new();
            let mut bufferviews: Vec<BufferView> = Vec::with_capacity(3);
            let mut bufferview_idx: usize = 0;

            let mut tile_availability_bits = tile_availability_bitstream.clone();
            let tile_availability =
                Self::create_availability(bufferview_idx, &mut tile_availability_bits);
            if tile_availability.constant.is_none() {
                Self::add_padding(&mut buffer_vec, 8);
                Self::add_bitstream(&mut buffer_vec, &mut bufferviews, tile_availability_bits);
                bufferview_idx += 1;
            }

            let mut content_availability_bits = content_availability_bitstream.clone();
            let content_availability =
                Self::create_availability(bufferview_idx, &mut content_availability_bits);
            if content_availability.constant.is_none() {
                Self::add_padding(&mut buffer_vec, 8);
                Self::add_bitstream(&mut buffer_vec, &mut bufferviews, content_availability_bits);
                bufferview_idx += 1;
            }

            let mut child_subtree_availability_bits = child_subtree_availability_bitstream.clone();
            let child_subtree_availability =
                Self::create_availability(bufferview_idx, &mut child_subtree_availability_bits);
            if child_subtree_availability.constant.is_none() {
                Self::add_padding(&mut buffer_vec, 8);
                Self::add_bitstream(
                    &mut buffer_vec,
                    &mut bufferviews,
                    child_subtree_availability_bits,
                );
            }

            Self::add_padding(&mut buffer_vec, 8);

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
            let mut subtree_json =
                serde_json::to_string(&subtree).expect("failed to serialize the subtree to json");
            let remainder = subtree_json.as_bytes().len() % 8;
            if remainder > 0 {
                let padding = 8 - remainder;
                for _i in 0..padding {
                    subtree_json += " ";
                }
            }
            let subtree_json_bytes = subtree_json.as_bytes();

            let mut subtree_bytes: Vec<u8> = Vec::new();
            subtree_bytes.extend(0x74627573u32.to_le_bytes());
            subtree_bytes.extend(1_u32.to_le_bytes());
            subtree_bytes.extend((subtree_json_bytes.len() as u64).to_le_bytes());
            subtree_bytes.extend((buffer_vec.len() as u64).to_le_bytes());
            subtree_bytes.extend(subtree_json_bytes);
            subtree_bytes.extend(buffer_vec);
            subtree_bytes
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

        #[allow(dead_code)]
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
        #[allow(dead_code)]
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
            let grid_for_level = SquareGrid::new(bbox, tile_width, epsg);

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

        /// Splits a tileset into several tilesets at the given level, to create
        /// [external tilesets](https://docs.ogc.org/cs/22-025r4/22-025r4.html#core-external-tilesets).
        /// The tile at `level` becomes the root tile of the new tileset.
        /// It modifies the tileset, by removing the tiles after `level` and replacing
        /// the `Tile.Content.uri` of the tiles on `level` to point to the child tilesets that are
        /// returned from this function. The child tileset URIs follow the the pattern of
        /// `tileset-<root tile ID>.json`, for example `tileset-7-0-0.json`.
        /// The child tileset URI is returned with the Tileset, as the first member of the tuple.
        pub fn split(&mut self, level: u16) -> Vec<(String, Tileset)> {
            let max_nr_tilesets = 4_usize.pow(level as u32);
            let mut child_tilesets: Vec<(String, Tileset)> = Vec::with_capacity(max_nr_tilesets);
            let mut q = VecDeque::new();
            q.push_back(&mut self.root);
            while let Some(tile) = q.pop_front() {
                if tile.id.level == level {
                    let filename =
                        format!("tileset-{}-{}-{}.json", tile.id.level, tile.id.x, tile.id.y);
                    // Create the new tileset
                    child_tilesets.push((
                        filename.clone(),
                        Tileset {
                            asset: Default::default(),
                            geometric_error: tile.geometric_error,
                            root: tile.clone(),
                            properties: None,
                            extensions_used: None,
                            extensions_required: None,
                            extensions: None,
                        },
                    ));
                    // Update the current tile to point to the new tileset
                    tile.content = Some(Content {
                        bounding_volume: None,
                        uri: filename,
                    });
                    tile.children = None;
                }
                if let Some(ref mut children) = tile.children {
                    for child in children.iter_mut() {
                        q.push_back(child);
                    }
                }
            }
            child_tilesets
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
        #[allow(dead_code)]
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
        #[allow(dead_code)]
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

        #[allow(dead_code)]
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
        #[allow(dead_code)]
        pub fn add_content(&mut self) {
            self.content = Some(Content {
                bounding_volume: Some(self.bounding_volume),
                uri: format!("t/{}.glb", self.id),
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

    #[derive(
        Serialize, Deserialize, Clone, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd,
    )]
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

    /// Format the BoundingVolume coordinates with appropriate precision in the JSON output.
    /// For Region bounding volumes (lat/lon in radians), we need 6+ decimal places for precision.
    /// 2 decimal places for radians gives ~0.01 rad = ~0.57° = ~63km precision (way too coarse!).
    /// 6 decimal places gives ~0.000001 rad = ~0.000057° = ~6.3m precision (appropriate).
    /// See https://wiki.openstreetmap.org/wiki/Precision_of_coordinates
    struct BoundingVolumeFormatter;

    impl serde_json::ser::Formatter for BoundingVolumeFormatter {
        fn write_f64<W>(&mut self, writer: &mut W, value: f64) -> std::io::Result<()>
        where
            W: ?Sized + Write,
        {
            // Use 6 decimal places for sufficient precision in Region bounding volumes (radians)
            // This ensures lat/lon coordinates have ~6.3m precision instead of ~63km
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
        /// The output is an array of 12 numbers that define an oriented bounding box in a
        /// right-handed 3-axis (x, y, z) Cartesian coordinate system where the z-axis is up.
        /// The first three elements define the x, y, and z values for the center of the box.
        /// The next three elements (with indices 3, 4, and 5) define the x-axis direction
        /// and half-length.
        /// The next three elements (indices 6, 7, and 8) define the y-axis direction and
        /// half-length.
        /// The last three elements (indices 9, 10, and 11) define the z-axis direction and
        /// half-length.
        #[allow(dead_code)]
        fn box_from_bbox(
            bbox: &Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            // Input CRS box dimensions and center
            let dx = bbox[3] - bbox[0];
            let dy = bbox[4] - bbox[1];
            let dz = bbox[5] - bbox[2];
            let center: [f64; 3] = [bbox[0] + dx * 0.5, bbox[1] + dy * 0.5, (bbox[2] + dz * 0.5)];

            // The center points on the side faces in the X and Y directions
            let x_max_pt: [f64; 3] = [bbox[3], center[1], center[2]];
            let y_max_pt: [f64; 3] = [center[0], bbox[4], center[2]];
            // let z_max_pt: [f64; 3] = [center[0], center[1], bbox[5] + magic_z_shift];

            // Determine X/Y axis orientation in target CRS (ECEF). We transform both endpoints
            // of a short vector in the X/Y direction in the input CRS. The lenght of this vector
            // should be short to minimize distortions caused by the CRS conversion. The distortion
            // can cause a tilt in the ECEF box. Assuming we have a input CRS unit of metres,
            // 1 meter is a good choice for the input vector length.
            // NB: be careful if input CRS is not in metres.
            let center_ecef = transformer.convert((center[0], center[1], center[2]))?;
            let pnx = transformer.convert((center[0] + 1.0, center[1], center[2]))?;
            let pny = transformer.convert((center[0], center[1] + 1.0, center[2]))?;

            // Compute the correct half lengths for X and Y vectors of the box
            let x_max_pt_ecef = transformer.convert((x_max_pt[0], x_max_pt[1], x_max_pt[2]))?;
            let y_max_pt_ecef = transformer.convert((y_max_pt[0], y_max_pt[1], y_max_pt[2]))?;

            let dx_ = ((x_max_pt_ecef.0 - center_ecef.0).powi(2)
                + (x_max_pt_ecef.1 - center_ecef.1).powi(2)
                + (x_max_pt_ecef.2 - center_ecef.2).powi(2))
            .sqrt();
            let dy_ = ((y_max_pt_ecef.0 - center_ecef.0).powi(2)
                + (y_max_pt_ecef.1 - center_ecef.1).powi(2)
                + (y_max_pt_ecef.2 - center_ecef.2).powi(2))
            .sqrt();

            // Compute half length X Y vectors of the ECEF box
            let s_to_unit_vx = ((pnx.0 - center_ecef.0).powi(2)
                + (pnx.1 - center_ecef.1).powi(2)
                + (pnx.2 - center_ecef.2).powi(2))
            .sqrt();
            let vx = (
                (pnx.0 - center_ecef.0) / s_to_unit_vx * dx_,
                (pnx.1 - center_ecef.1) / s_to_unit_vx * dx_,
                (pnx.2 - center_ecef.2) / s_to_unit_vx * dx_,
            );
            let s_to_unit_vy = ((pny.0 - center_ecef.0).powi(2)
                + (pny.1 - center_ecef.1).powi(2)
                + (pny.2 - center_ecef.2).powi(2))
            .sqrt();
            let vy = (
                (pny.0 - center_ecef.0) / s_to_unit_vy * dy_,
                (pny.1 - center_ecef.1) / s_to_unit_vy * dy_,
                (pny.2 - center_ecef.2) / s_to_unit_vy * dy_,
            );

            // Z unit vector in the ECEF box (before curvature correction)
            let dvz =
                (center_ecef.0.powi(2) + center_ecef.1.powi(2) + center_ecef.2.powi(2)).sqrt();
            let vz_unit = (
                center_ecef.0 / dvz,
                center_ecef.1 / dvz,
                center_ecef.2 / dvz,
            );

            // Calculate the height difference between the lower corners and the curved earth surface
            let r_earth: f64 = 6371000.0;
            let dxy_ = (dx_.powi(2) + dy.powi(2)).sqrt();
            let curvature_drop = (dxy_.powi(2) + r_earth.powi(2)).sqrt() - r_earth; // this is the h difference between the lower corners and the earth surface

            // Calculate the total correction that needs to be applied to the center point of the ECEF box
            let center_z_correction = (dz + curvature_drop) / 2.0 - 0.5 * dz;

            // Drop center_ecef
            let center_ecef_dropped = (
                center_ecef.0 - center_z_correction * vz_unit.0,
                center_ecef.1 - center_z_correction * vz_unit.1,
                center_ecef.2 - center_z_correction * vz_unit.2,
            );

            // Final box matrix, NB the Z half length vector are now also corrected for earth curvature
            Ok(Self::Box([
                center_ecef_dropped.0,
                center_ecef_dropped.1,
                center_ecef_dropped.2,
                vx.0,
                vx.1,
                vx.2,
                vy.0,
                vy.1,
                vy.2,
                vz_unit.0 * (curvature_drop + dz) / 2.0,
                vz_unit.1 * (curvature_drop + dz) / 2.0,
                vz_unit.2 * (curvature_drop + dz) / 2.0,
            ]))
        }

        #[allow(dead_code)]
        fn region_from_bbox(
            bbox: &Bbox,
            transformer: &Proj,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            let mut west = f64::INFINITY;
            let mut south = f64::INFINITY;
            let mut min_h = f64::INFINITY;
            let mut east = f64::NEG_INFINITY;
            let mut north = f64::NEG_INFINITY;
            let mut max_h = f64::NEG_INFINITY;

            for [x, y, z] in bbox_corners(bbox) {
                let (lon, lat, height) = transformer.convert((x, y, z))?;
                west = west.min(lon);
                south = south.min(lat);
                min_h = min_h.min(height);
                east = east.max(lon);
                north = north.max(lat);
                max_h = max_h.max(height);
            }

            Ok(BoundingVolume::Region([
                west.to_radians(),
                south.to_radians(),
                east.to_radians(),
                north.to_radians(),
                min_h,
                max_h,
            ]))
        }

        /// Cast to 2D WKT
        fn as_wkt(&self) -> String {
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

    fn bbox_corners(bbox: &Bbox) -> [[f64; 3]; 8] {
        [
            [bbox[0], bbox[1], bbox[2]],
            [bbox[0], bbox[1], bbox[5]],
            [bbox[0], bbox[4], bbox[2]],
            [bbox[0], bbox[4], bbox[5]],
            [bbox[3], bbox[1], bbox[2]],
            [bbox[3], bbox[1], bbox[5]],
            [bbox[3], bbox[4], bbox[2]],
            [bbox[3], bbox[4], bbox[5]],
        ]
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
        use std::fs;
        use std::path::PathBuf;
        use std::time::{SystemTime, UNIX_EPOCH};

        fn test_data_dir() -> PathBuf {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("resources")
                .join("data")
        }

        fn resource_path(name: &str) -> PathBuf {
            test_data_dir().join(name)
        }

        fn unique_test_dir(prefix: &str) -> PathBuf {
            let unique = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system time")
                .as_nanos();
            let path = std::env::temp_dir().join(format!("tyler-{prefix}-{unique}"));
            fs::create_dir_all(&path).expect("create test dir");
            path
        }

        fn assert_tile_bounding_volumes_are_regions(tile: &Tile) {
            assert!(
                matches!(tile.bounding_volume, BoundingVolume::Region(_)),
                "tile {} bounding volume should be a region",
                tile.id
            );
            if let Some(content) = &tile.content {
                if let Some(bounding_volume) = &content.bounding_volume {
                    assert!(
                        matches!(bounding_volume, BoundingVolume::Region(_)),
                        "tile {} content bounding volume should be a region",
                        tile.id
                    );
                }
            }
            if let Some(children) = &tile.children {
                for child in children {
                    assert_tile_bounding_volumes_are_regions(child);
                }
            }
        }

        fn assert_geometric_error_decreases_to_zero(tile: &Tile) {
            if let Some(children) = &tile.children {
                for child in children {
                    assert!(
                        tile.geometric_error >= child.geometric_error,
                        "tile {} geometricError {} should be >= child {} geometricError {}",
                        tile.id,
                        tile.geometric_error,
                        child.id,
                        child.geometric_error
                    );
                    assert_geometric_error_decreases_to_zero(child);
                }
            } else {
                assert!(
                    tile.geometric_error.abs() <= f64::EPSILON,
                    "leaf tile {} should have zero geometricError",
                    tile.id
                );
            }
        }

        #[test]
        fn test_implicittiling() {
            // 85162.9 447106.8 85562.9 447706.8
            // let bbox: crate::spatial_structs::Bbox =
            //     [85162.9, 447106.8, -10.7, 85962.9, 447906.8, 320.5];
            // let grid = crate::spatial_structs::SquareGrid::new(&bbox, 200, 7415, Some(10.0));
            let dataset_dir = unique_test_dir("implicittiling");
            let features_dir = dataset_dir.join("features");
            fs::create_dir_all(&features_dir).expect("create features dir");
            let metadata_path = dataset_dir.join("metadata.city.json");
            fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
            fs::copy(
                resource_path("3dbag_feature_x71.city.jsonl"),
                features_dir.join("sample.city.jsonl"),
            )
            .expect("copy feature");

            let mut world = crate::parser::World::new(
                &metadata_path,
                &features_dir,
                200,
                Some(vec![
                    crate::parser::CityObjectType::Building,
                    crate::parser::CityObjectType::BuildingPart,
                ]),
                None,
                None,
            )
            .unwrap();
            world.index_with_grid().unwrap();

            world.export_grid(false, None).unwrap();

            let quadtree = QuadTree::from_world(
                &world,
                crate::spatial_structs::QuadTreeCapacity::Vertices(15000),
            );
            quadtree.export(&world, None).unwrap();

            let source_crs = format!("EPSG:{}", world.crs.to_epsg().unwrap());
            let root_bbox = quadtree.bbox(&world.grid);
            let root_enu_frame =
                RootEnuFrame::from_bbox(&source_crs, &root_bbox).expect("root ENU frame");
            let tileset = Tileset::from_quadtree(
                &quadtree,
                &world,
                0.024,
                200,
                None,
                None,
                true,
                true,
                &root_enu_frame,
            );
            assert_tile_bounding_volumes_are_regions(&tileset.root);
            assert_geometric_error_decreases_to_zero(&tileset.root);
            assert!((tileset.geometric_error - tileset.root.geometric_error).abs() <= f64::EPSILON);
            let root_transform = tileset
                .root
                .transform
                .as_ref()
                .expect("root tile should have ENU-to-ECEF transform");
            assert!((root_transform.0[0] - 1.0).abs() > 1.0e-6);
            assert!(root_transform.0[12].abs() > 1.0);

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
            let bounding_volume = BoundingVolume::region_from_bbox(&bbox, &transformer)
                .expect("bbox should transform to EPSG:4979 region");
            let BoundingVolume::Region(region) = bounding_volume else {
                panic!("region_from_bbox should return a region");
            };

            let mut expected = [
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ];
            for [x, y, z] in bbox_corners(&bbox) {
                let (lon, lat, height) = transformer.convert((x, y, z)).unwrap();
                expected[0] = expected[0].min(lon.to_radians());
                expected[1] = expected[1].min(lat.to_radians());
                expected[2] = expected[2].max(lon.to_radians());
                expected[3] = expected[3].max(lat.to_radians());
                expected[4] = expected[4].min(height);
                expected[5] = expected[5].max(height);
            }

            for (actual, expected) in region.into_iter().zip(expected) {
                assert!(
                    (actual - expected).abs() < 1.0e-12,
                    "expected {expected}, got {actual}"
                );
            }
            assert!(region[0].abs() < std::f64::consts::PI);
            assert!(region[1].abs() < std::f64::consts::FRAC_PI_2);
            assert!(region[0] < region[2]);
            assert!(region[1] < region[3]);
            assert!(region[4] <= region[5]);
        }

        #[test]
        fn test_bbox_corners_returns_all_eight_corners() {
            let bbox: Bbox = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
            let corners = bbox_corners(&bbox);
            assert_eq!(corners.len(), 8);
            assert!(corners.contains(&[1.0, 2.0, 3.0]));
            assert!(corners.contains(&[1.0, 2.0, 6.0]));
            assert!(corners.contains(&[1.0, 5.0, 3.0]));
            assert!(corners.contains(&[1.0, 5.0, 6.0]));
            assert!(corners.contains(&[4.0, 2.0, 3.0]));
            assert!(corners.contains(&[4.0, 2.0, 6.0]));
            assert!(corners.contains(&[4.0, 5.0, 3.0]));
            assert!(corners.contains(&[4.0, 5.0, 6.0]));
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
