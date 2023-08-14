//! Spatial data structures for indexing the features.
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
use crate::parser::FeatureSet;
use log::{debug, warn};
use std::collections::VecDeque;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use morton_encoding::{morton_decode, morton_encode};

/// Quadtree
///
/// We don't expect that the quadtree has more than 65535 levels (u16).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QuadTree {
    pub id: QuadTreeNodeId,
    side_length: u64,
    pub children: Vec<QuadTree>,
    cells: Vec<CellId>,
    pub nr_items: usize,
}

impl QuadTree {
    pub fn from_world(world: &crate::parser::World, limit: QuadTreeCapacity) -> Self {
        Self::from_grid(&world.grid, limit)
    }

    fn from_grid(grid: &SquareGrid, limit: QuadTreeCapacity) -> Self {
        let mut merge_limit: usize = 0;
        let nr_cells = grid.length.pow(2) as f64;
        let max_level = (nr_cells.ln() / 4.0_f64.ln()).ceil() as u16;
        debug!("Calculated maximum level for quadtree: {}", &max_level);
        // Somehow, the morton_encoding::morton_encode needs [y,x] for creating a Z-curve.
        //  If [x,y] is used, it creates an N-curve.
        //  crate::spatial_structs::interleave(x,y) creates a Z-curve.
        let mut mortoncodes: Vec<u128> = grid
            .into_iter()
            .map(|(cellid, _)| morton_encode([cellid.row as u64, cellid.column as u64]))
            .collect();
        mortoncodes.sort();
        let tiles_morton: Vec<QuadTree> = mortoncodes
            .iter()
            .map(|mc| {
                let coords: [u64; 2] = morton_decode(*mc);
                coords
            })
            .map(|[x, y]| {
                let cellid = CellId {
                    row: y as usize,
                    column: x as usize,
                };
                let items: usize;
                match limit {
                    QuadTreeCapacity::Objects(l) => {
                        // Use the number of features as a limit
                        items = grid.cell(&cellid).feature_ids.len();
                        merge_limit = l;
                    }
                    QuadTreeCapacity::Vertices(l) => {
                        // Use the number of vertices as a limit
                        items = grid.cell(&cellid).nr_vertices;
                        merge_limit = l;
                    }
                }
                QuadTree {
                    id: QuadTreeNodeId::new(x as usize, y as usize, max_level),
                    side_length: grid.cellsize as u64,
                    children: Vec::new(),
                    cells: vec![cellid],
                    nr_items: items,
                }
            })
            .collect();
        Self::merge_tiles(0, tiles_morton, merge_limit)
    }

    fn merge_tiles(level: u16, tiles: Vec<QuadTree>, limit: usize) -> QuadTree {
        let len_tiles = tiles.len();
        if len_tiles > 4 {
            let q0: usize = len_tiles / 4;
            let q1: usize = q0 * 2;
            let q2: usize = q0 * 3;
            let next_level = level + 1;
            Self::merge_tiles(
                level,
                vec![
                    Self::merge_tiles(next_level, tiles[0..q0].to_vec(), limit),
                    Self::merge_tiles(next_level, tiles[q0..q1].to_vec(), limit),
                    Self::merge_tiles(next_level, tiles[q1..q2].to_vec(), limit),
                    Self::merge_tiles(next_level, tiles[q2..].to_vec(), limit),
                ],
                limit,
            )
        } else {
            let sum_items: usize = tiles.iter().map(|t| t.nr_items).sum();
            let mut cells: Vec<CellId> = Vec::new();
            for t in tiles.iter() {
                for c in t.cells.iter() {
                    cells.push(*c)
                }
            }
            let id = QuadTreeNodeId::new(tiles[0].id.x, tiles[0].id.y, level);
            let _id_string = id.to_string();
            // FIXME: this also adds the quadtree if sum_items == 0 so the parent will have 4
            //  children instead of 3. Probably should return Option<Quadtree>.
            //  Currently these empty tiles are removed in Tile.prune().
            if sum_items <= limit {
                QuadTree {
                    id,
                    side_length: tiles[0].side_length * 2,
                    children: vec![],
                    cells,
                    nr_items: sum_items,
                }
            } else {
                if tiles.len() % 4 != 0 {
                    warn!(
                        "number of children is not dividable by 4: {}, in level {}",
                        tiles.len(),
                        &level
                    );
                }
                QuadTree {
                    id,
                    side_length: tiles[0].side_length * 2,
                    children: tiles.clone(),
                    cells: vec![],
                    nr_items: sum_items,
                }
            }
        }
    }

    #[allow(dead_code)]
    fn collect_leaves_recurse<'collect>(&'collect self, leaves: &mut Vec<&'collect QuadTree>) {
        if !self.children.is_empty() {
            for child in self.children.iter() {
                child.collect_leaves_recurse(leaves);
            }
        } else {
            leaves.push(self);
        }
    }

    #[allow(dead_code)]
    pub fn collect_leaves(&self) -> Vec<&Self> {
        let mut leaves: Vec<&QuadTree> = Vec::new();
        self.collect_leaves_recurse(&mut leaves);
        leaves
    }

    pub fn bbox(&self, grid: &SquareGrid) -> Bbox {
        let minx = grid.origin[0] + (self.id.x * grid.cellsize as usize) as f64;
        let miny = grid.origin[1] + (self.id.y * grid.cellsize as usize) as f64;
        [
            minx,
            miny,
            grid.bbox[2],
            minx + self.side_length as f64,
            miny + self.side_length as f64,
            grid.bbox[5],
        ]
    }

    /// Compute the bounding box of all the features in the node
    pub fn node_content_bbox(
        &self,
        world: &crate::parser::World,
        arg_minz: Option<i32>,
        arg_maxz: Option<i32>,
    ) -> Bbox {
        let mut tile_content_bbox_qc = BboxQc([0, 0, 0, 0, 0, 0]);
        for cellid in self.cells() {
            let cell = world.grid.cell(cellid);
            if !cell.feature_ids.is_empty() {
                tile_content_bbox_qc = world.features[cell.feature_ids[0]].bbox_qc.clone();
                break;
            }
        }
        for cellid in self.cells() {
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
        tile_content_bbox_qc.to_bbox(&world.transform, arg_minz, arg_maxz)
    }

    /// Breadth-first search for a node.
    pub fn node(&self, id: &QuadTreeNodeId) -> Option<&QuadTree> {
        let mut q = VecDeque::new();
        q.push_back(self);

        while let Some(n) = q.pop_front() {
            if &n.id == id {
                return Some(n);
            } else {
                for child in &n.children {
                    q.push_back(child);
                }
            }
        }
        // Did not find the node
        None
    }

    pub fn cells(&self) -> Vec<&CellId> {
        let mut cellids: Vec<&CellId> = Vec::new();
        let mut q = VecDeque::new();
        q.push_back(self);
        while let Some(node) = q.pop_front() {
            if node.children.is_empty() {
                for cell in node.cells.iter() {
                    cellids.push(cell);
                }
            } else {
                for child in &node.children {
                    q.push_back(child);
                }
            }
        }
        cellids
    }

    /// The WKT of the 2D boundary of the current node.
    pub fn to_wkt(&self, grid: &SquareGrid) -> String {
        bbox_to_wkt(&self.bbox(grid))
    }

    pub fn export(
        &self,
        world: &crate::parser::World,
        output_dir: Option<&Path>,
    ) -> std::io::Result<()> {
        let mut q = VecDeque::new();
        q.push_back(self);
        let mut quadtree_level: u16 = self.id.level;
        let [outdir_quadtree, outdir_quadtree_content] = match output_dir {
            None => [Path::new(""), Path::new("")],
            Some(outdir) => [outdir, outdir],
        };
        let mut file_quadtree =
            File::create(outdir_quadtree.join(format!("quadtree_level-{quadtree_level}.tsv")))?;
        let mut file_quadtree_content = File::create(
            outdir_quadtree_content.join(format!("quadtree_content_level-{quadtree_level}.tsv")),
        )?;

        file_quadtree
            .write_all("node_id\tnode_level\tnr_items\twkt\n".as_bytes())
            .expect("cannot write quadtree header");
        file_quadtree_content
            .write_all("node_id\tnode_level\tnr_items\twkt\n".as_bytes())
            .expect("cannot write quadtree content header");

        while let Some(node) = q.pop_front() {
            if node.id.level != quadtree_level {
                quadtree_level = node.id.level;
                file_quadtree = File::create(
                    outdir_quadtree.join(format!("quadtree_level-{quadtree_level}.tsv")),
                )?;
                file_quadtree_content = File::create(
                    outdir_quadtree_content
                        .join(format!("quadtree_content_level-{quadtree_level}.tsv")),
                )?;
                file_quadtree
                    .write_all("node_id\tnode_level\tnr_items\twkt\n".as_bytes())
                    .expect("cannot write quadtree header");
                file_quadtree_content
                    .write_all("node_id\tnode_level\tnr_items\twkt\n".as_bytes())
                    .expect("cannot write quadtree content header");
            }
            let wkt = node.to_wkt(&world.grid);
            file_quadtree
                .write_all(
                    format!(
                        "{}\t{}\t{}\t{}\n",
                        node.id, node.id.level, node.nr_items, wkt
                    )
                    .as_bytes(),
                )
                .expect("cannot write quadtree node");
            let node_content_bbox = node.node_content_bbox(world, None, None);
            let wkt_content_bbox = bbox_to_wkt(&node_content_bbox);
            file_quadtree_content
                .write_all(
                    format!(
                        "{}\t{}\t{}\t{}\n",
                        node.id, node.id.level, node.nr_items, wkt_content_bbox
                    )
                    .as_bytes(),
                )
                .expect("cannot write quadtree node content");

            for child in &node.children {
                q.push_back(child);
            }
        }
        Ok(())
    }

    pub fn export_bincode(
        &self,
        name: Option<&str>,
        output_dir: Option<&Path>,
    ) -> bincode::Result<()> {
        let file_name: &str = name.unwrap_or("quadtree");
        let file = match output_dir {
            None => File::create(format!("{file_name}.bincode"))?,
            Some(outdir) => File::create(outdir.join(format!("{file_name}.bincode")))?,
        };
        bincode::serialize_into(file, self)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct QuadTreeNodeId {
    pub x: usize,
    pub y: usize,
    pub level: u16,
}

impl QuadTreeNodeId {
    pub fn new(x: usize, y: usize, level: u16) -> Self {
        Self { x, y, level }
    }
}

impl Display for QuadTreeNodeId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}/{}", self.level, self.x, self.y)
    }
}

/// We have these double enum, QuadTreeCapacity and QuadTreeCapacityType, because of
/// how the CLI arguments are parsed. In the quadtree, we need QuadTreeCapacity, because
/// it can hold both the leaf capacity and the capacity type. But clap can only parse
/// into unit variants (I think), so we take the the capacity and capacity type as
/// separate arguments.
#[derive(Debug)]
pub enum QuadTreeCapacity {
    Objects(usize),
    Vertices(usize),
}

/// The type of items to count for the quadtree leaf capacity.
#[derive(Debug, Default, Clone, clap::ValueEnum)]
pub enum QuadTreeCriteria {
    Objects,
    #[default]
    Vertices,
}

/// 64-bit mask
#[allow(dead_code)]
fn part1by1_64(number: &u64) -> u64 {
    let mut n = *number;
    n &= 0x00000000ffffffff; // binary: 11111111111111111111111111111111,                                len: 32
    n = (n | (n << 16)) & 0x0000FFFF0000FFFF; // binary: 1111111111111111000000001111111111111111,                        len: 40
    n = (n | (n << 8)) & 0x00FF00FF00FF00FF; // binary: 11111111000000001111111100000000111111110000000011111111,        len: 56
    n = (n | (n << 4)) & 0x0F0F0F0F0F0F0F0F; // binary: 111100001111000011110000111100001111000011110000111100001111,    len: 60
    n = (n | (n << 2)) & 0x3333333333333333; // binary: 11001100110011001100110011001100110011001100110011001100110011,  len: 62
    n = (n | (n << 1)) & 0x5555555555555555; // binary: 101010101010101010101010101010101010101010101010101010101010101, len: 63
    n
}

/// Computing Morton-code from 64bit integers.
///
/// Reference: https://github.com/trevorprater/pymorton
#[allow(dead_code)]
pub fn interleave(x: &u64, y: &u64) -> u64 {
    part1by1_64(x) | (part1by1_64(y) << 1)
}

#[allow(dead_code)]
fn unpart1by1_64(mortoncode: &u64) -> u64 {
    let mut n = *mortoncode;
    n &= 0x5555555555555555; // binary: 101010101010101010101010101010101010101010101010101010101010101, len: 63
    n = (n ^ (n >> 1)) & 0x3333333333333333; // binary: 11001100110011001100110011001100110011001100110011001100110011,  len: 62
    n = (n ^ (n >> 2)) & 0x0f0f0f0f0f0f0f0f; // binary: 111100001111000011110000111100001111000011110000111100001111,    len: 60
    n = (n ^ (n >> 4)) & 0x00ff00ff00ff00ff; // binary: 11111111000000001111111100000000111111110000000011111111,        len: 56
    n = (n ^ (n >> 8)) & 0x0000ffff0000ffff; // binary: 1111111111111111000000001111111111111111,                        len: 40
    n = (n ^ (n >> 16)) & 0x00000000ffffffff; // binary: 11111111111111111111111111111111,                                len: 32
    n
}

/// Computing `[x, y]` from a Morton-code.
///
/// Reference: https://github.com/trevorprater/pymorton
#[allow(dead_code)]
pub fn deinterleave(mortoncode: &u64) -> [u64; 2] {
    [
        unpart1by1_64(mortoncode),
        unpart1by1_64(&(*mortoncode >> 1)),
    ]
}

/// Represents a square grid with square cells.
/// The grid stores the feature-indices in its cells.
/// The `length` of the grid is the number of cells of one dimension, thus the total
/// number of cells is obtained by `length * length`.
///
/// Note the a 'column' in the grid is represented by the X-axis, and a 'row' by the Y-axis.
/// See [CellId] for details.
///
/// The EPSG code of the input features is also stored in the grid, because the grid is initialized
/// directly from the feature coordinates without reprojection. Often we need to reproject the grid
/// to another CRS, for instance in order to convert it to 3D Tiles.
///
/// ```shell
///  (column)     (column)
///   +----+       +----+
///   |    | +---+ |    | +------------------+
///   |  --+-+-> | |    | |Vec<usize> (cell) |
///   |    | +---+ |    | +------------------+
///   |    |       |    |
///   |    | +---+ |    | +------------------+
///   |  --+-+-> | |    | |Vec<usize> (cell) |
///   |    | +---+ |    | +------------------+
///   |    |       |    |
///   |    | +---+ |    | +------------------+
///   |  --+-+-> | |    | |Vec<usize> (cell) |
///   | ^  | +---+ |    | +------------------+
///   +-+--+       +----+
///     |
/// +---+------------------------+
/// |   |                        |
/// | Vec<Vec<Vec<usize>>> (row) |
/// +----------------------------+
///
/// (created with https://asciiflow.com)
/// ```
///
/// ## Examples
///
/// ```
/// let grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 4.0], 1);
/// let grid_idx = grid.locate_point(&[2.5, 1.5]);
/// assert_eq!(grid_idx, [3_u64, 2_u64]);
/// ```
///
#[derive(Debug, Serialize, Deserialize)]
pub struct SquareGrid {
    origin: [f64; 3],
    pub bbox: Bbox,
    pub length: usize,
    cellsize: u32,
    pub data: Vec<Vec<Cell>>,
    pub epsg: u16,
}

impl Display for SquareGrid {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SquareGrid (origin: {:?}, bbox: {:?}, length: {}, cellsize: {}, data: not-displayed, epsg: {})",
            self.origin, self.bbox, self.length, self.cellsize, self.epsg
        )
    }
}

impl SquareGrid {
    /// Creates a grid with cells of `cellsize`, that covers the `extent`.
    /// The grid and the cells are square.
    /// The grid center is the `extent` center.
    /// The grid is returned as an origin coordinate and the number of cells.
    pub fn new(extent: &Bbox, cellsize: u32, epsg: u16) -> Self {
        // We only compute 2D. Z is constant for the grid.
        let extent_center = [
            extent[0] + (extent[3] - extent[0]) / 2.0,
            extent[1] + (extent[4] - extent[1]) / 2.0,
        ];
        let dx = extent[3] - extent[0];
        let dy = extent[4] - extent[1];
        // The grid starting dimension is the longest edge of the rectangle
        let mut d = if dx > dy { dx } else { dy };
        // The grid is constructed as a square, from the center of the extent.
        // We don't need to add any buffer to the extent, because we already round up when we
        // compute the number of cells.
        // We need a grid that is has 2^n cells in one dimension, so that we can build
        // a 4^n cells quadtree.
        // Adjust the cellsize so that we can get a tightly fit square on the extent
        let mut cellsize_new = d;
        let cellsize_f64 = cellsize as f64;
        loop {
            let cn = cellsize_new / 2.0;
            if cn < cellsize_f64 {
                break;
            } else {
                cellsize_new /= 2.0;
            }
        }
        let d_cells = (d / cellsize_new).ceil() as usize;
        let cellsize = cellsize_new.ceil() as u32;
        debug!(
            "Computed grid cells dimension: {}, with cell size: {}",
            &d_cells, cellsize
        );
        // Compute new dimension from the calculated length
        d = d_cells as f64 * cellsize as f64;
        let origin = [
            extent_center[0] - d / 2.0,
            extent_center[1] - d / 2.0,
            extent[2],
        ];
        let bbox = [
            origin[0],
            origin[1],
            origin[2],
            origin[0] + d,
            origin[1] + d,
            extent[5],
        ];
        // A row-vector (x-axis) to store the column-vectors (y-axis).
        let mut row: Vec<Vec<Cell>> = Vec::with_capacity(d_cells);
        // For each column create a column vector that stores the cells and for each row in the
        // column create a cell to store the feature IDs.
        row.resize_with(d_cells, || {
            let mut column: Vec<Cell> = Vec::with_capacity(d_cells);
            column.resize(
                d_cells,
                Cell {
                    feature_ids: Vec::new(),
                    nr_vertices: 0,
                },
            );
            column
        });
        Self {
            origin,
            bbox,
            length: d_cells,
            cellsize,
            data: row,
            epsg,
        }
    }

    /// Returns the cell index (x, y) where the point is located.
    pub fn locate_point(&self, point: &[f64; 2]) -> CellId {
        let dx = point[0] - self.origin[0];
        let dy = point[1] - self.origin[1];
        let col_i = (dx / self.cellsize as f64).floor() as usize;
        let row_i = (dy / self.cellsize as f64).floor() as usize;
        CellId {
            row: row_i,
            column: col_i,
        }
    }

    /// Insert a point into the grid, by assigning it to the Cell where it is located.
    pub fn insert(&mut self, point: &[f64; 2], feature_id: usize) -> CellId {
        let cell_id = self.locate_point(point);
        let cell = self.cell_mut(&cell_id);
        cell.feature_ids.push(feature_id);
        cell_id
    }

    /// Return the Cells that intersect the Bounding Box.
    pub fn intersect_bbox(&self, bbox: &Bbox) -> Vec<CellId> {
        let mut cellids: Vec<CellId> = Vec::new();
        let [minx, miny, _, maxx, maxy, _] = *bbox;
        let min_cellid = self.locate_point(&[minx, miny]);
        let max_cellid = self.locate_point(&[maxx, maxy]);
        for column in min_cellid.column..=max_cellid.column {
            for row in min_cellid.row..=max_cellid.row {
                cellids.push(CellId { row, column });
            }
        }
        cellids
    }

    /// Exports the grid and the feature centroids into TSV files into the working directory.
    /// The grid is written to `grid.tsv`.
    /// If `feature_set` is provided, then the feature centroids are written to
    /// `features.tsv`.
    /// If `feature_set` is provided, `transform` must be provided too (and vica-versa).
    /// If `output_dir` is provided, the files are written there. Else they are written to the
    /// working directory.
    pub fn export(
        &self,
        feature_set: Option<&FeatureSet>,
        transform: Option<&crate::parser::Transform>,
        output_dir: Option<&Path>,
    ) -> std::io::Result<()> {
        let [file_grid_path, file_features_path] = match output_dir {
            None => [PathBuf::from("grid.tsv"), PathBuf::from("features.tsv")],
            Some(outdir) => [outdir.join("grid.tsv"), outdir.join("features.tsv")],
        };
        let mut file_grid = File::create(&file_grid_path)?;
        let mut file_features = File::create(&file_features_path)?;
        let root_wkt = format!(
            "POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
            minx = self.bbox[0],
            miny = self.bbox[1],
            maxx = self.bbox[3],
            maxy = self.bbox[4]
        );
        file_grid
            .write_all("cell_id\tnr_items\twkt\n".as_bytes())
            .expect("cannot write grid header");
        file_grid
            .write_all(format!("x-x\t0\t{}\n", root_wkt).as_bytes())
            .expect("cannot write grid line");
        file_features
            .write_all("fid\tcell_id\twkt\n".as_bytes())
            .expect("cannot write features header");
        for (cellid, cell) in self {
            let wkt = self.cell_to_wkt(&cellid);
            file_grid
                .write_all(format!("{}\t{}\t{}\n", &cellid, cell.nr_vertices, wkt).as_bytes())
                .expect("cannot write grid line");
            let mut cellbuffer = String::new();
            if let Some(fset) = feature_set {
                for fid in cell.feature_ids.iter() {
                    let f = &fset[*fid];
                    let centroid = f.centroid(transform.unwrap());
                    cellbuffer += format!(
                        "{}\t{}\tPOINT({} {})\n",
                        fid, &cellid, centroid[0], centroid[1]
                    )
                    .as_str();
                }
            }
            if feature_set.is_some() {
                file_features
                    .write_all(cellbuffer.as_bytes())
                    .expect("cannot write cell contents");
            }
        }
        if feature_set.is_none() {
            // Remove empty file
            std::fs::remove_file(file_features_path)?;
        }
        Ok(())
    }

    pub fn cell_to_wkt(&self, cellid: &CellId) -> String {
        let minx = self.origin[0] + (cellid.column * self.cellsize as usize) as f64;
        let miny = self.origin[1] + (cellid.row * self.cellsize as usize) as f64;
        bbox_to_wkt(&[
            minx,
            miny,
            self.bbox[2],
            minx + self.cellsize as f64,
            miny + self.cellsize as f64,
            self.bbox[5],
        ])
    }

    pub fn cell_bbox(&self, cellid: &CellId) -> Bbox {
        let minx = self.origin[0] + (cellid.column * self.cellsize as usize) as f64;
        let miny = self.origin[1] + (cellid.row * self.cellsize as usize) as f64;
        let minz = self.bbox[2];
        let maxx = minx + self.cellsize as f64;
        let maxy = miny + self.cellsize as f64;
        let maxz = self.bbox[5];
        [minx, miny, minz, maxx, maxy, maxz]
    }

    pub fn cell(&self, cell_id: &CellId) -> &Cell {
        &self.data[cell_id.column][cell_id.row]
    }

    pub fn cell_mut(&mut self, cell_id: &CellId) -> &mut Cell {
        &mut self.data[cell_id.column][cell_id.row]
    }
}

/// Returns a tuple of `(CellId, &Cell)` for each cell in column-major order.
impl<'squaregrid> IntoIterator for &'squaregrid SquareGrid {
    type Item = (CellId, &'squaregrid Cell);
    type IntoIter = SquareGridIterator<'squaregrid>;

    fn into_iter(self) -> Self::IntoIter {
        SquareGridIterator {
            row_index: 0,
            col_index: 0,
            items: &self.data,
        }
    }
}

pub struct SquareGridIterator<'squaregrid> {
    row_index: usize,
    col_index: usize,
    items: &'squaregrid Vec<Vec<Cell>>,
}

impl<'squaregrid> Iterator for SquareGridIterator<'squaregrid> {
    type Item = (CellId, &'squaregrid Cell);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(column) = self.items.get(self.col_index) {
            if let Some(cell) = column.get(self.row_index) {
                let item = Some((
                    CellId {
                        row: self.row_index,
                        column: self.col_index,
                    },
                    cell,
                ));
                self.row_index += 1;
                item
            } else {
                // We are at the end of the current column, so jump to the next
                self.col_index += 1;
                self.row_index = 0;
                self.next()
            }
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Ord, PartialOrd, Eq, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    pub feature_ids: Vec<usize>,
    pub nr_vertices: usize,
}

/// Grid cell identifier.
/// The identifier is (row, column).
///
/// ```shell
///   ^
/// Y |
///   | row-1     row-1
///   |
///   | row-0     row-0
///   | column-0  column-1
///   +--------------------->
///                       X
/// ```
///
#[derive(Copy, Clone, Hash, Debug, Ord, PartialOrd, PartialEq, Eq, Serialize, Deserialize)]
pub struct CellId {
    // A row is along the y-axis
    pub row: usize,
    // A column is along the x-axis
    pub column: usize,
}

impl Display for CellId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", &self.column, &self.row)
    }
}

/// 3D bounding box.
///
/// [min x, min y, min z, max x, max y, max z]
/// TODO: this must become a struct and have a .to_wkt() method
pub type Bbox = [f64; 6];

/// Serialize a 3D bounding box as 2D WKT Polygon.
pub fn bbox_to_wkt(bbox: &Bbox) -> String {
    format!(
        "POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
        minx = bbox[0],
        miny = bbox[1],
        maxx = bbox[3],
        maxy = bbox[4]
    )
}

/// 3D bounding box with quantized coordinates.
///
/// [min x, min y, min z, max x, max y, max z]
#[derive(Debug, Default, Clone, Ord, PartialOrd, Eq, PartialEq, Serialize, Deserialize)]
pub struct BboxQc(pub [i64; 6]); // This `pub [i64; 6]` makes the BboxQc constructor public

impl BboxQc {
    /// Compute the real-world coordinates from the quantized coordinates and the
    /// transformation properties.
    ///
    /// Optionally, the z-coordinate limits can be overriden by the `arg_minz` and
    /// `arg_maxz` arguments.
    pub fn to_bbox(
        &self,
        transform: &crate::parser::Transform,
        arg_minz: Option<i32>,
        arg_maxz: Option<i32>,
    ) -> Bbox {
        // Get the real-world coordinates for the extent
        let extent_rw_min = self.0[0..3]
            .iter()
            .enumerate()
            .map(|(i, qc)| (*qc as f64 * transform.scale[i]) + transform.translate[i]);
        let extent_rw_max = self.0[3..6]
            .iter()
            .enumerate()
            .map(|(i, qc)| (*qc as f64 * transform.scale[i]) + transform.translate[i]);
        let mut extent_rw: [f64; 6] = extent_rw_min
            .chain(extent_rw_max)
            .collect::<Vec<f64>>()
            .try_into()
            .expect("should be able to create an [f64; 6] from the extent_rw vector");
        if let Some(minz) = arg_minz {
            if extent_rw[2] < minz as f64 {
                debug!(
                "Setting min. z for the grid extent from provided value {}, instead of the computed value {}",
                minz, extent_rw[2]
            );
                extent_rw[2] = minz as f64
            }
        }
        if let Some(maxz) = arg_maxz {
            if extent_rw[5] > maxz as f64 {
                debug!(
                "Setting max. z for the grid extent from provided value {}, instead of the computed value {}",
                maxz, extent_rw[5]
            );
                extent_rw[5] = maxz as f64
            }
        }
        extent_rw
    }

    // Update with another bounding box, if the other is larger.
    pub fn update_with(&mut self, bbox_qc: &Self) {
        if bbox_qc.0[0] < self.0[0] {
            self.0[0] = bbox_qc.0[0]
        } else if bbox_qc.0[3] > self.0[3] {
            self.0[3] = bbox_qc.0[3]
        }
        if bbox_qc.0[1] < self.0[1] {
            self.0[1] = bbox_qc.0[1]
        } else if bbox_qc.0[4] > self.0[4] {
            self.0[4] = bbox_qc.0[4]
        }
        if bbox_qc.0[2] < self.0[2] {
            self.0[2] = bbox_qc.0[2]
        } else if bbox_qc.0[5] > self.0[5] {
            self.0[5] = bbox_qc.0[5]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use morton_encoding::morton_encode;

    #[test]
    fn test_intersect_bbox() {
        let extent = [195548.0, 538909.0, 0.0, 264268.0, 590410.0, 0.0];
        let grid = SquareGrid::new(&extent, 400, 7415);
        grid.export(None, None, None).unwrap();

        // Polygon ((248923.44474360189633444 601084.25658657902386039, 249381.04931766359368339 601093.95845033996738493, 249369.73047660905285738 601954.19037048425525427, 248923.44474360189633444 601084.25658657902386039))
        let bbox: Bbox = [248923.4, 601084.2, 0.0, 249381.0, 601954.1, 0.0];
        // 133-155, 133-156, 133-157, 134-155, 134-156, 134-157
        let expected = vec![
            CellId {
                row: 157,
                column: 133,
            },
            CellId {
                row: 157,
                column: 134,
            },
            CellId {
                row: 156,
                column: 133,
            },
            CellId {
                row: 156,
                column: 134,
            },
            CellId {
                row: 155,
                column: 133,
            },
            CellId {
                row: 155,
                column: 134,
            },
        ];
        let res = grid.intersect_bbox(&bbox);
        assert_eq!(expected.len(), res.len());
        assert!(res.iter().all(|res_cellid| expected.contains(res_cellid)))
    }

    #[test]
    fn test_create_grid() {
        let extent = [84372.91, 446316.814, -10.66, 171800.0, 472700.0, 52.882];
        let extent = [13603.33, 314127.708, -15.0, 268943.608, 612658.036, 400.0];
        println!("extent: {}", bbox_to_wkt(&extent));
        let grid = SquareGrid::new(&extent, 500, 7415);
        println!("grid: {}", bbox_to_wkt(&grid.bbox));
    }

    #[test]
    fn test_locate_point() {
        let grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 4.0], 1, 0);
        let cellid = grid.locate_point(&[2.5, 1.5]);
        println!("{}", cellid);
        assert_eq!(
            cellid,
            CellId {
                row: 1_usize,
                column: 2_usize,
            }
        );
    }

    #[test]
    fn test_morton_encode_rd() {
        let coords = vec![
            [84362.9, 446306.814],
            [84362.9, 446706.814],
            [84362.9, 447106.814],
            [84362.9, 447506.814],
            [84762.9, 446306.814],
            [84762.9, 446706.814],
            [84762.9, 447106.814],
            [84762.9, 447506.814],
            [85162.9, 446306.814],
            [85162.9, 446706.814],
            [85162.9, 447106.814],
            [85162.9, 447506.814],
            [85562.9, 446306.814],
            [85562.9, 446706.814],
            [85562.9, 447106.814],
            [85562.9, 447506.814],
        ];
        let row_col = vec![
            [0, 0],
            [0, 1],
            [0, 2],
            [0, 3],
            [1, 0],
            [1, 1],
            [1, 2],
            [1, 3],
            [2, 0],
            [2, 1],
            [2, 2],
            [2, 3],
            [3, 0],
            [3, 1],
            [3, 2],
            [3, 3],
        ];
        let grid = coords.iter().zip(row_col.iter());

        let mut morton: Vec<(u128, CellId)> = grid
            .map(|(coord, idx)| {
                let x: f64 = coord[0];
                let y: f64 = coord[1];
                let xu64 = x.floor() as u64;
                let yu64 = y.floor() as u64;
                let mc = morton_encode([xu64, yu64]);
                (
                    // interleave(&xu64, &yu64),
                    mc,
                    CellId {
                        row: idx[0],
                        column: idx[1],
                    },
                )
            })
            .collect();
        morton.sort_by_key(|k| k.0);

        let res_col_row_morton = [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
            [2, 0],
            [3, 0],
            [2, 1],
            [3, 1],
            [0, 2],
            [1, 2],
            [0, 3],
            [1, 3],
            [2, 2],
            [3, 2],
            [2, 3],
            [3, 3],
        ];
        let _res_coords_morton = [
            [84362.9, 446306.814],
            [84762.9, 446306.814],
            [84362.9, 446706.814],
            [84762.9, 446706.814],
            [85162.9, 446306.814],
            [85562.9, 446306.814],
            [85162.9, 446706.814],
            [85562.9, 446706.814],
            [84362.9, 447106.814],
            [84762.9, 447106.814],
            [84362.9, 447506.814],
            [84762.9, 447506.814],
            [85162.9, 447106.814],
            [85562.9, 447106.814],
            [85162.9, 447506.814],
            [85562.9, 447506.814],
        ];

        for (i, (mc, cellid)) in morton.iter().enumerate() {
            let res_cellid = CellId {
                column: res_col_row_morton[i][0],
                row: res_col_row_morton[i][1],
            };
            println!(
                "correct: {} computed: {} morton: {}",
                res_cellid, cellid, mc
            );
        }

        for (i, (_mc, cellid)) in morton.iter().enumerate() {
            let res_cellid = CellId {
                column: res_col_row_morton[i][0],
                row: res_col_row_morton[i][1],
            };
            assert_eq!(*cellid, res_cellid);
        }
    }

    #[test]
    fn test_morton_encode() {
        let mut cells_rowwise: Vec<(u64, u64)> = Vec::new();
        for x in 0..4_u64 {
            for y in 0..4u64 {
                cells_rowwise.push((x, y));
            }
        }
        let mut mortoncodes: Vec<u64> = cells_rowwise
            .iter()
            .map(|cell| interleave(&cell.0, &cell.1))
            .collect();
        mortoncodes.sort();
        let cells_morton: Vec<[u64; 2]> = mortoncodes.iter().map(deinterleave).collect();

        let expected: Vec<[u64; 2]> = vec![
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
            [2, 0],
            [3, 0],
            [2, 1],
            [3, 1],
            [0, 2],
            [1, 2],
            [0, 3],
            [1, 3],
            [2, 2],
            [3, 2],
            [2, 3],
            [3, 3],
        ];

        assert_eq!(cells_morton, expected)
    }

    #[test]
    fn test_quadtree_construction() {
        let mut feature_set: FeatureSet = Vec::new();
        let mut grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 1.0], 1, 0);
        for x in 0..4_u64 {
            for y in 0..4u64 {
                for f in 0..5 {
                    feature_set.push(crate::parser::Feature {
                        centroid_qc: [0, 0],
                        nr_vertices: 0,
                        path_jsonl: Default::default(),
                        bbox_qc: BboxQc([0, 0, 0, 0, 0, 0]),
                    });
                    let xc: f64 = format!("{}.{}", &x, &f).parse().unwrap();
                    grid.insert(&[xc, y as f64], f as usize);
                }
            }
        }
        let _ = QuadTree::from_grid(&grid, QuadTreeCapacity::Objects(20));
    }

    #[test]
    fn test_quadtree_leaves() {
        let mut feature_set: FeatureSet = Vec::new();
        let mut grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 1.0], 1, 0);
        for x in 0..4_u64 {
            for y in 0..4u64 {
                for f in 0..5 {
                    feature_set.push(crate::parser::Feature {
                        centroid_qc: [0, 0],
                        nr_vertices: 0,
                        path_jsonl: Default::default(),
                        bbox_qc: BboxQc([0, 0, 0, 0, 0, 0]),
                    });
                    let xc: f64 = format!("{}.{}", &x, &f).parse().unwrap();
                    grid.insert(&[xc, y as f64], f as usize);
                }
            }
        }
        let qtree = QuadTree::from_grid(&grid, QuadTreeCapacity::Objects(20));
        let leaves: Vec<&QuadTree> = QuadTree::collect_leaves(&qtree);
        for tile in leaves {
            println!("{}", tile.id);
        }
    }

    #[test]
    fn test_quadtree_node() {
        let mut feature_set: FeatureSet = Vec::new();
        let mut grid = SquareGrid::new(&[0.0, 0.0, 0.0, 16.0, 16.0, 1.0], 1, 0);
        for x in 0..16_u64 {
            for y in 0..16u64 {
                for f in 0..5 {
                    feature_set.push(crate::parser::Feature {
                        centroid_qc: [0, 0],
                        nr_vertices: 0,
                        path_jsonl: Default::default(),
                        bbox_qc: BboxQc([0, 0, 0, 0, 0, 0]),
                    });
                    let xc: f64 = format!("{}.{}", &x, &f).parse().unwrap();
                    grid.insert(&[xc, y as f64], f as usize);
                }
            }
        }
        let qtree = QuadTree::from_grid(&grid, QuadTreeCapacity::Objects(20));
        let _leaves: Vec<&QuadTree> = QuadTree::collect_leaves(&qtree);
        let n = qtree.node(&QuadTreeNodeId::new(0, 0, 2));
        if n.is_some() {
            println!("{:?}", &n);
        } else {
            println!("did not find node");
        }
    }
}
