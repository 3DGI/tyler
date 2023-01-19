//! Spatial data structures for indexing the features.
use crate::parser::Feature;
use crate::FeatureSet;
use log::{debug, error};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::prelude::*;

/// Quadtree
///
/// We don't expect that the quadtree has more than 65535 levels (u16).
#[derive(Clone, Debug)]
pub struct QuadTree {
    x: usize,
    y: usize,
    z: u16,
    side_length: u64,
    pub children: Vec<QuadTree>,
    pub cells: Vec<CellId>,
    pub nr_items: usize,
}

impl QuadTree {
    pub fn from_grid(grid: &SquareGrid, feature_set: &FeatureSet, limit: QuadTreeLimit) -> Self {
        let mut merge_limit: usize = 0;
        let nr_cells = grid.length.pow(2) as f64;
        let max_level = (nr_cells.ln() / 4.0_f64.ln()).ceil() as u16;
        debug!("Calculated maximum level for quadtree: {}", &max_level);
        let mut mortoncodes: Vec<u64> = grid
            .into_iter()
            .map(|(cellid, _)| interleave(&(cellid.column as u64), &(cellid.row as u64)))
            .collect();
        mortoncodes.sort();
        let tiles_morton: Vec<QuadTree> = mortoncodes
            .iter()
            .map(|mc| deinterleave(mc))
            .map(|[x, y]| {
                let cellid = CellId {
                    row: y as usize,
                    column: x as usize,
                };
                let items: usize;
                match limit {
                    QuadTreeLimit::Objects(l) => {
                        // Use the number of features as a limit
                        items = grid.cell(&cellid).feature_ids.len();
                        merge_limit = l;
                    }
                    QuadTreeLimit::Vertices(l) => {
                        // Use the number of vertices as a limit
                        items = grid.cell(&cellid).nr_vertices;
                        merge_limit = l;
                    }
                }
                QuadTree {
                    x: x as usize,
                    y: y as usize,
                    z: max_level,
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
            if sum_items <= limit {
                QuadTree {
                    x: tiles[0].x,
                    y: tiles[0].y,
                    z: level,
                    side_length: tiles[0].side_length * 2,
                    children: vec![],
                    cells,
                    nr_items: sum_items,
                }
            } else {
                if tiles.len() % 4 != 0 {
                    error!(
                        "number of children is not divisable by 4: {}, in level {}",
                        tiles.len(),
                        &level
                    );
                }
                QuadTree {
                    x: tiles[0].x,
                    y: tiles[0].y,
                    z: level,
                    side_length: tiles[0].side_length * 2,
                    children: tiles.clone(),
                    cells: vec![],
                    nr_items: sum_items,
                }
            }
        }
    }

    fn visit_leaves_helper(quadtree: &QuadTree) {
        if !quadtree.children.is_empty() {
            Self::visit_leaves_helper(&quadtree.children[0]);
            Self::visit_leaves_helper(&quadtree.children[1]);
            Self::visit_leaves_helper(&quadtree.children[2]);
            Self::visit_leaves_helper(&quadtree.children[3]);
        } else {
            println!("leaf {}, items: {}", quadtree.id(), quadtree.nr_items);
        }
    }

    pub fn visit_leaves(&self) {
        Self::visit_leaves_helper(self);
    }

    fn collect_leaves_recurse<'collect>(&'collect self, leaves: &mut Vec<&'collect QuadTree>) {
        if !self.children.is_empty() {
            for child in self.children.iter() {
                child.collect_leaves_recurse(leaves);
            }
        } else {
            leaves.push(self);
        }
    }

    pub fn collect_leaves(&self) -> Vec<&Self> {
        let mut leaves: Vec<&QuadTree> = Vec::new();
        self.collect_leaves_recurse(&mut leaves);
        leaves
    }

    // pub fn collect_features(&self, grid: &SquareGrid, feature_set: &FeatureSet) -> Vec<(String, Vec<&Feature>)> {
    //     let mut features: Vec<(String, Vec<&Feature>)> = Vec::new();
    //     self.collect_features_recurse(&mut features, grid,feature_set);
    //     features
    // }
    //
    // fn collect_features_recurse<'collect>(&'collect self, features: &mut Vec<(String, Vec<&'collect Feature>)>, grid: &SquareGrid, feature_set: &'collect Vec<crate::parser::Feature>) {
    //     if !self.children.is_empty() {
    //         for child in self.children.iter() {
    //             child.collect_features_recurse(features, grid, feature_set);
    //         }
    //     } else {
    //         let mut features_in_leaf: Vec<&Feature> = Vec::with_capacity(self.nr_items);
    //         for cellid in &self.cells {
    //             for fid in grid.cell(cellid) {
    //                 features_in_leaf.push(&feature_set[*fid]);
    //             }
    //         }
    //         features.push((self.id(), features_in_leaf));
    //     }
    // }

    pub fn id(&self) -> String {
        format!("{}/{}/{}", self.z, self.x, self.y)
    }

    pub fn bbox(&self, grid: &SquareGrid) -> crate::Bbox {
        let minx = grid.origin[0] + (self.x * grid.cellsize as usize) as f64;
        let miny = grid.origin[1] + (self.y * grid.cellsize as usize) as f64;
        [
            minx,
            miny,
            grid.bbox[2],
            minx + self.side_length as f64,
            miny + self.side_length as f64,
            grid.bbox[5],
        ]
    }
}

#[derive(Debug)]
pub enum QuadTreeLimit {
    Objects(usize),
    Vertices(usize),
}

/// 64-bit mask
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
pub fn interleave(x: &u64, y: &u64) -> u64 {
    part1by1_64(x) | (part1by1_64(y) << 1)
}

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
#[derive(Debug)]
pub struct SquareGrid {
    origin: [f64; 3],
    pub bbox: crate::Bbox,
    pub length: usize,
    cellsize: u16,
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
    /// The grid origin is the `extent` origin.
    /// The grid is returned as an origin coordinate and the number of cells.
    pub fn new(extent: &crate::Bbox, cellsize: u16, epsg: u16, buffer: Option<f64>) -> Self {
        // Add some buffer to the extent, to make sure all points will be within the grid.
        let buffer: f64 = buffer.unwrap_or(0.0);
        // Add the buffer to the computed extent
        let extent_with_buffer = [
            extent[0] - buffer,
            extent[1] - buffer,
            extent[2] - buffer,
            extent[3] + buffer,
            extent[4] + buffer,
            extent[5] + buffer,
        ];
        let dx = extent_with_buffer[3] - extent_with_buffer[0];
        let dy = extent_with_buffer[4] - extent_with_buffer[1];
        // The grid dimension is the longest edge of the rectangle, so we get a square
        let mut d = if dx > dy { dx } else { dy };
        // We need a grid that is has 2^n cells in one dimension, so that we can build
        // a 4^n cells quadtree.
        let d_cells = 2_usize.pow((d / cellsize as f64).log2().ceil() as u32);
        debug!("Computed grid cells dimension: {}", &d_cells);
        let origin = [
            extent_with_buffer[0],
            extent_with_buffer[1],
            extent_with_buffer[2],
        ];
        // Compute new dimension from the calculated length
        d = d_cells as f64 * cellsize as f64;
        let bbox = [
            origin[0],
            origin[1],
            origin[2],
            origin[0] + d,
            origin[1] + d,
            extent_with_buffer[5],
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

    pub fn insert(&mut self, point: &[f64; 2], feature_id: usize) -> CellId {
        let cell_id = self.locate_point(point);
        let cell = self.cell_mut(&cell_id);
        cell.feature_ids.push(feature_id);
        cell_id
    }

    /// Exports the grid and the feature centroids into TSV files into the working directory.
    /// Two files are created, `grid.tsv` and `features.tsv`.
    pub fn export(
        &self,
        feature_set: &crate::FeatureSet,
        cm: &crate::parser::CityJSONMetadata,
    ) -> std::io::Result<()> {
        let mut file_grid = File::create("grid.tsv")?;
        let mut file_features = File::create("features.tsv")?;

        let root_wkt = format!(
            "POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
            minx = self.bbox[0],
            miny = self.bbox[1],
            maxx = self.bbox[3],
            maxy = self.bbox[4]
        );
        file_grid
            .write_all(format!("x-x\t0\t{}\n", root_wkt).as_bytes())
            .expect("cannot write grid line");
        for (cellid, cell) in self {
            let wkt = self.cell_to_wkt(&cellid);
            file_grid
                .write_all(format!("{}\t{}\t{}\n", &cellid, cell.nr_vertices, wkt).as_bytes())
                .expect("cannot write grid line");
            let mut cellbuffer = String::new();
            for fid in cell.feature_ids.iter() {
                let f = &feature_set[*fid];
                let centroid = f.centroid(cm);
                cellbuffer += format!(
                    "{}\t{}\tPOINT({} {})\n",
                    fid, &cellid, centroid[0], centroid[1]
                )
                .as_str();
            }
            file_features
                .write_all(cellbuffer.as_bytes())
                .expect("cannot write cell contents");
        }
        Ok(())
    }

    fn cell_to_wkt(&self, cellid: &CellId) -> String {
        let minx = self.origin[0] + (cellid.column * self.cellsize as usize) as f64;
        let miny = self.origin[1] + (cellid.row * self.cellsize as usize) as f64;
        format!(
            "POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
            minx = minx,
            miny = miny,
            maxx = minx + self.cellsize as f64,
            maxy = miny + self.cellsize as f64
        )
    }

    pub fn cell_bbox(&self, cellid: &CellId) -> crate::Bbox {
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

#[derive(Debug, Clone, Ord, PartialOrd, Eq, PartialEq)]
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
#[derive(Copy, Clone, Debug, Ord, PartialOrd, PartialEq, Eq)]
pub struct CellId {
    // A row is along the y-axis
    pub row: usize,
    // A column is along the x-axis
    pub column: usize,
}

impl Display for CellId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", &self.row, &self.column)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_grid() {
        let extent = [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
        let grid = SquareGrid::new(&extent, 500, 7415, Some(0.0));
        println!("grid: {:?}", grid);
    }

    #[test]
    fn test_locate_point() {
        let grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 4.0], 1, 0, Some(0.0));
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
        let cells_morton: Vec<[u64; 2]> = mortoncodes.iter().map(|mc| deinterleave(mc)).collect();

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
        let mut feature_set: crate::FeatureSet = Vec::new();
        let mut grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 1.0], 1, 0, None);
        for x in 0..4_u64 {
            for y in 0..4u64 {
                for f in 0..5 {
                    feature_set.push(crate::parser::Feature {
                        centroid_quantized: [0, 0],
                        nr_vertices: 0,
                        path_jsonl: Default::default(),
                        bbox_quantized: [0, 0, 0, 0, 0, 0],
                    });
                    let xc: f64 = format!("{}.{}", &x, &f).parse().unwrap();
                    grid.insert(&[xc, y as f64], f as usize);
                }
            }
        }
        let qtree = QuadTree::from_grid(&grid, &feature_set, QuadTreeLimit::Objects(20));
        qtree.visit_leaves();
    }

    #[test]
    fn test_quadtree_leaves() {
        let mut feature_set: crate::FeatureSet = Vec::new();
        let mut grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 1.0], 1, 0, None);
        for x in 0..4_u64 {
            for y in 0..4u64 {
                for f in 0..5 {
                    feature_set.push(crate::parser::Feature {
                        centroid_quantized: [0, 0],
                        nr_vertices: 0,
                        path_jsonl: Default::default(),
                        bbox_quantized: [0, 0, 0, 0, 0, 0],
                    });
                    let xc: f64 = format!("{}.{}", &x, &f).parse().unwrap();
                    grid.insert(&[xc, y as f64], f as usize);
                }
            }
        }
        let qtree = QuadTree::from_grid(&grid, &feature_set, QuadTreeLimit::Objects(20));
        let leaves: Vec<&QuadTree> = QuadTree::collect_leaves(&qtree);
        for tile in leaves {
            println!("{}", tile.id());
        }
    }
}
