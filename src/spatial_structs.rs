//! Spatial data structures for indexing the features.
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::prelude::*;

pub fn morton_encode(x: &f64, y: &f64) -> u128 {
    1
}

/// Represents a square grid with square cells.
/// The grid stores the feature-indices in its cells.
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
    origin: [f64; 2],
    pub nr_cells: usize,
    cellsize: u16,
    pub data: Vec<Vec<Cell>>,
}

impl Display for SquareGrid {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "SquareGrid (origin: {:?}, nr_cells: {}, cellsize: {}, data: not-displayed)",
               self.origin, self.nr_cells, self.cellsize)
    }
}

impl SquareGrid {
    /// Creates a grid with cells of `cellsize`, that covers the `extent`.
    /// The grid and the cells are square.
    /// The grid origin is the `extent` origin.
    /// The grid is returned as an origin coordinate and the number of cells.
    pub fn new(extent: &[f64; 6], cellsize: u16) -> Self {
        // Add some buffer to the extent, to make sure all points will be within the grid.
        // We are assuming quantized, metric coordinates with a scaling factor of 0.001, thus
        // 10 units translates to 10mm.
        let buffer = 10_f64;
        let dx = (extent[3] - extent[0]).abs() + buffer * 2.0;
        let dy = (extent[4] - extent[1]).abs() + buffer * 2.0;
        let gridsize = if dx > dy { dx } else { dy };
        let nr_cells = (gridsize / cellsize as f64).ceil() as usize;
        // A row-vector (x-axis) to store the column-vectors (y-axis).
        let mut row: Vec<Vec<Vec<usize>>> = Vec::with_capacity(nr_cells);
        // For each column create a column vector that stores the cells and for each row in the
        // column create a cell to store the feature IDs.
        // We create the vectors with length nr_cells+1, because we use `ceil` to determine the
        // point location in the grid.
        row.resize_with(nr_cells+1, || {
            let mut column: Vec<Vec<usize>> = Vec::with_capacity(nr_cells);
            column.resize(nr_cells+1, Vec::new());
            column
        });
        Self {
            origin: [extent[0] - buffer, extent[1] - buffer],
            nr_cells,
            cellsize,
            data: row,
        }
    }

    /// Returns the cell index (x, y) where the point is located.
    fn locate_point(&self, point: &[f64; 2]) -> CellId {
        let dx = point[0] - self.origin[0];
        let dy = point[1] - self.origin[1];
        let row_i = (dx / self.cellsize as f64).ceil() as usize;
        let col_i = (dy / self.cellsize as f64).ceil() as usize;
        [row_i, col_i]
    }

    pub fn insert(&mut self, point: &[f64; 2], feature_id: usize) -> CellId {
        let cell_id = self.locate_point(point);
        self.data[cell_id[0]][cell_id[1]].push(feature_id);
        cell_id
    }

    /// Exports the grid and the feature centroids into TSV files into the working directory.
    /// Two files are created, `grid.tsv` and `features.tsv`.
    pub fn export(&self, feature_set: &crate::FeatureSet, cm: &crate::parser::CityJSONMetadata) -> std::io::Result<()> {
        let mut file_grid = File::create("grid.tsv")?;
        let mut file_features = File::create("features.tsv")?;

        for (col_i, col) in self.data.iter().enumerate() {
            for (row_i, cell) in col.iter().enumerate() {
                let wkt = self.cell_to_wkt(&row_i, &col_i);
                file_grid.write_all(format!("{}-{}\t{}\n", &row_i, &col_i, wkt).as_bytes()).expect("cannot write grid line");
                let mut cellbuffer = String::new();
                for fid in cell {
                    let f = &feature_set[*fid];
                    let centroid = f.centroid(cm);
                    cellbuffer += format!("{}\t{}-{}\tPOINT({} {})\n", fid, &row_i, &col_i, centroid[0], centroid[1]).as_str();
                }
                file_features.write_all(cellbuffer.as_bytes()).expect("cannot write cell contents");
            }
        }
        Ok(())
    }

    fn cell_to_wkt(&self, row_i: &usize, col_i: &usize) -> String {
        let minx = self.origin[0] + (row_i * self.cellsize as usize) as f64;
        let miny = self.origin[1] + (col_i * self.cellsize as usize) as f64;
        format!("POLYGON(({minx} {miny}, {maxx} {miny}, {maxx} {maxy}, {minx} {maxy}, {minx} {miny}))",
                minx=minx, miny=miny,
                maxx=minx + self.cellsize as f64, maxy=miny + self.cellsize as f64)
    }
}

type Cell = Vec<usize>;
pub type CellId = [usize; 2];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_grid() {
        let extent = [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
        let grid = SquareGrid::new(&extent, 500);
        println!("grid: {:?}", grid);
    }

    #[test]
    fn test_locate_point() {
        let grid = SquareGrid::new(&[0.0, 0.0, 0.0, 4.0, 4.0, 4.0], 1);
        let grid_idx = grid.locate_point(&[2.5, 1.5]);
        assert_eq!(grid_idx, [3_usize, 2_usize]);
    }
}