pub fn morton_encode(x: &f64, y: &f64) -> u128 {
    1
}

/// Represents a square grid with square cells.
#[derive(Debug)]
struct SquareGrid {
    origin: [f64; 2],
    nr_cells: u64,
    cellsize: u16
}

impl SquareGrid {
    // Returns the cell index (x, y) where the point is located.
    fn locate_point(&self, point: &[f64; 2]) -> [u64; 2] {
        let dx = point[0] - self.origin[0];
        let dy = point[1] - self.origin[1];
        [(dx / self.cellsize as f64).ceil() as u64, (dy / self.cellsize as f64).ceil() as u64]
    }
}

/// Creates a grid with cells of `cellsize`, that covers the `extent`.
/// The grid and the cells are square.
/// The grid origin is the `extent` origin.
/// The grid is returned as an origin coordinate and the number of cells.
fn create_grid(extent: &[f64; 6], cellsize: u16) -> SquareGrid {
    // Add some buffer to the extent, to make sure all points will be within the grid.
    // We are assuming quantized, metric coordinates with a scaling factor of 0.001, thus
    // 10 units translates to 10mm.
    let buffer = 10_f64;
    let dx = (extent[3] - extent[0]).abs() + buffer*2.0;
    let dy = (extent[4] - extent[1]).abs() + buffer*2.0;
    let gridsize = if dx > dy { dx } else { dy };
    let nr_cells = (gridsize / cellsize as f64).ceil() as u64;
    SquareGrid {
        origin: [extent[0] - buffer, extent[1] - buffer],
        nr_cells,
        cellsize
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_grid() {
        let extent = [84995.279, 446316.813, -5.333, 85644.748, 446996.132, 52.881];
        let grid = create_grid(&extent, 500);
        println!("grid: {:?}", grid);
    }

    #[test]
    fn test_locate_point() {
        let grid = SquareGrid {origin: [0.0, 0.0], nr_cells: 3, cellsize: 1};
        let grid_idx = grid.locate_point(&[2.5, 1.5]);
        assert_eq!(grid_idx, [3_u64, 2_u64]);
    }
}