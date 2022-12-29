pub fn morton_encode(x: &f64, y: &f64) -> u128 {
    1
}

/// Creates a grid with cells of `cellsize`, that covers the `extent`.
/// The grid and the cells are square.
/// The grid origin is the `extent` origin.
/// The grid is returned as an origin coordinate and the number of cells.
fn create_grid(extent: &[f64; 6], cellsize: u16) -> ([f64; 2], u64) {
    // Add some buffer to the extent, to make sure all points will be within the grid.
    // We are assuming quantized, metric coordinates with a scaling factor of 0.001, thus
    // 10 units translates to 10mm.
    let buffer = 10_f64;
    let dx = (extent[3] - extent[0]).abs() + buffer*2.0;
    let dy = (extent[4] - extent[1]).abs() + buffer*2.0;
    let gridsize = if dx > dy { dx } else { dy };
    let nr_cells = (gridsize / cellsize as f64).ceil() as u64;
    ([extent[0] - buffer, extent[1] - buffer], nr_cells)
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
}