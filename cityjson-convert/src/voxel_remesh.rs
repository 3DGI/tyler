//! Rust voxel reconstruction and checked native simplification.
//!
//! Reconstruction ports meshoptimizer v1.3, commit
//! `9e1f07b159d3cb777f1c67ed31fc11fd117986f4`, under the MIT licence retained in
//! `voxel_remesh/LICENSE.md`. The checked simplifier uses upstream meshopt's
//! native position-update API independently of the Rust reconstruction.

mod grid;
mod reconstruct;
mod tables;

use anyhow::{anyhow, Result};
pub use grid::{VoxelExtent, VoxelGrid};
use meshopt::{ffi, SimplifyOptions};

bitflags::bitflags! {
    /// Options for voxel reconstruction.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub struct RemeshOptions: u32 {
        /// Wrap surfaces in a two-sided shell instead of filling closed interiors.
        const Shell = 1;
        /// Solve output positions to approximate source surfaces.
        const Solve = 2;
    }
}

fn validate_mesh(indices: &[u32], positions: &[[f32; 3]]) -> Result<()> {
    if !indices.len().is_multiple_of(3) {
        return Err(anyhow!("indices must describe complete triangles"));
    }
    if positions.len() > u32::MAX as usize
        || indices
            .iter()
            .any(|&index| index as usize >= positions.len())
    {
        return Err(anyhow!("vertex index out of range"));
    }
    if positions.iter().flatten().any(|value| !value.is_finite()) {
        return Err(anyhow!("positions must be finite"));
    }
    Ok(())
}

fn allocate<T: Clone>(count: usize, initial: T) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|error| anyhow!(error.to_string()))?;
    values.resize(count, initial);
    Ok(values)
}

/// Indexed reconstruction on an explicitly defined processing domain.
/// Cell provenance belongs to this unsimplified mesh; indices may share positions.
#[derive(Debug)]
pub struct RemeshResult {
    pub positions: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
    pub grid_dimensions: [u32; 3],
    pub grid: VoxelGrid,
    pub extent: VoxelExtent,
    vertex_cells: Vec<usize>,
}

impl RemeshResult {
    /// Global lattice identity of the cell that produced this vertex.
    /// Identities remain comparable across processing extents on the same grid.
    #[must_use]
    pub fn vertex_cell(&self, vertex: usize) -> Option<[i64; 3]> {
        let flat = *self.vertex_cells.get(vertex)?;
        let stride = (self.grid_dimensions[0] as usize).div_ceil(64) * 64;
        let row = flat / stride;
        Some([
            self.extent.min[0] + i64::try_from(flat % stride).ok()?,
            self.extent.min[1] + i64::try_from(row % self.grid_dimensions[1] as usize).ok()?,
            self.extent.min[2] + i64::try_from(row / self.grid_dimensions[1] as usize).ok()?,
        ])
    }

    /// Consumes the reconstruction and releases cell provenance.
    #[must_use]
    pub fn into_mesh(self) -> (Vec<[f32; 3]>, Vec<u32>) {
        (self.positions, self.indices)
    }
}

/// Reconstructs an indexed mesh using an explicit lattice and processing extent.
/// Builds occupancy once, samples twice and only visits nonuniform output cubes.
/// The processing extent is not a clipping boundary. Solid reconstruction requires
/// full input coverage; regional shell reconstruction discards out-of-domain samples.
///
/// # Errors
/// Rejects malformed meshes, invalid grids/options, unusable coordinate extents,
/// unsupported cropped solid domains, index overflows and allocation failures.
pub fn remesh(
    indices: &[u32],
    positions: &[[f32; 3]],
    grid: VoxelGrid,
    extent: VoxelExtent,
    options: RemeshOptions,
) -> Result<RemeshResult> {
    validate_mesh(indices, positions)?;
    grid.validate()?;
    if options.bits() & !RemeshOptions::all().bits() != 0 {
        return Err(anyhow!("invalid remeshing options"));
    }
    reconstruct::reconstruct(indices, positions, grid, extent, options)
}

/// Simplifies and updates positions in place, returning the reported error.
/// `vertex_locks`, when supplied, has one byte per position with native
/// `meshopt_SimplifyVertex_*` flags (1 locks a vertex's position and connectivity).
///
/// # Errors
/// Rejects malformed meshes, invalid triangle targets/tolerances and invalid locks.
pub fn simplify_positions_with_update(
    indices: &mut Vec<u32>,
    positions: &mut [[f32; 3]],
    vertex_locks: Option<&[u8]>,
    target_index_count: usize,
    target_error: f32,
    options: SimplifyOptions,
) -> Result<f32> {
    validate_mesh(indices, positions)?;
    if target_index_count > indices.len()
        || !target_index_count.is_multiple_of(3)
        || !target_error.is_finite()
        || target_error < 0.0
        || options.bits() & !SimplifyOptions::all().bits() != 0
    {
        return Err(anyhow!("invalid simplification target, error or options"));
    }
    if vertex_locks.is_some_and(|locks| {
        locks.len() != positions.len() || locks.iter().any(|lock| lock & !7 != 0)
    }) {
        return Err(anyhow!("invalid simplification vertex locks"));
    }
    if indices.is_empty() {
        return Ok(0.0);
    }
    let mut error = 0.0;
    // SAFETY: valid initialized index/float3 buffers, no attribute streams, and
    // legal targets. The output index count cannot exceed the input count.
    let count = unsafe {
        ffi::meshopt_simplifyWithUpdate(
            indices.as_mut_ptr(),
            indices.len(),
            positions.as_mut_ptr().cast(),
            positions.len(),
            12,
            std::ptr::null_mut(),
            0,
            std::ptr::null(),
            0,
            vertex_locks.map_or(std::ptr::null(), <[u8]>::as_ptr),
            target_index_count,
            target_error,
            options.bits(),
            &raw mut error,
        )
    };
    indices.truncate(count);
    validate_mesh(indices, positions)?;
    if !error.is_finite() || error < 0.0 {
        return Err(anyhow!("native simplification returned an invalid error"));
    }
    Ok(error)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    /// Purpose: validate indexed reconstruction, wide packed rows, global cell
    /// identities across regional extents, cropped-solid rejection and vertex locks.
    #[test]
    fn remesh_checks_inputs_and_sizes_output() {
        check_inputs_and_indexing();
        check_wide_grids_and_regional_cells();
        check_simplifier_locks();
    }

    fn check_inputs_and_indexing() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let frame = VoxelGrid {
            origin: [0.0; 3],
            voxel_size: 0.2,
        };
        let extent = frame.cover(&positions).unwrap();
        for indices in [&[0, 1][..], &[0, 1, 3][..]] {
            assert!(remesh(indices, &positions, frame, extent, RemeshOptions::Solve).is_err());
        }
        for size in [
            0.0,
            -1.0,
            f32::NAN,
            f32::INFINITY,
            f32::from_bits(1),
            1.0e-8,
        ] {
            assert!(VoxelGrid::centered(&positions, size)
                .and_then(|grid| grid.cover(&positions))
                .is_err());
        }
        assert!(remesh(
            &[0, 1, 2],
            &positions,
            frame,
            extent,
            RemeshOptions::from_bits_retain(1 << 31)
        )
        .is_err());
        let mut invalid = positions;
        invalid[0][0] = f32::NAN;
        assert!(remesh(&[0, 1, 2], &invalid, frame, extent, RemeshOptions::Solve).is_err());
        invalid[0][0] = f32::MAX;
        assert!(remesh(&[0, 1, 2], &invalid, frame, extent, RemeshOptions::Solve).is_err());
        let bad_extent = VoxelExtent {
            min: [i64::MIN; 3],
            max: [i64::MAX; 3],
        };
        assert!(remesh(
            &[0, 1, 2],
            &positions,
            frame,
            bad_extent,
            RemeshOptions::Solve
        )
        .is_err());
        let mut invalid_frame = frame;
        invalid_frame.origin[0] = f64::NAN;
        assert!(remesh(
            &[0, 1, 2],
            &positions,
            invalid_frame,
            extent,
            RemeshOptions::Solve
        )
        .is_err());
        assert!(remesh(&[], &[], frame, extent, RemeshOptions::Solve)
            .unwrap()
            .indices
            .is_empty());
        let output = remesh(&[0, 1, 2], &positions, frame, extent, RemeshOptions::Solve).unwrap();
        assert!(!output.indices.is_empty());
        assert_eq!(output.indices.len() % 3, 0);
        assert_eq!(output.grid_dimensions, [8, 8, 4]);
        assert!(output
            .indices
            .iter()
            .all(|&index| (index as usize) < output.positions.len()));
        assert!(output
            .positions
            .iter()
            .flatten()
            .all(|value| value.is_finite()));
        assert!(output.vertex_cell(output.positions.len()).is_none());
        let negative = positions.map(|point| point.map(|value| value - 1.25));
        let negative_extent = frame.cover(&negative).unwrap();
        let negative_mesh = remesh(
            &[0, 1, 2],
            &negative,
            frame,
            negative_extent,
            RemeshOptions::Solve,
        )
        .unwrap();
        assert!(!negative_mesh.indices.is_empty());
        for vertex in 0..negative_mesh.positions.len() {
            assert!(negative_mesh
                .vertex_cell(vertex)
                .unwrap()
                .iter()
                .all(|&cell| cell < 0));
        }
    }

    fn check_wide_grids_and_regional_cells() {
        for axis in 0..3 {
            let mut elongated = [[0.0; 3]; 4];
            elongated[1][axis] = 1100.0;
            elongated[2][(axis + 1) % 3] = 4.0;
            elongated[3][axis] = 1100.0;
            elongated[3][(axis + 1) % 3] = 4.0;
            let frame = VoxelGrid {
                origin: [0.0; 3],
                voxel_size: 1.0,
            };
            let extent = frame.cover(&elongated).unwrap();
            let indices = [0, 1, 2, 2, 1, 3];
            for options in [
                RemeshOptions::Solve,
                RemeshOptions::Solve | RemeshOptions::Shell,
            ] {
                let output = remesh(&indices, &elongated, frame, extent, options).unwrap();
                assert_eq!(output.grid_dimensions[axis], 1103);
                assert_eq!(output.grid_dimensions[(axis + 1) % 3], 7);
                assert_eq!(output.grid_dimensions[(axis + 2) % 3], 4);
                assert!(output
                    .positions
                    .iter()
                    .flatten()
                    .all(|value| value.is_finite()));
                assert!(output.positions.iter().any(|point| point[axis] > 1099.0));
                for point in output.positions {
                    for (coordinate, &value) in point.iter().enumerate() {
                        let maximum = if coordinate == axis {
                            1100.0
                        } else if coordinate == (axis + 1) % 3 {
                            4.0
                        } else {
                            0.0
                        };
                        // A fixed shared lattice may place the source extremum
                        // on a cell edge; reconstruction can move within that cell.
                        assert!((-0.001..=maximum + 1.001).contains(&value));
                    }
                }
            }
            // The same original triangles and lattice must produce exactly the
            // same local reconstruction in the overlap, regardless of extent.
            let full = remesh(
                &indices,
                &elongated,
                frame,
                extent,
                RemeshOptions::Shell | RemeshOptions::Solve,
            )
            .unwrap();
            let reference: BTreeMap<_, _> = full
                .positions
                .iter()
                .enumerate()
                .map(|(index, &point)| (full.vertex_cell(index).unwrap(), point))
                .collect();
            let mut cropped = extent;
            cropped.min[axis] = 300;
            cropped.max[axis] = 700;
            assert!(remesh(&indices, &elongated, frame, cropped, RemeshOptions::Solve).is_err());
            let regional = remesh(
                &indices,
                &elongated,
                frame,
                cropped,
                RemeshOptions::Shell | RemeshOptions::Solve,
            )
            .unwrap();
            assert!(!regional.indices.is_empty());
            for (vertex, &position) in regional.positions.iter().enumerate() {
                assert_eq!(
                    reference[&regional.vertex_cell(vertex).unwrap()].map(f32::to_bits),
                    position.map(f32::to_bits),
                );
            }
        }
    }

    fn check_simplifier_locks() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut indices = vec![0, 1, 2];
        let mut points = positions;
        for (target, error) in [(4, 0.1), (6, 0.1), (0, f32::NAN), (0, -1.0)] {
            assert!(simplify_positions_with_update(
                &mut indices,
                &mut points,
                None,
                target,
                error,
                SimplifyOptions::ErrorAbsolute
            )
            .is_err());
        }
        for locks in [&[1, 1][..], &[8, 1, 1][..]] {
            assert!(simplify_positions_with_update(
                &mut indices,
                &mut points,
                Some(locks),
                0,
                1.0,
                SimplifyOptions::ErrorAbsolute
            )
            .is_err());
        }
        assert!(simplify_positions_with_update(
            &mut indices,
            &mut points,
            Some(&[1, 1, 1]),
            0,
            1.0,
            SimplifyOptions::ErrorAbsolute
        )
        .unwrap()
        .is_finite());
        assert_eq!(indices, [0, 1, 2]);
        assert_eq!(points, positions);
    }
}
