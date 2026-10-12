use anyhow::{bail, ensure, Context, Result};
use bvh::aabb::{Aabb, Bounded};
use bvh::ball::Ball;
use bvh::bounding_hierarchy::BHShape;
use bvh::bvh::Bvh;
use bvh::point_query::PointDistance;
use meshopt::{generate_vertex_remap, remap_index_buffer, remap_vertex_buffer};
use nalgebra::Point3;

use crate::voxel_remesh::{self, RemeshOptions, VoxelGrid};
use crate::{RemeshReport, VoxelRemeshOptions};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Triangle {
    pub positions: [[f32; 3]; 3],
    pub feature_id: u32,
}

struct SourceTriangle {
    points: [Point3<f64>; 3],
    bounds: Aabb<f64, 3>,
    feature_id: u32,
    /// Source order is sorted by object ID before BVH construction.
    ordinal: usize,
    node: usize,
}

impl Bounded<f64, 3> for SourceTriangle {
    fn aabb(&self) -> Aabb<f64, 3> {
        self.bounds
    }
}

impl BHShape<f64, 3> for SourceTriangle {
    fn set_bh_node_index(&mut self, index: usize) {
        self.node = index;
    }

    fn bh_node_index(&self) -> usize {
        self.node
    }
}

impl PointDistance<f64, 3> for SourceTriangle {
    fn distance_squared(&self, point: Point3<f64>) -> f64 {
        // Project onto the triangle plane, then use barycentric coordinates to
        // distinguish its interior from the three closest-edge regions.
        let [a, b, c] = self.points;
        let ab = b - a;
        let ac = c - a;
        let normal = ab.cross(&ac);
        let normal_squared = normal.norm_squared();
        if normal_squared > 0.0 {
            let projection = point - normal * ((point - a).dot(&normal) / normal_squared);
            let ap = projection - a;
            let barycentric_b = ap.cross(&ac).dot(&normal) / normal_squared;
            let barycentric_c = ab.cross(&ap).dot(&normal) / normal_squared;
            if barycentric_b >= 0.0 && barycentric_c >= 0.0 && barycentric_b + barycentric_c <= 1.0
            {
                return (point - projection).norm_squared();
            }
        }
        [(a, b), (b, c), (c, a)]
            .into_iter()
            .map(|(start, end)| {
                let edge = end - start;
                let length = edge.norm_squared();
                let t = if length > 0.0 {
                    ((point - start).dot(&edge) / length).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (point - (start + edge * t)).norm_squared()
            })
            .fold(f64::INFINITY, f64::min)
    }
}

struct SourceSurface {
    triangles: Vec<SourceTriangle>,
    bvh: Bvh<f64, 3>,
}

impl SourceSurface {
    fn new(input: &[Triangle]) -> Self {
        let mut triangles = input
            .iter()
            .enumerate()
            .map(|(ordinal, triangle)| {
                let points = triangle.positions.map(|p| Point3::from(p.map(f64::from)));
                let mut bounds = Aabb::empty();
                for point in &points {
                    bounds.grow_mut(point);
                }
                SourceTriangle {
                    points,
                    bounds,
                    feature_id: triangle.feature_id,
                    ordinal,
                    node: 0,
                }
            })
            .collect::<Vec<_>>();
        let bvh = Bvh::build(&mut triangles);
        Self { triangles, bvh }
    }

    fn owner(&self, positions: [[f32; 3]; 3]) -> Result<u32> {
        let points = positions.map(|p| Point3::from(p.map(f64::from)));
        let centroid = Point3::from((points[0].coords + points[1].coords + points[2].coords) / 3.0);
        let (nearest, distance) = self
            .bvh
            .nearest_to(centroid, &self.triangles)
            .context("no source surface available for feature attribution")?;
        // Expand by roundoff so exact ties on coincident surfaces are included.
        let radius = distance + f64::EPSILON * distance.max(1.0) * 8.0;
        let owner = self
            .bvh
            .traverse(&Ball::new(centroid, radius), &self.triangles)
            .into_iter()
            .min_by(|a, b| {
                a.distance_squared(centroid)
                    .total_cmp(&b.distance_squared(centroid))
                    .then_with(|| a.ordinal.cmp(&b.ordinal))
            })
            .unwrap_or(nearest);
        Ok(owner.feature_id)
    }
}

fn voxel_size(error: f32) -> Result<f32> {
    ensure!(
        error.is_finite() && error > 0.0,
        "remeshing geometric error must be finite and positive"
    );
    let target = f64::from(error) / (2.0 * 3.0_f64.sqrt());
    let size = num_traits::ToPrimitive::to_f32(&target)
        .context("remeshing voxel size cannot be represented")?;
    // Round toward finer quality so the voxel diagonal stays within half the
    // requested error budget even when the size is rounded to float.
    let size = if f64::from(size) > target {
        size.next_down()
    } else {
        size
    };
    ensure!(
        size > 0.0 && (2.0 / size).is_finite(),
        "remeshing error is too small for float voxel coordinates"
    );
    Ok(size)
}

pub(crate) fn remesh(
    input: &[Triangle],
    options: VoxelRemeshOptions,
) -> Result<(Vec<Triangle>, RemeshReport)> {
    let size = voxel_size(options.geometric_error)?;
    ensure!(!input.is_empty(), "cannot remesh an empty source mesh");
    let source_positions = input
        .iter()
        .flat_map(|triangle| triangle.positions)
        .collect::<Vec<_>>();
    let grid = VoxelGrid::centered(&source_positions, size)?;
    let extent = grid.cover(&source_positions)?;
    let source_indices = (0..source_positions.len())
        .map(u32::try_from)
        .collect::<Result<Vec<_>, _>>()?;
    let result = voxel_remesh::remesh(
        &source_indices,
        &source_positions,
        grid,
        extent,
        RemeshOptions::Solve,
    )?;
    drop(source_positions);
    drop(source_indices);
    let grid_dimensions = result.grid_dimensions;
    let (reconstructed_positions, reconstructed_indices) = result.into_mesh();
    let remeshed_triangles = reconstructed_indices.len() / 3;
    ensure!(
        !reconstructed_indices.is_empty(),
        "voxel remeshing produced empty geometry"
    );
    // Weld coincident solved vertices for simplification, without expanding the
    // indexed reconstruction into triangle corners.
    let (count, remap) =
        generate_vertex_remap(&reconstructed_positions, Some(&reconstructed_indices));
    let mut positions = remap_vertex_buffer(&reconstructed_positions, count, &remap);
    let mut indices = remap_index_buffer(
        Some(&reconstructed_indices),
        reconstructed_positions.len(),
        &remap,
    );
    drop(reconstructed_positions);
    drop(reconstructed_indices);
    drop(remap);
    let target = indices.len().min(12);
    let simplification_error = voxel_remesh::simplify_positions_with_update(
        &mut indices,
        &mut positions,
        None,
        target,
        options.geometric_error / 2.0,
        meshopt::SimplifyOptions::ErrorAbsolute,
    )?;
    let surface = SourceSurface::new(input);
    let mut triangles = Vec::with_capacity(indices.len() / 3);
    for triangle in indices.chunks_exact(3) {
        let positions = [
            positions[triangle[0] as usize],
            positions[triangle[1] as usize],
            positions[triangle[2] as usize],
        ];
        let feature_id = surface.owner(positions)?;
        triangles.push(Triangle {
            positions,
            feature_id,
        });
    }
    if triangles.is_empty() {
        bail!("simplification produced empty geometry");
    }
    let effective_voxel_size = f64::from(size);
    let remeshing_error_estimate = 3.0_f64.sqrt() * effective_voxel_size;
    let simplification_error = f64::from(simplification_error);
    let report = RemeshReport {
        requested_geometric_error: f64::from(options.geometric_error),
        grid_dimensions,
        effective_voxel_size,
        remeshing_error_estimate,
        simplification_error,
        encoding_error_estimate: 0.0,
        estimated_geometric_error: remeshing_error_estimate + simplification_error,
        input_triangles: input.len(),
        remeshed_triangles,
        output_triangles: triangles.len(),
    };
    Ok((triangles, report))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Purpose: check quality accounting. Every representable metric budget must
    /// produce its requested voxel size without extent-dependent clamping.
    #[test]
    fn budgets_choose_voxel_size_without_clamping() {
        for error in [0.1, 2.0, 40.0] {
            let size = voxel_size(error).unwrap();
            let diagonal = f64::from(size) * 3.0_f64.sqrt();
            assert!(diagonal <= f64::from(error) / 2.0);
            assert!(diagonal >= f64::from(error) / 2.0 * (1.0 - f64::from(f32::EPSILON)));
        }
        for error in [0.0, -1.0, f32::NAN, f32::INFINITY, f32::from_bits(1)] {
            assert!(voxel_size(error).is_err());
            assert!(remesh(
                &[],
                VoxelRemeshOptions {
                    geometric_error: error
                }
            )
            .is_err());
        }
    }

    /// Purpose: check categorical transfer. Two separated surfaces and coincident
    /// ties must choose the nearest triangle, breaking ties by stable source order.
    #[test]
    fn attribution_uses_triangle_distance_and_stable_ties() {
        let a = Triangle {
            positions: [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0]],
            feature_id: 7,
        };
        let b = Triangle {
            positions: a.positions.map(|p| [p[0], p[1], 5.0]),
            feature_id: 3,
        };
        let surface = SourceSurface::new(&[a, b]);
        assert_eq!(surface.owner(a.positions).unwrap(), 7);
        assert_eq!(surface.owner(b.positions).unwrap(), 3);
        let duplicate = Triangle { feature_id: 9, ..a };
        let surface = SourceSurface::new(&[a, duplicate]);
        assert_eq!(surface.owner(a.positions).unwrap(), 7);
        assert!(
            (surface.triangles[0].distance_squared(Point3::new(2.0, 2.0, 1.0)) - 1.0).abs()
                < f64::EPSILON
        );
        assert!(
            (surface.triangles[0].distance_squared(Point3::new(-1.0, 0.0, 0.0)) - 1.0).abs()
                < f64::EPSILON
        );
    }
}
