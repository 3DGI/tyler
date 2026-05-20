use anyhow::{bail, Context, Result};
use cityjson_lib::cityjson_types::prelude::CityObjectHandle;
use cityjson_lib::cityjson_types::v2_0::{GeometryType, VertexIndex};
use cityjson_lib::ops::Transformer;
use cityjson_lib::CityModel;
use earcutr::earcut;

use crate::GeographicClipRegion;

const CLIP_PLANE_EPSILON: f64 = 1.0e-12;
const GEOGRAPHIC_CLIP_INTERSECTION_ITERATIONS: usize = 64;

#[derive(Clone, Debug, Default)]
pub struct TriangleMeshOptions {
    pub clip_bbox: Option<[f64; 6]>,
    pub clip_geographic_region: Option<GeographicClipRegion>,
}

#[derive(Clone, Debug, Default)]
pub struct TriangleMesh {
    pub objects: Vec<ObjectMesh>,
}

#[derive(Clone, Debug)]
pub struct ObjectMesh {
    pub handle: CityObjectHandle,
    pub object_id: String,
    pub feature_type: String,
    pub triangles: Vec<SourceTriangle>,
}

#[derive(Clone, Copy, Debug)]
pub struct SourceTriangle {
    pub source_positions: [[f64; 3]; 3],
}

#[derive(Clone, Copy, Debug)]
struct ClipVertex {
    source_position: [f64; 3],
    clip_position: [f64; 3],
}

struct GeographicClipVolume {
    transformer: Option<CachedProjTransform>,
    west: f64,
    south: f64,
    east: f64,
    north: f64,
}

enum ClipVolume {
    SourceBbox([f64; 6]),
    GeographicRegion(GeographicClipVolume),
}

#[derive(Clone, Debug)]
struct CachedProjTransform {
    transformer: Transformer,
}

impl CachedProjTransform {
    fn new(source_crs: &str, target_crs: &'static str) -> Result<Self> {
        let transformer = cityjson_lib::ops::transformer(source_crs, target_crs)
            .with_context(|| format!("failed to create {source_crs} to {target_crs} transform"))?;
        Ok(Self { transformer })
    }

    fn convert(&self, point: [f64; 3]) -> Result<[f64; 3]> {
        self.transformer.transform(point).map_err(Into::into)
    }
}

pub fn build_triangle_mesh(
    model: &CityModel,
    options: &TriangleMeshOptions,
) -> Result<TriangleMesh> {
    let clip_volume = ClipVolume::from_options(options)?;
    let mut objects = Vec::new();

    for (handle, cityobject) in model.cityobjects().iter() {
        let Some(geometry_handles) = cityobject.geometry() else {
            continue;
        };

        let mut object_mesh = ObjectMesh {
            handle,
            object_id: cityobject.id().to_string(),
            feature_type: cityobject.type_cityobject().to_string(),
            triangles: Vec::new(),
        };

        for geometry_handle in geometry_handles {
            let geometry = model.resolve_geometry(*geometry_handle)?;
            let Some(boundary) = geometry.geometry().boundaries() else {
                continue;
            };
            match geometry.geometry().type_geometry() {
                GeometryType::MultiSurface | GeometryType::CompositeSurface => {
                    for surface in boundary.to_nested_multi_or_composite_surface()? {
                        add_surface(
                            &mut object_mesh.triangles,
                            &surface,
                            model,
                            clip_volume.as_ref(),
                        )?;
                    }
                }
                GeometryType::Solid => {
                    for shell in boundary.to_nested_solid()? {
                        for surface in shell {
                            add_surface(
                                &mut object_mesh.triangles,
                                &surface,
                                model,
                                clip_volume.as_ref(),
                            )?;
                        }
                    }
                }
                GeometryType::MultiSolid | GeometryType::CompositeSolid => {
                    for solid in boundary.to_nested_multi_or_composite_solid()? {
                        for shell in solid {
                            for surface in shell {
                                add_surface(
                                    &mut object_mesh.triangles,
                                    &surface,
                                    model,
                                    clip_volume.as_ref(),
                                )?;
                            }
                        }
                    }
                }
                _ => {}
            }
        }

        if !object_mesh.triangles.is_empty() {
            objects.push(object_mesh);
        }
    }

    Ok(TriangleMesh { objects })
}

fn add_surface(
    triangles: &mut Vec<SourceTriangle>,
    surface: &[Vec<u32>],
    model: &CityModel,
    clip_volume: Option<&ClipVolume>,
) -> Result<()> {
    if surface.is_empty() || surface[0].len() < 3 {
        return Ok(());
    }

    let exterior = &surface[0];
    let mut source_positions = Vec::new();
    let mut flat_coords = Vec::new();
    let mut hole_indices = Vec::new();
    let mut vertex_count = 0usize;

    for (ring_idx, ring) in surface.iter().enumerate() {
        if ring.len() < 3 {
            continue;
        }
        if ring_idx > 0 {
            hole_indices.push(vertex_count);
        }
        for &vertex_id in ring {
            let vertex = model
                .get_vertex(VertexIndex::new(vertex_id))
                .ok_or_else(|| anyhow::anyhow!("missing vertex {vertex_id}"))?;
            source_positions.push(vertex.to_array());
            vertex_count += 1;
        }
    }

    if source_positions.len() < 3 {
        return Ok(());
    }

    if surface.len() == 1 && exterior.len() == 3 {
        let source_triangle = [
            source_positions[0],
            source_positions[1],
            source_positions[2],
        ];
        push_triangle(triangles, source_triangle, clip_volume)?;
        return Ok(());
    }

    let (drop_axis, surface_normal) =
        projection_axis_and_normal(&source_positions[..exterior.len()]);
    for pos in &source_positions {
        match drop_axis {
            0 => {
                flat_coords.push(pos[1]);
                flat_coords.push(pos[2]);
            }
            1 => {
                flat_coords.push(pos[0]);
                flat_coords.push(pos[2]);
            }
            _ => {
                flat_coords.push(pos[0]);
                flat_coords.push(pos[1]);
            }
        }
    }

    let (flat_coords, hole_indices, source_index_map) =
        dedupe_polygon_rings(&flat_coords, &hole_indices);
    let triangulated =
        earcut(&flat_coords, &hole_indices, 2).context("Failed to triangulate surface")?;
    if triangulated.len() < 3 {
        return Ok(());
    }
    // earcutr normalizes the outer ring to positive 2D winding. Map that
    // known 2D orientation back to the selected 3D projection once per surface.
    let reverse_earcut_winding = should_reverse_earcut_winding(drop_axis, surface_normal);
    for tri in triangulated.chunks_exact(3) {
        let orig0 = source_index_map[tri[0]];
        let orig1 = source_index_map[tri[1]];
        let orig2 = source_index_map[tri[2]];
        let mut source_triangle = [
            source_positions[orig0],
            source_positions[orig1],
            source_positions[orig2],
        ];
        if reverse_earcut_winding {
            source_triangle.swap(1, 2);
        }
        push_triangle(triangles, source_triangle, clip_volume)?;
    }

    Ok(())
}

/// Strip consecutive bit-identical vertices from each ring before
/// triangulation. earcutr can spin forever on polygons that contain inner rings
/// whose vertices repeat.
#[allow(clippy::float_cmp)]
fn dedupe_polygon_rings(
    flat_coords: &[f64],
    hole_indices: &[usize],
) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
    let mut ring_offsets: Vec<usize> = std::iter::once(0)
        .chain(hole_indices.iter().copied())
        .collect();
    ring_offsets.push(flat_coords.len() / 2);

    let mut new_flat: Vec<f64> = Vec::with_capacity(flat_coords.len());
    let mut new_holes: Vec<usize> = Vec::with_capacity(hole_indices.len());
    let mut index_map: Vec<usize> = Vec::with_capacity(flat_coords.len() / 2);

    for ring_idx in 0..ring_offsets.len() - 1 {
        let start_vertex = ring_offsets[ring_idx];
        let end_vertex = ring_offsets[ring_idx + 1];
        if ring_idx > 0 {
            new_holes.push(new_flat.len() / 2);
        }

        let ring_start_in_new_flat = new_flat.len();
        for orig_idx in start_vertex..end_vertex {
            let x = flat_coords[orig_idx * 2];
            let y = flat_coords[orig_idx * 2 + 1];
            if new_flat.len() >= ring_start_in_new_flat + 2 {
                if let Some(&[px, py]) = new_flat.last_chunk::<2>() {
                    if px == x && py == y {
                        continue;
                    }
                }
            }

            new_flat.push(x);
            new_flat.push(y);
            index_map.push(orig_idx);
        }
    }

    (new_flat, new_holes, index_map)
}

fn should_reverse_earcut_winding(drop_axis: usize, surface_normal: Option<[f64; 3]>) -> bool {
    let Some(surface_normal) = surface_normal else {
        return false;
    };
    let positive_projected_winding_axis = match drop_axis {
        1 => -1.0,
        _ => 1.0,
    };
    positive_projected_winding_axis * surface_normal[drop_axis] < 0.0
}

fn push_triangle(
    triangles: &mut Vec<SourceTriangle>,
    source_positions: [[f64; 3]; 3],
    clip_volume: Option<&ClipVolume>,
) -> Result<()> {
    let clipped_triangles = if let Some(clip_volume) = clip_volume {
        clip_triangle_to_volume(source_positions, clip_volume)?
    } else {
        vec![source_positions]
    };

    for triangle in clipped_triangles {
        if !is_degenerate_source_triangle(triangle) {
            triangles.push(SourceTriangle {
                source_positions: triangle,
            });
        }
    }
    Ok(())
}

impl ClipVolume {
    fn from_options(options: &TriangleMeshOptions) -> Result<Option<Self>> {
        if let Some(region) = &options.clip_geographic_region {
            return Self::geographic_region(region).map(Some);
        }
        Ok(options.clip_bbox.map(Self::SourceBbox))
    }

    fn geographic_region(region: &GeographicClipRegion) -> Result<Self> {
        let source_crs = canonical_epsg_crs(&region.source_crs)?;
        let transformer = if source_crs == "EPSG:4979" {
            None
        } else {
            Some(CachedProjTransform::new(&source_crs, "EPSG:4979")?)
        };
        Ok(Self::GeographicRegion(GeographicClipVolume {
            transformer,
            west: region.west,
            south: region.south,
            east: region.east,
            north: region.north,
        }))
    }

    fn clip_position(&self, source_position: [f64; 3]) -> Result<[f64; 3]> {
        match self {
            Self::SourceBbox(_) => Ok(source_position),
            Self::GeographicRegion(region) => {
                let geographic = if let Some(transformer) = &region.transformer {
                    transformer
                        .convert(source_position)
                        .context("failed to project position to EPSG:4979 for clipping")?
                } else {
                    source_position
                };
                Ok(geographic)
            }
        }
    }

    fn planes(&self) -> ([(usize, f64, bool); 6], usize) {
        match self {
            Self::SourceBbox(bbox) => (
                [
                    (0, bbox[0], false),
                    (0, bbox[3], true),
                    (1, bbox[1], false),
                    (1, bbox[4], true),
                    (2, bbox[2], false),
                    (2, bbox[5], true),
                ],
                6,
            ),
            Self::GeographicRegion(region) => (
                [
                    (0, region.west, false),
                    (0, region.east, true),
                    (1, region.south, false),
                    (1, region.north, true),
                    (0, 0.0, true),
                    (0, 0.0, true),
                ],
                4,
            ),
        }
    }

    fn intersect_edge(
        &self,
        start: ClipVertex,
        end: ClipVertex,
        axis: usize,
        boundary: f64,
    ) -> Result<Option<ClipVertex>> {
        match self {
            Self::SourceBbox(_) => self.intersect_edge_linear(start, end, axis, boundary),
            Self::GeographicRegion(_) => self.intersect_edge_geographic(start, end, axis, boundary),
        }
    }

    fn intersect_edge_linear(
        &self,
        start: ClipVertex,
        end: ClipVertex,
        axis: usize,
        boundary: f64,
    ) -> Result<Option<ClipVertex>> {
        let delta = end.clip_position[axis] - start.clip_position[axis];
        if delta.abs() <= CLIP_PLANE_EPSILON {
            return Ok(None);
        }
        let t = (boundary - start.clip_position[axis]) / delta;
        self.interpolate_clip_vertex(start, end, t).map(Some)
    }

    fn intersect_edge_geographic(
        &self,
        start: ClipVertex,
        end: ClipVertex,
        axis: usize,
        boundary: f64,
    ) -> Result<Option<ClipVertex>> {
        let mut low_t = 0.0;
        let mut high_t = 1.0;
        let mut low_distance = start.clip_position[axis] - boundary;
        let high_distance = end.clip_position[axis] - boundary;

        if low_distance.abs() <= CLIP_PLANE_EPSILON {
            return Ok(Some(start));
        }
        if high_distance.abs() <= CLIP_PLANE_EPSILON {
            return Ok(Some(end));
        }
        if low_distance.signum() == high_distance.signum() {
            bail!(
                "geographic clip edge does not bracket boundary on axis {axis}: {low_distance} to {high_distance}"
            );
        }

        for _ in 0..GEOGRAPHIC_CLIP_INTERSECTION_ITERATIONS {
            let mid_t = f64::midpoint(low_t, high_t);
            let midpoint = self.interpolate_clip_vertex(start, end, mid_t)?;
            let mid_distance = midpoint.clip_position[axis] - boundary;
            if mid_distance.abs() <= CLIP_PLANE_EPSILON || (high_t - low_t).abs() <= f64::EPSILON {
                return Ok(Some(midpoint));
            }
            if low_distance.signum() == mid_distance.signum() {
                low_t = mid_t;
                low_distance = mid_distance;
            } else {
                high_t = mid_t;
            }
        }

        let t = f64::midpoint(low_t, high_t);
        self.interpolate_clip_vertex(start, end, t).map(Some)
    }

    fn interpolate_clip_vertex(
        &self,
        start: ClipVertex,
        end: ClipVertex,
        t: f64,
    ) -> Result<ClipVertex> {
        let source_position = [
            start.source_position[0] + (end.source_position[0] - start.source_position[0]) * t,
            start.source_position[1] + (end.source_position[1] - start.source_position[1]) * t,
            start.source_position[2] + (end.source_position[2] - start.source_position[2]) * t,
        ];
        Ok(ClipVertex {
            source_position,
            clip_position: self.clip_position(source_position)?,
        })
    }
}

fn clip_triangle_to_volume(
    triangle: [[f64; 3]; 3],
    clip_volume: &ClipVolume,
) -> Result<Vec<[[f64; 3]; 3]>> {
    let mut polygon = vec![
        ClipVertex {
            source_position: triangle[0],
            clip_position: clip_volume.clip_position(triangle[0])?,
        },
        ClipVertex {
            source_position: triangle[1],
            clip_position: clip_volume.clip_position(triangle[1])?,
        },
        ClipVertex {
            source_position: triangle[2],
            clip_position: clip_volume.clip_position(triangle[2])?,
        },
    ];

    let (planes, plane_count) = clip_volume.planes();
    for (axis, boundary, keep_less_equal) in planes.into_iter().take(plane_count) {
        polygon =
            clip_polygon_against_plane(polygon, clip_volume, axis, boundary, keep_less_equal)?;
        if polygon.len() < 3 {
            return Ok(Vec::new());
        }
    }

    let mut triangles = Vec::with_capacity(polygon.len().saturating_sub(2));
    for index in 1..polygon.len() - 1 {
        triangles.push([
            polygon[0].source_position,
            polygon[index].source_position,
            polygon[index + 1].source_position,
        ]);
    }
    Ok(triangles)
}

fn clip_polygon_against_plane(
    polygon: Vec<ClipVertex>,
    clip_volume: &ClipVolume,
    axis: usize,
    boundary: f64,
    keep_less_equal: bool,
) -> Result<Vec<ClipVertex>> {
    if polygon.is_empty() {
        return Ok(polygon);
    }

    let mut clipped = Vec::new();
    for index in 0..polygon.len() {
        let current = polygon[index];
        let next = polygon[(index + 1) % polygon.len()];
        let current_distance = current.clip_position[axis] - boundary;
        let next_distance = next.clip_position[axis] - boundary;
        let current_inside = if keep_less_equal {
            current_distance <= CLIP_PLANE_EPSILON
        } else {
            current_distance >= -CLIP_PLANE_EPSILON
        };
        let next_inside = if keep_less_equal {
            next_distance <= CLIP_PLANE_EPSILON
        } else {
            next_distance >= -CLIP_PLANE_EPSILON
        };

        if current_inside && next_inside {
            clipped.push(next);
            continue;
        }
        if current_inside != next_inside {
            if let Some(intersection) = clip_volume.intersect_edge(current, next, axis, boundary)? {
                clipped.push(intersection);
            }
        }
        if !current_inside && next_inside {
            clipped.push(next);
        }
    }

    Ok(clipped)
}

fn canonical_epsg_crs(value: &str) -> Result<String> {
    if let Some(code) = value.strip_prefix("EPSG:") {
        let parsed = code.parse::<u32>().context("invalid EPSG code")?;
        return Ok(format!("EPSG:{parsed}"));
    }

    let code = value
        .rsplit(['/', ':'])
        .find(|part| !part.is_empty())
        .ok_or_else(|| anyhow::anyhow!("could not extract EPSG code from {value}"))?;
    let parsed = code
        .parse::<u32>()
        .with_context(|| format!("invalid EPSG code in {value}"))?;
    Ok(format!("EPSG:{parsed}"))
}

fn projection_axis_and_normal(positions: &[[f64; 3]]) -> (usize, Option<[f64; 3]>) {
    let normal = compute_polygon_normal(positions);
    if let Some(normal) = normal {
        return (dominant_axis(normal), Some(normal));
    }

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for pos in positions {
        for (axis, coordinate) in pos.iter().enumerate() {
            min[axis] = min[axis].min(*coordinate);
            max[axis] = max[axis].max(*coordinate);
        }
    }

    let drop_axis = (0..3)
        .min_by(|&lhs, &rhs| {
            (max[lhs] - min[lhs])
                .partial_cmp(&(max[rhs] - min[rhs]))
                .unwrap()
        })
        .unwrap_or(2);
    (drop_axis, None)
}

fn compute_polygon_normal(positions: &[[f64; 3]]) -> Option<[f64; 3]> {
    if positions.len() < 3 {
        return None;
    }

    let mut normal = [0.0_f64; 3];
    for (current, next) in positions
        .iter()
        .zip(positions.iter().cycle().skip(1))
        .take(positions.len())
    {
        normal[0] += (current[1] - next[1]) * (current[2] + next[2]);
        normal[1] += (current[2] - next[2]) * (current[0] + next[0]);
        normal[2] += (current[0] - next[0]) * (current[1] + next[1]);
    }

    let length = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    (length > f64::EPSILON).then_some([normal[0] / length, normal[1] / length, normal[2] / length])
}

fn dominant_axis(normal: [f64; 3]) -> usize {
    (0..3)
        .max_by(|&lhs, &rhs| normal[lhs].abs().partial_cmp(&normal[rhs].abs()).unwrap())
        .unwrap_or(2)
}

fn triangle_normal(points: [[f64; 3]; 3]) -> [f64; 3] {
    let u = [
        points[1][0] - points[0][0],
        points[1][1] - points[0][1],
        points[1][2] - points[0][2],
    ];
    let v = [
        points[2][0] - points[0][0],
        points[2][1] - points[0][1],
        points[2][2] - points[0][2],
    ];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn is_degenerate_source_triangle(points: [[f64; 3]; 3]) -> bool {
    let cross = triangle_normal(points);
    let area_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
    area_sq <= 1.0e-18
}
