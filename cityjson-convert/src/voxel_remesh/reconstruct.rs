//! Surface reconstruction ported from meshoptimizer v1.3 (MIT; see LICENSE.md).
#![forbid(unsafe_code)]

use super::grid::{Layout, Occupancy, SetBits};
use super::tables::TABLES;
use super::{allocate, RemeshOptions, RemeshResult, VoxelExtent, VoxelGrid};
use anyhow::{anyhow, Result};
use num_traits::ToPrimitive;

#[derive(Clone, Copy, Debug, Default)]
struct Voxel {
    point: [f32; 3],
    weight: f32,
    diagonal: [f32; 3],
    off_diagonal: [f32; 3],
    linear: [f32; 3],
    constant: f32,
    octants: u8,
}

impl Voxel {
    fn accumulate(
        &mut self,
        point: [f32; 3],
        octant: u8,
        weight: f32,
        normal: [f32; 3],
        solve: bool,
    ) {
        for (sum, value) in self.point.iter_mut().zip(point) {
            *sum += value * weight;
        }
        self.weight += weight;
        self.octants |= octant;
        if !solve {
            return;
        }
        let [a, b, c] = normal;
        let d = -(a * point[0] + b * point[1] + c * point[2]);
        let [aw, bw, cw, dw] = [a * weight, b * weight, c * weight, d * weight];
        self.diagonal[0] += a * aw;
        self.diagonal[1] += b * bw;
        self.diagonal[2] += c * cw;
        self.off_diagonal[0] += a * bw;
        self.off_diagonal[1] += a * cw;
        self.off_diagonal[2] += b * cw;
        self.linear[0] += a * dw;
        self.linear[1] += b * dw;
        self.linear[2] += c * dw;
        self.constant += d * dw;
    }

    fn minimum(&self) -> Option<[f32; 3]> {
        let lambda = 3e-2;
        let rw = lambda * self.weight;
        let [a00, a11, a22] = self.diagonal.map(|value| value + rw);
        let [a10, a20, a21] = self.off_diagonal;
        let rhs =
            std::array::from_fn::<_, 3, _>(|axis| lambda * self.point[axis] - self.linear[axis]);
        let epsilon = 1e-6 * self.weight;
        // LDL^T decomposition of the regularized quadric.
        let d0 = a00;
        let l10 = a10 / d0;
        let l20 = a20 / d0;
        let d1 = a11 - a10 * l10;
        let dl21 = a21 - a20 * l10;
        let l21 = dl21 / d1;
        let d2 = a22 - a20 * l20 - dl21 * l21;
        if d0.abs() <= epsilon || d1.abs() <= epsilon || d2.abs() <= epsilon {
            return None;
        }
        let y0 = rhs[0];
        let y1 = rhs[1] - l10 * y0;
        let y2 = rhs[2] - l20 * y0 - l21 * y1;
        let z0 = y0 / d0;
        let z1 = y1 / d1;
        let z2 = y2 / d2;
        let rz = z2;
        let ry = z1 - l21 * rz;
        Some([z0 - l10 * ry - l20 * rz, ry, rz])
    }

    fn solve(&mut self, cell: [i64; 3], scale: f32, enabled: bool) -> Result<()> {
        let mut point = self.point.map(|value| value / self.weight);
        if enabled {
            if let Some(solution) = self.minimum() {
                let distance = (solution[0] - point[0]).powi(2)
                    + (solution[1] - point[1]).powi(2)
                    + (solution[2] - point[2]).powi(2);
                if distance < 3.0 / (scale * scale) {
                    let size = 1.0 / scale;
                    for axis in 0..3 {
                        let corner = cell[axis].to_f32().expect("cell coordinate fits f32") * size;
                        point[axis] = solution[axis].clamp(corner, corner + size);
                    }
                }
            }
        }
        if point.iter().any(|value| !value.is_finite()) {
            return Err(anyhow!(
                "voxel reconstruction produced non-finite positions"
            ));
        }
        self.point = point;
        Ok(())
    }

    fn error(&self, point: [f32; 3]) -> f32 {
        let [x, y, z] = point;
        let [a00, a11, a22] = self.diagonal;
        let [a10, a20, a21] = self.off_diagonal;
        let [b0, b1, b2] = self.linear;
        let rx = (b0 + a10 * y) * 2.0 + a00 * x;
        let ry = (b1 + a21 * z) * 2.0 + a11 * y;
        let rz = (b2 + a20 * x) * 2.0 + a22 * z;
        (self.constant + rx * x + ry * y + rz * z).abs()
    }
}

/// Sampling depends on the original triangle and lattice, never the processing extent.
fn sample_triangles(
    indices: &[u32],
    positions: &[[f32; 3]],
    frame: VoxelGrid,
    layout: Layout,
    mut sample: impl FnMut(usize, [f32; 3], u8, f32, [f32; 3]),
) -> Result<()> {
    let scale = 1.0 / frame.voxel_size;
    for triangle in indices.chunks_exact(3) {
        let va = positions[triangle[0] as usize];
        let vb = positions[triangle[1] as usize];
        let vc = positions[triangle[2] as usize];
        let edge = std::array::from_fn::<_, 3, _>(|axis| vb[axis] - va[axis]);
        let side = std::array::from_fn::<_, 3, _>(|axis| vc[axis] - va[axis]);
        let other = std::array::from_fn::<_, 3, _>(|axis| vc[axis] - vb[axis]);
        let length_squared =
            |value: [f32; 3]| value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
        let max_edge = length_squared(edge)
            .max(length_squared(side))
            .max(length_squared(other))
            .sqrt();
        let normal = [
            edge[1] * side[2] - edge[2] * side[1],
            edge[2] * side[0] - edge[0] * side[2],
            edge[0] * side[1] - edge[1] * side[0],
        ];
        let area = length_squared(normal).sqrt();
        if area == 0.0 {
            continue;
        }
        let count = max_edge * scale * 2.0;
        if !area.is_finite() || !count.is_finite() || f64::from(count) > f64::from(u32::MAX - 2) {
            return Err(anyhow!(
                "triangle or voxel size is unsuitable for float sampling"
            ));
        }
        let samples = count
            .to_u32()
            .expect("validated sample count fits u32")
            .max(1);
        let reciprocal = 1.0 / samples.to_f32().expect("sample count fits f32");
        let weight = area
            / ((samples + 1).to_f32().expect("sample count fits f32")
                * (samples + 2).to_f32().expect("sample count fits f32"));
        let normal = normal.map(|value| value * (1.0 / area));
        let start = frame.relative(va)?;
        for u in 0..=samples {
            for v in 0..=samples - u {
                let su = u.to_f32().expect("sample index fits f32") * reciprocal;
                let sv = v.to_f32().expect("sample index fits f32") * reciprocal;
                let point =
                    std::array::from_fn(|axis| start[axis] + su * edge[axis] + sv * side[axis]);
                let half = frame.half_cell(point)?;
                let cell = half.map(|value| value.div_euclid(2));
                if let Some(flat) = layout.flat(cell) {
                    let octant =
                        u8::try_from((half[0] & 1) | ((half[1] & 1) << 1) | ((half[2] & 1) << 2))
                            .expect("three octant bits fit u8");
                    sample(flat, point, 1 << octant, weight, normal);
                }
            }
        }
    }
    Ok(())
}

fn octant_decider(corners: &[Option<u32>; 8], voxels: &[Voxel]) -> usize {
    for (corner, vertex) in corners.iter().enumerate() {
        if let Some(vertex) = vertex {
            if voxels[*vertex as usize].octants & (1 << (7 - corner)) != 0 {
                return 0;
            }
        }
    }
    1
}

fn quadric_decider(cube: usize, corners: &[Option<u32>; 8], voxels: &[Voxel], size: f32) -> usize {
    let primary = TABLES.triangles[cube][0];
    let quad = (primary[0] << 4) | (primary[1] & 15);
    let vertices: [_; 4] = std::array::from_fn(|index| {
        let corner = usize::from((quad >> (12 - index * 4)) & 15);
        &voxels[corners[corner].expect("triangulation only references surface corners") as usize]
    });
    let [a, b, c, d] = vertices;
    // Preserve summation order on both sides of thin sheets.
    let combine = |av, bv, cv, dv| (av + dv) + (bv + cv);
    let sum = Voxel {
        diagonal: std::array::from_fn(|i| {
            combine(a.diagonal[i], b.diagonal[i], c.diagonal[i], d.diagonal[i])
        }),
        off_diagonal: std::array::from_fn(|i| {
            combine(
                a.off_diagonal[i],
                b.off_diagonal[i],
                c.off_diagonal[i],
                d.off_diagonal[i],
            )
        }),
        linear: std::array::from_fn(|i| {
            combine(a.linear[i], b.linear[i], c.linear[i], d.linear[i])
        }),
        constant: combine(a.constant, b.constant, c.constant, d.constant),
        weight: combine(a.weight, b.weight, c.weight, d.weight),
        ..Voxel::default()
    };
    let midpoint = |left: &Voxel, right: &Voxel| {
        std::array::from_fn(|i| (left.point[i] + right.point[i]) * 0.5)
    };
    let error0 = sum.error(midpoint(b, c));
    let error1 = sum.error(midpoint(a, d));
    usize::from(error0 > sum.weight * (0.1 * 0.1 * size * size) && error1 < error0 * (0.85 * 0.85))
}

/// Skip uniform 64-cell spans and enumerate only cubes crossing the surface.
fn visit_cubes(
    grid: &Occupancy,
    mut visit: impl FnMut(usize, [usize; 8]) -> Result<()>,
) -> Result<()> {
    let layout = grid.layout;
    let stride = layout.row_words * 64;
    for z in 0..layout.dimensions[2] - 1 {
        for y in 0..layout.dimensions[1] - 1 {
            let row = y + layout.dimensions[1] * z;
            let rows = [
                row,
                row + 1,
                row + layout.dimensions[1],
                row + layout.dimensions[1] + 1,
            ];
            if rows.iter().all(|&row| !grid.has_surface(row)) {
                continue;
            }
            for word in 0..layout.row_words {
                let first = word * 64;
                if first >= layout.dimensions[0] - 1 {
                    break;
                }
                let left = rows.map(|row| grid.occupied(row * layout.row_words + word));
                let right = std::array::from_fn::<_, 4, _>(|i| {
                    let next = if word + 1 < layout.row_words {
                        grid.occupied(rows[i] * layout.row_words + word + 1)
                    } else {
                        0
                    };
                    (left[i] >> 1) | ((next & 1) << 63)
                });
                let any_left = left.iter().fold(0, |a, b| a | b);
                let any_right = right.iter().fold(0, |a, b| a | b);
                let full_left = left.iter().fold(u64::MAX, |a, b| a & b);
                let full_right = right.iter().fold(u64::MAX, |a, b| a & b);
                let count = (layout.dimensions[0] - 1 - first).min(64);
                let candidates =
                    (any_left | any_right) & !(full_left & full_right) & (u64::MAX >> (64 - count));
                for bit in SetBits(candidates) {
                    let cube = (0..4).fold(0_usize, |cube, i| {
                        cube | ((((left[i] >> bit) & 1) as usize) << (2 * i))
                            | ((((right[i] >> bit) & 1) as usize) << (2 * i + 1))
                    });
                    let x = first + bit;
                    let corners =
                        std::array::from_fn(|corner| x + (corner & 1) + stride * rows[corner >> 1]);
                    visit(cube, corners)?;
                }
            }
        }
    }
    Ok(())
}

fn polygonize(
    grid: &Occupancy,
    voxels: &[Voxel],
    options: RemeshOptions,
    scale: f32,
) -> Result<Vec<u32>> {
    let mut capacity = 0_usize;
    visit_cubes(grid, |cube, _| {
        capacity = capacity
            .checked_add(usize::from(TABLES.counts[cube]) * 3)
            .ok_or_else(|| anyhow!("triangle index count overflow"))?;
        Ok(())
    })?;
    let mut indices = Vec::new();
    indices
        .try_reserve_exact(capacity)
        .map_err(|error| anyhow!(error.to_string()))?;
    visit_cubes(grid, |cube, cells| {
        let corners = cells.map(|cell| grid.vertex(cell));
        let variant = match TABLES.modes[cube] {
            1 => octant_decider(&corners, voxels),
            2 if options.contains(RemeshOptions::Solve) => {
                quadric_decider(cube, &corners, voxels, 1.0 / scale)
            }
            _ => 0,
        };
        for &triangle in &TABLES.triangles[cube][variant] {
            if triangle == 0 {
                break;
            }
            for shift in [8, 4, 0] {
                indices.push(
                    corners[usize::from((triangle >> shift) & 15)]
                        .ok_or_else(|| anyhow!("triangulation referenced an interior cell"))?,
                );
            }
        }
        Ok(())
    })?;
    Ok(indices)
}

pub(super) fn reconstruct(
    indices: &[u32],
    positions: &[[f32; 3]],
    frame: VoxelGrid,
    extent: VoxelExtent,
    options: RemeshOptions,
) -> Result<RemeshResult> {
    let layout = Layout::new(extent)?;
    if !options.contains(RemeshOptions::Shell) && !positions.is_empty() {
        let required = frame.cover(positions)?;
        if (0..3).any(|axis| {
            extent.min[axis] > required.min[axis] || extent.max[axis] < required.max[axis]
        }) {
            return Err(anyhow!("solid reconstruction requires a padded extent covering the full input; cropped domains require shell mode until interior boundary classification is supported"));
        }
    }
    if indices.is_empty() {
        return Ok(RemeshResult {
            positions: Vec::new(),
            indices: Vec::new(),
            grid_dimensions: layout
                .dimensions
                .map(|value| u32::try_from(value).expect("grid axis fits u32")),
            grid: frame,
            extent,
            vertex_cells: Vec::new(),
        });
    }
    let mut grid = Occupancy::new(layout)?;
    sample_triangles(indices, positions, frame, layout, |flat, _, _, _, _| {
        grid.surface[flat / 64] |= 1 << (flat % 64);
    })?;
    let count = grid.pack()?;
    let mut voxels = allocate(count, Voxel::default())?;
    let solve = options.contains(RemeshOptions::Solve);
    sample_triangles(
        indices,
        positions,
        frame,
        layout,
        |flat, point, octant, weight, normal| {
            let vertex = grid
                .vertex(flat)
                .expect("sampling passes visit identical cells");
            voxels[vertex as usize].accumulate(point, octant, weight, normal, solve);
        },
    )?;
    if !options.contains(RemeshOptions::Shell) {
        grid.solidify()?;
    }
    let scale = 1.0 / frame.voxel_size;
    for (voxel, flat) in voxels.iter_mut().zip(grid.cells()) {
        voxel.solve(layout.cell(flat), scale, solve)?;
    }
    let indices = polygonize(&grid, &voxels, options, scale)?;
    // Occupancy/rank and quadrics are no longer needed when materializing positions.
    let mut vertex_cells = Vec::new();
    vertex_cells
        .try_reserve_exact(count)
        .map_err(|error| anyhow!(error.to_string()))?;
    vertex_cells.extend(grid.cells());
    drop(grid);
    let mut positions = Vec::new();
    positions
        .try_reserve_exact(count)
        .map_err(|error| anyhow!(error.to_string()))?;
    for voxel in voxels {
        let position = std::array::from_fn(|axis| {
            (f64::from(voxel.point[axis]) + frame.origin[axis])
                .to_f32()
                .unwrap_or(f32::INFINITY)
        });
        if position.iter().any(|value| !value.is_finite()) {
            return Err(anyhow!("output positions cannot be represented as float"));
        }
        positions.push(position);
    }
    Ok(RemeshResult {
        positions,
        indices,
        grid_dimensions: layout
            .dimensions
            .map(|value| u32::try_from(value).expect("grid axis fits u32")),
        grid: frame,
        extent,
        vertex_cells,
    })
}
