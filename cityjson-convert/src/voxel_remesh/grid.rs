//! Grid coordinates and compact occupancy, derived from meshoptimizer v1.3.
#![forbid(unsafe_code)]

use super::allocate;
use anyhow::{anyhow, Result};
use num_traits::ToPrimitive;

/// A voxel lattice in the same Cartesian frame as the supplied positions.
/// Reuse the same origin and size for all regions that must share cell identities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VoxelGrid {
    pub origin: [f64; 3],
    pub voxel_size: f32,
}

/// Half-open global cell range, including an empty padding cell on each side.
/// This is a processing domain, not a clipping or ownership boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct VoxelExtent {
    pub min: [i64; 3],
    pub max: [i64; 3],
}

impl VoxelGrid {
    pub(super) fn validate(self) -> Result<()> {
        if self.origin.iter().any(|value| !value.is_finite())
            || !self.voxel_size.is_finite()
            || self.voxel_size <= 0.0
            || !(2.0 / self.voxel_size).is_finite()
        {
            return Err(anyhow!(
                "grid origin and positive voxel size must be finite"
            ));
        }
        Ok(())
    }

    /// Creates a centred lattice for a standalone whole mesh.
    /// Tiled callers should instead supply an explicit shared origin.
    ///
    /// # Errors
    /// Rejects invalid positions, voxel sizes and unrepresentable dimensions.
    pub fn centered(positions: &[[f32; 3]], voxel_size: f32) -> Result<Self> {
        let mut grid = Self {
            origin: [0.0; 3],
            voxel_size,
        };
        grid.validate()?;
        if positions.is_empty() {
            return Err(anyhow!("cannot centre an empty mesh"));
        }
        let mut min = [f32::INFINITY; 3];
        let mut max = [f32::NEG_INFINITY; 3];
        for position in positions {
            for axis in 0..3 {
                if !position[axis].is_finite() {
                    return Err(anyhow!("positions must be finite"));
                }
                min[axis] = min[axis].min(position[axis]);
                max[axis] = max[axis].max(position[axis]);
            }
        }
        for axis in 0..3 {
            let extent = f64::from(max[axis]) - f64::from(min[axis]);
            let dimension = (extent / f64::from(voxel_size)).ceil() + 3.0;
            if dimension > f64::from(1 << 23) {
                return Err(anyhow!("voxel grid exceeds float coordinate precision"));
            }
            let width = (dimension.max(4.0) - 2.0) * f64::from(voxel_size);
            // Keep standalone grid placement at float precision, like the source algorithm.
            grid.origin[axis] = f64::from(
                (f64::from(min[axis]) - (width - extent) * 0.5)
                    .to_f32()
                    .unwrap_or(f32::INFINITY),
            );
        }
        grid.validate()?;
        Ok(grid)
    }

    /// Returns a padded processing extent covering the positions on this lattice.
    ///
    /// # Errors
    /// Rejects empty/non-finite positions and coordinates without half-cell precision.
    pub fn cover(self, positions: &[[f32; 3]]) -> Result<VoxelExtent> {
        self.validate()?;
        if positions.is_empty() {
            return Err(anyhow!("cannot cover an empty mesh"));
        }
        let mut min = [i64::MAX; 3];
        let mut max = [i64::MIN; 3];
        for &position in positions {
            let relative = self.relative(position)?;
            let half = self.half_cell(relative)?;
            for axis in 0..3 {
                let cell = half[axis].div_euclid(2);
                min[axis] = min[axis].min(cell - 1);
                max[axis] = max[axis].max(cell + 2);
            }
        }
        for axis in 0..3 {
            max[axis] = max[axis].max(min[axis] + 4);
        }
        let extent = VoxelExtent { min, max };
        Layout::new(extent)?;
        Ok(extent)
    }

    pub(super) fn relative(self, position: [f32; 3]) -> Result<[f32; 3]> {
        let relative = std::array::from_fn(|axis| {
            (f64::from(position[axis]) - self.origin[axis])
                .to_f32()
                .unwrap_or(f32::INFINITY)
        });
        if relative.iter().any(|value| !value.is_finite()) {
            return Err(anyhow!("positions must be finite in the grid frame"));
        }
        Ok(relative)
    }

    #[inline]
    pub(super) fn half_cell(self, position: [f32; 3]) -> Result<[i64; 3]> {
        let scaled = position.map(|value| value * (2.0 / self.voxel_size));
        if scaled
            .iter()
            .any(|value| !value.is_finite() || value.abs() > 16_777_216.0)
        {
            return Err(anyhow!(
                "grid origin is too far from positions to retain half-cell precision"
            ));
        }
        Ok(scaled.map(|value| {
            // The checked half-cell range fits in i32. Most standalone grids
            // have positive coordinates; only negative fractions need flooring.
            let integer = value.to_i32().expect("validated half-cell range fits i32");
            let truncated = integer.to_f32().expect("i32 fits f32");
            i64::from(integer - i32::from(value < 0.0 && value < truncated))
        }))
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Layout {
    pub extent: VoxelExtent,
    pub dimensions: [usize; 3],
    pub row_words: usize,
    pub rows: usize,
    pub words: usize,
}

impl Layout {
    pub fn new(extent: VoxelExtent) -> Result<Self> {
        let mut dimensions = [0; 3];
        for (axis, dimension) in dimensions.iter_mut().enumerate() {
            let length = extent.max[axis]
                .checked_sub(extent.min[axis])
                .ok_or_else(|| anyhow!("cell extent overflows"))?;
            if !(4..=1 << 23).contains(&length)
                || extent.min[axis].abs_diff(0) > 1 << 23
                || extent.max[axis].abs_diff(0) > (1 << 23) + 3
            {
                return Err(anyhow!(
                    "processing extent must have 4..2^23 cells per axis within float precision"
                ));
            }
            *dimension = usize::try_from(length)
                .map_err(|error| anyhow!(format!("cell extent cannot be addressed: {error}")))?;
        }
        let row_words = dimensions[0].div_ceil(64);
        let rows = dimensions[1]
            .checked_mul(dimensions[2])
            .ok_or_else(|| anyhow!("grid row count overflow"))?;
        let words = row_words
            .checked_mul(rows)
            .filter(|&count| count < (isize::MAX as usize) / 8)
            .ok_or_else(|| anyhow!("grid byte count overflow"))?;
        Ok(Self {
            extent,
            dimensions,
            row_words,
            rows,
            words,
        })
    }

    pub fn flat(self, cell: [i64; 3]) -> Option<usize> {
        let mut local = [0; 3];
        for axis in 0..3 {
            let coordinate =
                usize::try_from(cell[axis].checked_sub(self.extent.min[axis])?).ok()?;
            // The outer layer stays empty. Samples outside are skipped, never clamped to it.
            if coordinate < 1 || coordinate >= self.dimensions[axis] - 1 {
                return None;
            }
            local[axis] = coordinate;
        }
        Some(local[0] + self.row_words * 64 * (local[1] + self.dimensions[1] * local[2]))
    }

    pub fn cell(self, flat: usize) -> [i64; 3] {
        let stride = self.row_words * 64;
        let row = flat / stride;
        [
            self.extent.min[0] + i64::try_from(flat % stride).expect("grid X fits i64"),
            self.extent.min[1] + i64::try_from(row % self.dimensions[1]).expect("grid Y fits i64"),
            self.extent.min[2] + i64::try_from(row / self.dimensions[1]).expect("grid Z fits i64"),
        ]
    }

    fn mask(self, word: usize) -> u64 {
        let first = word * 64;
        let end = self.dimensions[0] - 1;
        if first >= end {
            return 0;
        }
        let count = (end - first).min(64);
        let mask = u64::MAX >> (64 - count);
        if first == 0 {
            mask & !1
        } else {
            mask
        }
    }
}

pub(super) struct Occupancy {
    pub layout: Layout,
    pub surface: Vec<u64>,
    pub interior: Vec<u64>,
    pub rank: Vec<u32>,
}

impl Occupancy {
    pub fn new(layout: Layout) -> Result<Self> {
        Ok(Self {
            layout,
            surface: allocate(layout.words, 0)?,
            interior: Vec::new(),
            rank: Vec::new(),
        })
    }

    pub fn pack(&mut self) -> Result<usize> {
        self.rank = allocate(self.layout.words + 1, 0)?;
        let mut count = 0_u32;
        for (word, &bits) in self.surface.iter().enumerate() {
            self.rank[word] = count;
            count = count
                .checked_add(bits.count_ones())
                .ok_or_else(|| anyhow!("surface vertex count exceeds u32 indices"))?;
        }
        self.rank[self.layout.words] = count;
        Ok(count as usize)
    }

    pub fn vertex(&self, flat: usize) -> Option<u32> {
        let word = flat / 64;
        let bit = flat % 64;
        let mask = 1_u64 << bit;
        if self.surface[word] & mask == 0 {
            return None;
        }
        Some(self.rank[word] + (self.surface[word] & (mask - 1)).count_ones())
    }

    pub fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        self.surface
            .iter()
            .enumerate()
            .flat_map(|(word, &bits)| SetBits(bits).map(move |bit| word * 64 + bit))
    }

    pub fn occupied(&self, word: usize) -> u64 {
        self.surface[word] | self.interior.get(word).copied().unwrap_or(0)
    }

    pub fn has_surface(&self, row: usize) -> bool {
        let start = row * self.layout.row_words;
        self.rank[start] != self.rank[start + self.layout.row_words]
    }

    /// The original row flood fill, operating on 64 cells at a time.
    pub fn solidify(&mut self) -> Result<()> {
        let layout = self.layout;
        self.interior = allocate(layout.words, 0)?;
        for z in 1..layout.dimensions[2] - 1 {
            for y in 1..layout.dimensions[1] - 1 {
                let row = y + layout.dimensions[1] * z;
                if !self.has_surface(row) {
                    continue;
                }
                for x in 0..layout.row_words {
                    let word = row * layout.row_words + x;
                    self.interior[word] = !self.surface[word] & layout.mask(x);
                }
            }
        }
        let mut queue = allocate(layout.rows, 0_usize)?;
        for (row, entry) in queue.iter_mut().enumerate() {
            *entry = row;
        }
        let mut queued = allocate(layout.rows, true)?;
        while let Some(row) = queue.pop() {
            queued[row] = false;
            if self.has_surface(row) {
                self.clear_horizontal(row);
            }
            let y = row % layout.dimensions[1];
            let z = row / layout.dimensions[1];
            for (yn, zn) in [
                (y.wrapping_sub(1), z),
                (y + 1, z),
                (y, z.wrapping_sub(1)),
                (y, z + 1),
            ] {
                if yn < 1
                    || yn >= layout.dimensions[1] - 1
                    || zn < 1
                    || zn >= layout.dimensions[2] - 1
                {
                    continue;
                }
                let neighbor = yn + layout.dimensions[1] * zn;
                if !self.has_surface(neighbor) {
                    continue;
                }
                let mut changed = false;
                for x in 0..layout.row_words {
                    let from = row * layout.row_words + x;
                    let to = neighbor * layout.row_words + x;
                    let remaining = self.interior[to] & self.occupied(from);
                    changed |= remaining != self.interior[to];
                    self.interior[to] = remaining;
                }
                if changed && !queued[neighbor] {
                    queued[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
        Ok(())
    }

    fn clear_horizontal(&mut self, row: usize) {
        let start = row * self.layout.row_words;
        let end = start + self.layout.row_words;
        let mut carry = true;
        for word in start..end {
            let open = !self.surface[word];
            let mut outside = !self.occupied(word) | u64::from(carry) & open;
            let mut path = open;
            for shift in [1, 2, 4, 8, 16, 32] {
                outside |= (outside << shift) & path;
                path &= path << shift;
            }
            self.interior[word] &= !outside;
            carry = outside >> 63 != 0;
        }
        carry = true;
        for word in (start..end).rev() {
            let open = !self.surface[word];
            let mut outside = !self.occupied(word) | (u64::from(carry) << 63) & open;
            let mut path = open;
            for shift in [1, 2, 4, 8, 16, 32] {
                outside |= (outside >> shift) & path;
                path &= path >> shift;
            }
            self.interior[word] &= !outside;
            carry = outside & 1 != 0;
        }
    }
}

pub(super) struct SetBits(pub u64);
impl Iterator for SetBits {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        if self.0 == 0 {
            return None;
        }
        let bit = self.0.trailing_zeros() as usize;
        self.0 &= self.0 - 1;
        Some(bit)
    }
}
