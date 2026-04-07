//! Transform alignment — reproject and re-quantize coordinates from a source CRS/transform
//! to a reference CRS/transform.
//!
//! Used when merging features from different sources (e.g. buildings from CityJSONL and
//! trees from GeoParquet) into a single tiling dataset.

use anyhow::{Context, Result};
use log::debug;

use crate::parser::Transform;
use crate::proj::Proj;

/// Reproject a single 3D point from `src_epsg` to `dst_epsg` using PROJ.
/// Returns the reprojected `(x, y, z)`.
fn reproject_point(proj: &Proj, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
    proj.convert((x, y, z))
        .map_err(|e| anyhow::anyhow!("PROJ reprojection failed: {e}"))
}

/// Quantize a real-world coordinate to an integer using the given transform.
///
/// `quantized = round((real_world - translate) / scale)`
#[inline]
fn quantize(value: f64, translate: f64, scale: f64) -> i64 {
    ((value - translate) / scale).round() as i64
}

/// Quantize a 3D real-world vertex into quantized integers for CityJSON.
pub fn quantize_vertex(vertex: &[f64; 3], transform: &Transform) -> [i64; 3] {
    [
        quantize(vertex[0], transform.translate[0], transform.scale[0]),
        quantize(vertex[1], transform.translate[1], transform.scale[1]),
        quantize(vertex[2], transform.translate[2], transform.scale[2]),
    ]
}

/// Pipeline for transforming real-world 3D vertices from a source CRS into
/// quantized CityJSON integers in a reference CRS/transform.
pub struct TransformAligner {
    /// PROJ transformer, `None` when source and destination EPSG are the same.
    proj: Option<Proj>,
    /// The reference (destination) transform used for quantization.
    ref_transform: Transform,
}

impl TransformAligner {
    /// Create a new aligner.
    ///
    /// * `src_epsg` — EPSG code of the source data.
    /// * `ref_epsg` — EPSG code of the reference (destination) data.
    /// * `ref_transform` — The reference transform (scale + translate) to quantize into.
    pub fn new(src_epsg: u16, ref_epsg: u16, ref_transform: Transform) -> Result<Self> {
        let proj = if src_epsg != ref_epsg {
            debug!(
                "CRS differs: source EPSG:{} → reference EPSG:{}; will reproject",
                src_epsg, ref_epsg
            );
            let from = format!("EPSG:{src_epsg}");
            let to = format!("EPSG:{ref_epsg}");
            Some(
                Proj::new_known_crs(&from, &to, None)
                    .context("creating PROJ transformer for CRS alignment")?,
            )
        } else {
            debug!("Source and reference CRS are the same (EPSG:{}), no reprojection needed", src_epsg);
            None
        };

        Ok(Self {
            proj,
            ref_transform,
        })
    }

    /// Transform a real-world 3D vertex: optionally reproject and then quantize.
    pub fn align_and_quantize(&self, vertex: &[f64; 3]) -> Result<[i64; 3]> {
        let aligned = if let Some(ref proj) = self.proj {
            let (x, y, z) = reproject_point(proj, vertex[0], vertex[1], vertex[2])?;
            [x, y, z]
        } else {
            *vertex
        };
        Ok(quantize_vertex(&aligned, &self.ref_transform))
    }

    /// Get a reference to the destination transform.
    #[allow(dead_code)]
    pub fn ref_transform(&self) -> &Transform {
        &self.ref_transform
    }
}
