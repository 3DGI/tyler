#![allow(
    clippy::manual_is_multiple_of,
    clippy::needless_as_bytes,
    clippy::needless_borrows_for_generic_args,
    clippy::to_string_trait_impl,
    clippy::too_many_arguments,
    clippy::unnecessary_fallible_conversions,
    clippy::unnecessary_to_owned,
    clippy::unnecessary_unwrap,
    clippy::useless_vec
)]

pub mod gltf_writer;
mod proj;

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use anyhow::Result;
use cityjson_lib::CityModel;
use log::info;

#[derive(Clone, Debug, Default)]
pub enum GeometryPlacement {
    #[default]
    SourceCoordinates,
    EcefRelative {
        source_crs: String,
        origin: [f64; 3],
    },
    Enu {
        source_crs: String,
        ecef_origin: [f64; 3],
        east: [f64; 3],
        north: [f64; 3],
        up: [f64; 3],
    },
}

#[derive(Clone, Debug)]
pub struct GeographicClipRegion {
    pub source_crs: String,
    pub west: f64,
    pub south: f64,
    pub east: f64,
    pub north: f64,
}

#[derive(Clone, Debug)]
#[allow(clippy::struct_excessive_bools)]
pub struct ExportOptions {
    pub native_glb_color: String,
    pub metadata_class_name: String,
    pub feature_type_colors: BTreeMap<String, String>,
    pub geometry_placement: GeometryPlacement,
    pub clip_bbox: Option<[f64; 6]>,
    pub clip_geographic_region: Option<GeographicClipRegion>,
    pub smooth_normals: bool,
    pub quantize_geometry: bool,
    pub meshopt_compression: bool,
}

impl Default for ExportOptions {
    fn default() -> Self {
        Self {
            native_glb_color: "#FFC0CB".to_string(),
            metadata_class_name: "cityobject".to_string(),
            feature_type_colors: BTreeMap::new(),
            geometry_placement: GeometryPlacement::default(),
            clip_bbox: None,
            clip_geographic_region: None,
            smooth_normals: false,
            quantize_geometry: true,
            meshopt_compression: true,
        }
    }
}

/// Converts a `CityJSON` model to a GLB file.
///
/// # Errors
///
/// Returns an error when the output directory cannot be created or GLB writing
/// fails.
pub fn convert_to_glb<P: AsRef<Path>>(
    model: &CityModel,
    output: P,
    options: &ExportOptions,
) -> Result<()> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    info!("Writing GLB output to {}", output.display());
    gltf_writer::write_city_model_glb(model, output, options)
}
