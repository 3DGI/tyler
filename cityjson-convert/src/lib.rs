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
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::Result;
use cityjson_lib::CityModel;

pub type JsonExportOptions = cityjson_lib::json::WriteOptions;

#[derive(Clone, Copy, Debug)]
pub struct CityJsonSeqExportOptions {
    pub scale: [f64; 3],
}

impl Default for CityJsonSeqExportOptions {
    fn default() -> Self {
        Self {
            scale: [0.001, 0.001, 0.001],
        }
    }
}

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
    gltf_writer::write_city_model_glb(model, output, options)
}

/// Converts a `CityJSON` model to a `CityJSON` file.
///
/// # Errors
///
/// Returns an error when the output directory cannot be created or `CityJSON`
/// writing fails.
pub fn convert_to_cityjson<P: AsRef<Path>>(
    model: &CityModel,
    output: P,
    options: &JsonExportOptions,
) -> Result<()> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file = File::create(output)?;
    let mut buffer = Vec::new();
    cityjson_lib::json::to_writer_with_options(&mut buffer, model, *options)?;
    let mut root: serde_json::Value = serde_json::from_slice(&buffer)?;
    if root.get("type").and_then(serde_json::Value::as_str) == Some("CityJSONFeature") {
        if let Some(root) = root.as_object_mut() {
            root.insert(
                "type".to_string(),
                serde_json::Value::String("CityJSON".to_string()),
            );
            root.insert(
                "version".to_string(),
                serde_json::Value::String("2.0".to_string()),
            );
            root.remove("id");
        }
    }
    if options.pretty {
        serde_json::to_writer_pretty(&mut file, &root)?;
    } else {
        serde_json::to_writer(&mut file, &root)?;
    }
    file.write_all(b"\n")?;
    Ok(())
}

/// Converts feature models to a `CityJSONSeq` file.
///
/// # Errors
///
/// Returns an error when the output directory cannot be created or `CityJSONSeq`
/// writing fails.
pub fn convert_to_cityjsonseq<P: AsRef<Path>>(
    base_root: &CityModel,
    feature_models: &[CityModel],
    output: P,
    options: &CityJsonSeqExportOptions,
) -> Result<()> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file = File::create(output)?;
    write_cityjsonseq(&mut file, base_root, feature_models, options)
}

/// Writes feature models to a `CityJSONSeq` writer.
///
/// # Errors
///
/// Returns an error when `CityJSONSeq` writing fails.
pub fn write_cityjsonseq<W: Write>(
    writer: &mut W,
    base_root: &CityModel,
    feature_models: &[CityModel],
    options: &CityJsonSeqExportOptions,
) -> Result<()> {
    cityjson_lib::json::write_cityjsonseq_auto_transform_refs(
        writer,
        base_root,
        feature_models.iter(),
        options.scale,
    )?;
    Ok(())
}
