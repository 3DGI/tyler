pub mod gltf_writer;
pub mod gpkg_writer;
pub mod obj_writer;
mod remeshing;
pub mod tabular;
#[path = "triangle-mesh.rs"]
mod triangle_mesh;
pub mod tsv_writer;
pub mod voxel_remesh;

pub use gpkg_writer::{
    convert_to_gpkg, write_metadata_gpkg, write_tiled_metadata_gpkg, GpkgExportOptions,
};
pub use tabular::{
    semantic_primitive_geometry, tabulate_addresses, tabulate_cityobject_hierarchy,
    tabulate_cityobjects, tabulate_model_metadata, tabulate_semantic_assignments,
    tabulate_semantic_hierarchy, tabulate_semantic_primitives, tabulate_semantics,
    tabulate_tiled_metadata, AddressRow, AddressRowRef, AddressTable, CityObjectHierarchyTable,
    CityObjectRow, CityObjectTable, ColumnOrigin, ColumnSchema, HierarchyRow, IdList, LogicalType,
    MetadataRow, MetadataRowRef, MetadataTable, PrimitiveType, SemanticAssignmentRow,
    SemanticAssignmentTable, SemanticHierarchyRow, SemanticHierarchyTable,
    SemanticPrimitiveGeometry, SemanticPrimitiveRow, SemanticPrimitiveRowRef,
    SemanticPrimitiveTable, SemanticRow, SemanticRowRef, SemanticTable, TableSchema, TileMetadata,
    TiledMetadataTable, Value,
};
pub use tsv_writer::{
    convert_to_tsv, write_addresses_tsv, write_cityobject_hierarchy_tsv, write_cityobjects_tsv,
    write_metadata_tsv, write_semantic_hierarchy_tsv, write_semantics_tsv,
    write_tiled_metadata_tsv, TsvExportOptions, TsvWriteOptions,
};

use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::Result;
use cityjson_lib::CityModel;

/// Optional reconstruction of all input geometry in a single voxel grid.
/// Coordinates and the error budget must be in Cartesian metres. Library callers
/// select source `LoDs` before conversion; remeshing never selects them implicitly.
#[derive(Clone, Copy, Debug)]
pub struct VoxelRemeshOptions {
    /// Requested quality budget, split equally between remeshing and simplification.
    /// This is an estimate, not a certified surface-distance bound.
    pub geometric_error: f32,
}

/// Statistics for the experimental remesh/simplify/attribute-transfer pipeline.
#[derive(Clone, Debug)]
pub struct RemeshReport {
    pub requested_geometric_error: f64,
    /// Padded rectangular grid dimensions, one per Cartesian axis.
    pub grid_dimensions: [u32; 3],
    pub effective_voxel_size: f64,
    pub remeshing_error_estimate: f64,
    pub simplification_error: f64,
    pub encoding_error_estimate: f64,
    /// Selection estimate, including encoding, rather than a guaranteed bound.
    pub estimated_geometric_error: f64,
    pub input_triangles: usize,
    pub remeshed_triangles: usize,
    pub output_triangles: usize,
}

/// GLB export statistics for CLI reporting and future parent-tile generation.
#[derive(Clone, Debug, Default)]
pub struct GlbExportReport {
    pub remeshing: Option<RemeshReport>,
    pub output_bytes: u64,
}

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
    pub voxel_remesh: Option<VoxelRemeshOptions>,
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
            voxel_remesh: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct ObjExportOptions {
    pub clip_bbox: Option<[f64; 6]>,
}

/// Converts a `CityJSON` model to a GLB file.
///
/// # Errors
///
/// Returns an error for invalid remeshing budgets, remeshing
/// with clipping, empty remeshed output, directory creation failures, or GLB writing failures.
pub fn convert_to_glb<P: AsRef<Path>>(
    model: &CityModel,
    output: P,
    options: &ExportOptions,
) -> Result<GlbExportReport> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    gltf_writer::write_city_model_glb(model, output, options)
}

/// Converts a `CityJSON` model to a Wavefront OBJ file.
///
/// # Errors
///
/// Returns an error when the output directory cannot be created or OBJ writing
/// fails.
pub fn convert_to_obj<P: AsRef<Path>>(
    model: &CityModel,
    output: P,
    options: &ObjExportOptions,
) -> Result<()> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    obj_writer::write_city_model_obj(model, output, options)
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
