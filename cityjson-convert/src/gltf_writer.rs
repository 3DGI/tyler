#![allow(clippy::too_many_lines)]

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::convert::TryFrom;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::proj::Proj;
use crate::{ExportOptions, GeographicClipRegion, GeometryPlacement};
use anyhow::{bail, Context, Result};
use cityjson_lib::cityjson::v2_0::{AttributeValue, CityObject, GeometryType, VertexIndex};
use cityjson_lib::CityModel;
use earcutr::earcut;
use gltf::json;
use log::debug;
use meshopt::{
    encode_index_buffer, encode_vertex_buffer, generate_vertex_remap,
    optimize_overdraw_in_place_decoder, optimize_vertex_cache, optimize_vertex_fetch,
    quantize_snorm, remap_index_buffer, remap_vertex_buffer, DecodePosition,
};
use serde_json::{json as json_value, Map as JsonMap, Value as JsonValue};

const GLTF_VERSION: &str = "2.0";
const OVERDRAW_THRESHOLD: f32 = 1.05;
const QUANTIZATION_EXTENSION: &str = "KHR_mesh_quantization";
const MESHOPT_EXTENSION: &str = "EXT_meshopt_compression";
const MESH_FEATURES_EXTENSION: &str = "EXT_mesh_features";
const STRUCTURAL_METADATA_EXTENSION: &str = "EXT_structural_metadata";
const QUANTIZED_POSITION_STRIDE: usize = std::mem::size_of::<QuantizedPosition>();
const QUANTIZED_NORMAL_STRIDE: usize = std::mem::size_of::<QuantizedNormal>();
const CLIP_PLANE_EPSILON: f64 = 1.0e-12;
const GEOGRAPHIC_CLIP_INTERSECTION_ITERATIONS: usize = 64;

thread_local! {
    static PROJ_TRANSFORM_CACHE: RefCell<HashMap<(String, String), Proj>> = RefCell::new(HashMap::new());
}

/// Parse hex color string (#RRGGBB) to RGBA f32 array [R, G, B, A]
fn hex_to_rgba(hex: &str) -> Result<[f32; 4], anyhow::Error> {
    if hex.len() != 7 || !hex.starts_with('#') {
        return Err(anyhow::anyhow!(
            "Invalid hex color format: expected #RRGGBB"
        ));
    }
    let hex_digits = &hex[1..];
    let r = u8::from_str_radix(&hex_digits[0..2], 16)?;
    let g = u8::from_str_radix(&hex_digits[2..4], 16)?;
    let b = u8::from_str_radix(&hex_digits[4..6], 16)?;
    Ok([
        f32::from(r) / 255.0,
        f32::from(g) / 255.0,
        f32::from(b) / 255.0,
        1.0,
    ])
}

/// Create default PBR material matching pg2b3dm's structure
fn create_default_material(base_color: &str) -> Result<json::Material, anyhow::Error> {
    let base_color_rgba = hex_to_rgba(base_color)?;

    let roughness_factor = 128.0 / 255.0;
    let metallic_factor = 0.0;

    Ok(json::Material {
        name: None,
        extensions: None,
        extras: Option::default(),
        pbr_metallic_roughness: json::material::PbrMetallicRoughness {
            base_color_factor: json::material::PbrBaseColorFactor(base_color_rgba),
            metallic_factor: json::material::StrengthFactor(metallic_factor),
            roughness_factor: json::material::StrengthFactor(roughness_factor),
            base_color_texture: None,
            metallic_roughness_texture: None,
            extensions: None,
            extras: Option::default(),
        },
        normal_texture: None,
        occlusion_texture: None,
        emissive_texture: None,
        emissive_factor: json::material::EmissiveFactor([0.0, 0.0, 0.0]),
        alpha_mode: json::validation::Checked::Valid(json::material::AlphaMode::Opaque),
        alpha_cutoff: None,
        double_sided: true,
    })
}

/// Writes a `CityJSON` model as a binary glTF file.
///
/// # Errors
///
/// Returns an error when the model geometry cannot be read or triangulated, or
/// when the output GLB cannot be created.
pub fn write_city_model_glb<P: AsRef<Path>>(
    model: &CityModel,
    output_path: P,
    options: &ExportOptions,
) -> Result<()> {
    let coordinate_transform = CoordinateTransform::from_model(model, options)?;
    let clip_volume = ClipVolume::from_options(options)?;
    let mut collector =
        MeshCollector::new(coordinate_transform, options.smooth_normals, clip_volume);
    collector.add_model(model)?;
    let processed = collector.finish()?;
    debug!(
        "Processed {} vertices and {} indices for the output GLB",
        processed.vertex_count(),
        processed.index_count()
    );
    processed.write_glb(output_path, options)
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct Vertex {
    position: [f32; 3],
    normal: [f32; 3],
    feature_id: u32,
}

impl DecodePosition for Vertex {
    fn decode_position(&self) -> [f32; 3] {
        self.position
    }
}

#[derive(Clone, Copy, Debug)]
struct Bounds {
    min: [f32; 3],
    max: [f32; 3],
}

impl Bounds {
    fn empty() -> Self {
        Self {
            min: [f32::INFINITY; 3],
            max: [f32::NEG_INFINITY; 3],
        }
    }

    fn add_point(&mut self, point: [f32; 3]) {
        for (axis, coordinate) in point.iter().enumerate() {
            self.min[axis] = self.min[axis].min(*coordinate);
            self.max[axis] = self.max[axis].max(*coordinate);
        }
    }

    fn from_vertices(vertices: &[Vertex]) -> Option<Self> {
        let mut bounds = Self::empty();
        let mut has_vertices = false;
        for vertex in vertices {
            bounds.add_point(vertex.position);
            has_vertices = true;
        }
        has_vertices.then_some(bounds)
    }

    fn center(&self) -> [f32; 3] {
        [
            f32::midpoint(self.min[0], self.max[0]),
            f32::midpoint(self.min[1], self.max[1]),
            f32::midpoint(self.min[2], self.max[2]),
        ]
    }
}

enum IndexBuffer {
    U16(Vec<u16>),
    U32(Vec<u32>),
}

impl IndexBuffer {
    fn byte_stride(&self) -> usize {
        match self {
            Self::U16(_) => std::mem::size_of::<u16>(),
            Self::U32(_) => std::mem::size_of::<u32>(),
        }
    }

    fn component_type(&self) -> json::validation::Checked<json::accessor::GenericComponentType> {
        match self {
            Self::U16(_) => json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U16,
            )),
            Self::U32(_) => json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U32,
            )),
        }
    }

    fn byte_length(&self) -> usize {
        match self {
            Self::U16(indices) => indices.len() * std::mem::size_of::<u16>(),
            Self::U32(indices) => indices.len() * std::mem::size_of::<u32>(),
        }
    }

    fn count(&self) -> usize {
        match self {
            Self::U16(indices) => indices.len(),
            Self::U32(indices) => indices.len(),
        }
    }

    fn max_value(&self) -> u32 {
        match self {
            Self::U16(indices) => indices.iter().copied().max().map_or(0, u32::from),
            Self::U32(indices) => indices.iter().copied().max().unwrap_or(0),
        }
    }

    fn as_u32_vec(&self) -> Vec<u32> {
        match self {
            Self::U16(indices) => indices.iter().copied().map(u32::from).collect(),
            Self::U32(indices) => indices.clone(),
        }
    }
}

enum FeatureIdBuffer {
    U8(Vec<u8>),
    U16(Vec<u16>),
}

impl FeatureIdBuffer {
    fn from_feature_ids(feature_ids: &[u32]) -> Result<Self> {
        let max_value = feature_ids.iter().copied().max().unwrap_or(0);
        if u8::try_from(max_value).is_ok() {
            Ok(Self::U8(
                feature_ids
                    .iter()
                    .copied()
                    .map(|feature_id| {
                        u8::try_from(feature_id)
                            .expect("feature ID should fit in u8 after max check")
                    })
                    .collect(),
            ))
        } else if u16::try_from(max_value).is_ok() {
            Ok(Self::U16(
                feature_ids
                    .iter()
                    .copied()
                    .map(|feature_id| {
                        u16::try_from(feature_id)
                            .expect("feature ID should fit in u16 after max check")
                    })
                    .collect(),
            ))
        } else {
            bail!(
                "feature ID attribute stream exceeds glTF limits: max feature ID {max_value} does not fit in u16"
            );
        }
    }

    fn byte_stride(&self) -> usize {
        match self {
            Self::U8(_) => std::mem::size_of::<u8>(),
            Self::U16(_) => std::mem::size_of::<u16>(),
        }
    }

    fn component_type(&self) -> json::validation::Checked<json::accessor::GenericComponentType> {
        match self {
            Self::U8(_) => json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U8,
            )),
            Self::U16(_) => json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U16,
            )),
        }
    }

    fn max_value(&self) -> u32 {
        match self {
            Self::U8(feature_ids) => feature_ids.iter().copied().max().map_or(0, u32::from),
            Self::U16(feature_ids) => feature_ids.iter().copied().max().map_or(0, u32::from),
        }
    }

    fn count(&self) -> usize {
        match self {
            Self::U8(feature_ids) => feature_ids.len(),
            Self::U16(feature_ids) => feature_ids.len(),
        }
    }

    fn meshopt_byte_stride() -> usize {
        4
    }
}

#[derive(Clone, Debug)]
enum MetadataValue {
    Bool(bool),
    Int(i32),
    Float(f32),
    String(String),
}

#[derive(Clone, Debug)]
struct FeatureRecord {
    object_id: String,
    feature_type: String,
    attributes: BTreeMap<String, MetadataValue>,
}

#[derive(Clone, Copy, Debug)]
struct SourceTriangle {
    feature_id: u32,
    source_positions: [[f64; 3]; 3],
    local_positions: [[f32; 3]; 3],
    face_normal: [f32; 3],
}

#[derive(Clone, Copy, Debug)]
struct ClipVertex {
    source_position: [f64; 3],
    clip_position: [f64; 3],
    barycentric: [f32; 3],
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

#[derive(Default)]
struct RawPrimitiveMesh {
    triangles: Vec<SourceTriangle>,
    source_vertex_normals: HashMap<SmoothVertexKey, [f32; 3]>,
}

#[derive(Default)]
struct BuiltPrimitiveMesh {
    vertices: Vec<Vertex>,
    indices: Vec<u32>,
}

struct ProcessedPrimitiveMesh {
    feature_type: String,
    vertices: Vec<Vertex>,
    indices: Vec<u32>,
}

struct ProcessedScene {
    primitives: Vec<ProcessedPrimitiveMesh>,
    features: Vec<FeatureRecord>,
    center: [f32; 3],
    node_translation_base: [f32; 3],
    bounds: Option<Bounds>,
}

struct MeshCollector {
    features: Vec<FeatureRecord>,
    primitives: BTreeMap<String, RawPrimitiveMesh>,
    coordinate_transform: CoordinateTransform,
    smooth_normals: bool,
    clip_volume: Option<ClipVolume>,
}

struct CoordinateTransform {
    placement: CoordinatePlacement,
    node_translation_base: [f32; 3],
}

enum CoordinatePlacement {
    SourceCoordinates,
    EcefRelative {
        vertex_transformer: Option<CachedProjTransform>,
        origin: [f64; 3],
    },
    Enu {
        vertex_transformer: Option<CachedProjTransform>,
        ecef_origin: [f64; 3],
        east: [f64; 3],
        north: [f64; 3],
        up: [f64; 3],
    },
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct SmoothVertexKey {
    position_bits: [u32; 3],
    feature_id: u32,
}

#[derive(Clone, Debug)]
struct CachedProjTransform {
    key: (String, String),
}

impl CachedProjTransform {
    fn new(source_crs: String, target_crs: &'static str) -> Self {
        Self {
            key: (source_crs, target_crs.to_owned()),
        }
    }

    fn convert(&self, point: [f64; 3]) -> Result<[f64; 3]> {
        PROJ_TRANSFORM_CACHE.with(|cache| {
            let mut cache = cache.borrow_mut();
            if !cache.contains_key(&self.key) {
                let transformer = Proj::new_known_crs(&self.key.0, &self.key.1, None)
                    .with_context(|| {
                        format!(
                            "failed to create {} to {} transform",
                            self.key.0, self.key.1
                        )
                    })?;
                cache.insert(self.key.clone(), transformer);
            }
            let transformer = cache
                .get(&self.key)
                .expect("cached PROJ transform should have been inserted");
            let transformed = transformer.convert((point[0], point[1], point[2]))?;
            Ok([transformed.0, transformed.1, transformed.2])
        })
    }
}

#[derive(Clone, Copy, Debug)]
struct VertexAccessors {
    positions: json::Index<json::Accessor>,
    normals: json::Index<json::Accessor>,
    feature_ids: Option<json::Index<json::Accessor>>,
}

#[derive(Clone, Debug)]
struct PrimitiveEncoding {
    feature_type: String,
    primitive: json::mesh::Primitive,
    feature_count: usize,
}

#[derive(Clone, Debug)]
struct StructuralMetadataColumn {
    property: JsonValue,
    property_table_entry: JsonValue,
}

#[derive(Clone, Debug)]
struct StructuralMetadataExtension {
    class_name: String,
    columns: BTreeMap<String, StructuralMetadataColumn>,
    feature_count: usize,
}

#[derive(Clone, Debug)]
struct MeshoptBufferView {
    buffer: u32,
    byte_offset: u64,
    byte_length: u64,
    byte_stride: u64,
    count: u64,
    mode: &'static str,
    filter: Option<&'static str>,
}

#[derive(Default)]
struct BufferBuilder {
    bytes: Vec<u8>,
    buffer_views: Vec<json::buffer::View>,
    accessors: Vec<json::Accessor>,
    meshopt_views: Vec<Option<MeshoptBufferView>>,
    fallback_buffer_length: usize,
}

struct EncodedGlb {
    root: json::Root,
    bin_buffer: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct QuantizedPosition {
    position: [i16; 4],
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct QuantizedNormal {
    normal: [i8; 4],
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct PaddedFeatureIdU8 {
    feature_id: u8,
    _padding: [u8; 3],
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct PaddedFeatureIdU16 {
    feature_id: u16,
    _padding: u16,
}

fn quantize_position_component(value: f32, scale: f32) -> i16 {
    i16::try_from(quantize_snorm(value / scale, 16))
        .expect("quantized position component should fit in i16")
}

fn quantize_normal_component(value: f32) -> i8 {
    i8::try_from(quantize_snorm(value, 8)).expect("quantized normal component should fit in i8")
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct FloatPosition {
    position: [f32; 3],
}

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct FloatNormal {
    normal: [f32; 3],
}

struct QuantizedPrimitiveMesh {
    feature_type: String,
    positions: Vec<QuantizedPosition>,
    normals: Vec<QuantizedNormal>,
    feature_ids: Vec<u32>,
    indices: IndexBuffer,
    position_bounds: QuantizedBounds,
}

#[derive(Clone, Copy, Debug)]
struct QuantizedBounds {
    min: [i16; 3],
    max: [i16; 3],
}

struct QuantizedScene {
    primitives: Vec<QuantizedPrimitiveMesh>,
    position_scale: f32,
    center: [f32; 3],
    features: Vec<FeatureRecord>,
}

impl SmoothVertexKey {
    fn new(position: [f32; 3], feature_id: u32) -> Self {
        Self {
            position_bits: position.map(f32::to_bits),
            feature_id,
        }
    }
}

impl CoordinateTransform {
    fn from_model(_model: &CityModel, options: &ExportOptions) -> Result<Self> {
        let placement = match &options.geometry_placement {
            GeometryPlacement::SourceCoordinates => CoordinatePlacement::SourceCoordinates,
            GeometryPlacement::EcefRelative { source_crs, origin } => {
                let source_crs = canonical_epsg_crs(source_crs)?;
                let vertex_transformer = source_to_ecef_transformer(&source_crs);
                CoordinatePlacement::EcefRelative {
                    vertex_transformer,
                    origin: *origin,
                }
            }
            GeometryPlacement::Enu {
                source_crs,
                ecef_origin,
                east,
                north,
                up,
            } => {
                let source_crs = canonical_epsg_crs(source_crs)?;
                let vertex_transformer = source_to_ecef_transformer(&source_crs);
                CoordinatePlacement::Enu {
                    vertex_transformer,
                    ecef_origin: *ecef_origin,
                    east: *east,
                    north: *north,
                    up: *up,
                }
            }
        };

        Ok(Self {
            placement,
            node_translation_base: [0.0; 3],
        })
    }

    fn transform_position(&self, position: [f64; 3]) -> Result<[f32; 3]> {
        match &self.placement {
            CoordinatePlacement::SourceCoordinates => Ok([
                f64_to_f32_checked(position[0], "x", None)?,
                f64_to_f32_checked(position[1], "y", None)?,
                f64_to_f32_checked(position[2], "z", None)?,
            ]),
            CoordinatePlacement::EcefRelative {
                vertex_transformer,
                origin,
            } => {
                let ecef = transform_to_ecef(position, vertex_transformer.as_ref())?;
                Ok([
                    f64_to_f32_checked(ecef[0] - origin[0], "ecef relative x", None)?,
                    f64_to_f32_checked(ecef[1] - origin[1], "ecef relative y", None)?,
                    f64_to_f32_checked(ecef[2] - origin[2], "ecef relative z", None)?,
                ])
            }
            CoordinatePlacement::Enu {
                vertex_transformer,
                ecef_origin,
                east,
                north,
                up,
            } => {
                let ecef = transform_to_ecef(position, vertex_transformer.as_ref())?;
                let delta = [
                    ecef[0] - ecef_origin[0],
                    ecef[1] - ecef_origin[1],
                    ecef[2] - ecef_origin[2],
                ];
                Ok([
                    f64_to_f32_checked(dot(delta, *east), "enu east", None)?,
                    f64_to_f32_checked(dot(delta, *north), "enu north", None)?,
                    f64_to_f32_checked(dot(delta, *up), "enu up", None)?,
                ])
            }
        }
    }
}

impl ClipVolume {
    fn from_options(options: &ExportOptions) -> Result<Option<Self>> {
        if let Some(region) = &options.clip_geographic_region {
            return Self::geographic_region(region).map(Some);
        }

        Ok(options.clip_bbox.map(Self::SourceBbox))
    }

    fn geographic_region(region: &GeographicClipRegion) -> Result<Self> {
        let source_crs = canonical_epsg_crs(&region.source_crs)?;
        let transformer = (source_crs != "EPSG:4979")
            .then(|| CachedProjTransform::new(source_crs.clone(), "EPSG:4979"));
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
                let [x, y, z] = source_position;
                let geographic = if let Some(transformer) = &region.transformer {
                    transformer
                        .convert([x, y, z])
                        .context("failed to project position to EPSG:4979 for clipping")?
                } else {
                    [x, y, z]
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

    #[allow(clippy::cast_possible_truncation)]
    fn interpolate_clip_vertex(
        &self,
        start: ClipVertex,
        end: ClipVertex,
        t: f64,
    ) -> Result<ClipVertex> {
        let tf32 = t as f32;
        let source_position = [
            start.source_position[0] + (end.source_position[0] - start.source_position[0]) * t,
            start.source_position[1] + (end.source_position[1] - start.source_position[1]) * t,
            start.source_position[2] + (end.source_position[2] - start.source_position[2]) * t,
        ];
        Ok(ClipVertex {
            source_position,
            clip_position: self.clip_position(source_position)?,
            barycentric: [
                start.barycentric[0] + (end.barycentric[0] - start.barycentric[0]) * tf32,
                start.barycentric[1] + (end.barycentric[1] - start.barycentric[1]) * tf32,
                start.barycentric[2] + (end.barycentric[2] - start.barycentric[2]) * tf32,
            ],
        })
    }
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

fn source_to_ecef_transformer(source_crs: &str) -> Option<CachedProjTransform> {
    (source_crs != "EPSG:4978")
        .then(|| CachedProjTransform::new(source_crs.to_owned(), "EPSG:4978"))
}

fn transform_to_ecef(
    point: [f64; 3],
    vertex_transformer: Option<&CachedProjTransform>,
) -> Result<[f64; 3]> {
    let [x, y, z] = point;
    if let Some(transformer) = vertex_transformer {
        let transformed = transformer
            .convert([x, y, z])
            .context("failed to project position to ECEF")?;
        Ok(transformed)
    } else {
        Ok(point)
    }
}

fn dot(lhs: [f64; 3], rhs: [f64; 3]) -> f64 {
    lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
}

#[allow(clippy::cast_possible_truncation)]
fn f64_to_f32_checked(value: f64, axis: &str, vertex_id: Option<u32>) -> Result<f32> {
    if !value.is_finite() || value < f64::from(f32::MIN) || value > f64::from(f32::MAX) {
        if let Some(vertex_id) = vertex_id {
            bail!("vertex {vertex_id} {axis} coordinate is outside the f32 range");
        }
        bail!("{axis} coordinate is outside the f32 range");
    }
    {
        Ok(value as f32)
    }
}

impl MeshCollector {
    fn new(
        coordinate_transform: CoordinateTransform,
        smooth_normals: bool,
        clip_volume: Option<ClipVolume>,
    ) -> Self {
        Self {
            features: Vec::new(),
            primitives: BTreeMap::new(),
            coordinate_transform,
            smooth_normals,
            clip_volume,
        }
    }

    fn add_model(&mut self, model: &CityModel) -> Result<()> {
        for (object_id, cityobject) in model.cityobjects().iter() {
            let Some(geometry_handles) = cityobject.geometry() else {
                continue;
            };
            let feature_type = cityobject.type_cityobject().to_string();
            let mut feature_index = None;
            for geometry_handle in geometry_handles {
                let geometry = model.resolve_geometry(*geometry_handle)?;
                let feature_index = *feature_index.get_or_insert_with(|| {
                    let feature_index =
                        u32::try_from(self.features.len()).expect("feature count within u32 range");
                    self.features.push(FeatureRecord {
                        object_id: object_id.to_string(),
                        feature_type: feature_type.clone(),
                        attributes: Self::collect_feature_attributes(model, cityobject),
                    });
                    feature_index
                });
                match geometry.geometry().type_geometry() {
                    GeometryType::MultiSurface | GeometryType::CompositeSurface => {
                        let Some(boundary) = geometry.geometry().boundaries() else {
                            continue;
                        };
                        for surface in boundary.to_nested_multi_or_composite_surface()? {
                            self.add_surface(&feature_type, feature_index, &surface, model)?;
                        }
                    }
                    GeometryType::Solid => {
                        let Some(boundary) = geometry.geometry().boundaries() else {
                            continue;
                        };
                        for shell in boundary.to_nested_solid()? {
                            for surface in shell {
                                self.add_surface(&feature_type, feature_index, &surface, model)?;
                            }
                        }
                    }
                    GeometryType::MultiSolid | GeometryType::CompositeSolid => {
                        let Some(boundary) = geometry.geometry().boundaries() else {
                            continue;
                        };
                        for solid in boundary.to_nested_multi_or_composite_solid()? {
                            for shell in solid {
                                for surface in shell {
                                    self.add_surface(
                                        &feature_type,
                                        feature_index,
                                        &surface,
                                        model,
                                    )?;
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
        }

        Ok(())
    }

    fn finish(self) -> Result<ProcessedScene> {
        ProcessedScene::from_collector(self)
    }

    fn add_surface(
        &mut self,
        feature_type: &str,
        feature_id: u32,
        surface: &[Vec<u32>],
        model: &CityModel,
    ) -> Result<()> {
        if surface.is_empty() {
            return Ok(());
        }
        let exterior = &surface[0];
        if exterior.len() < 3 {
            return Ok(());
        }

        let mut source_positions: Vec<[f64; 3]> = Vec::new();
        let mut local_positions: Vec<[f32; 3]> = Vec::new();
        let mut flat_coords: Vec<f64> = Vec::new();
        let mut hole_indices: Vec<usize> = Vec::new();
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
                let source_position = vertex.to_array();
                let local_position = self
                    .coordinate_transform
                    .transform_position(source_position)?;
                source_positions.push(source_position);
                local_positions.push(local_position);
                vertex_count += 1;
            }
        }

        if local_positions.len() < 3 {
            return Ok(());
        }

        let drop_axis = Self::find_projection_axis(&local_positions[..exterior.len()]);
        for pos in &local_positions {
            match drop_axis {
                0 => {
                    flat_coords.push(f64::from(pos[1]));
                    flat_coords.push(f64::from(pos[2]));
                }
                1 => {
                    flat_coords.push(f64::from(pos[0]));
                    flat_coords.push(f64::from(pos[2]));
                }
                _ => {
                    flat_coords.push(f64::from(pos[0]));
                    flat_coords.push(f64::from(pos[1]));
                }
            }
        }

        let triangulated =
            earcut(&flat_coords, &hole_indices, 2).context("Failed to triangulate surface")?;
        if triangulated.len() < 3 {
            return Ok(());
        }

        let primitive = self.primitives.entry(feature_type.to_string()).or_default();
        for tri in triangulated.chunks_exact(3) {
            let source_triangle = [
                source_positions[tri[0]],
                source_positions[tri[1]],
                source_positions[tri[2]],
            ];
            let local_triangle = [
                local_positions[tri[0]],
                local_positions[tri[1]],
                local_positions[tri[2]],
            ];
            primitive.add_triangle(
                feature_id,
                source_triangle,
                local_triangle,
                Self::compute_face_normal,
                self.smooth_normals,
            );
        }

        Ok(())
    }

    fn collect_feature_attributes(
        model: &CityModel,
        cityobject: &CityObject<cityjson_lib::cityjson::resources::storage::OwnedStringStorage>,
    ) -> BTreeMap<String, MetadataValue> {
        let mut attributes = Self::extract_attributes(cityobject);
        if !attributes.is_empty() {
            return attributes;
        }

        let Some(parent_handle) = cityobject
            .parents()
            .and_then(|parents| parents.first())
            .copied()
        else {
            return attributes;
        };
        let Some(parent) = model.cityobjects().get(parent_handle) else {
            return attributes;
        };
        attributes = Self::extract_attributes(parent);
        attributes
    }

    fn extract_attributes(
        cityobject: &CityObject<cityjson_lib::cityjson::resources::storage::OwnedStringStorage>,
    ) -> BTreeMap<String, MetadataValue> {
        let mut attributes = BTreeMap::new();
        let Some(cityjson_attributes) = cityobject.attributes() else {
            return attributes;
        };

        for (name, value) in cityjson_attributes.iter() {
            let Some(value) = Self::convert_attribute_value(value) else {
                continue;
            };
            attributes.insert(name.clone(), value);
        }

        attributes
    }

    fn convert_attribute_value(
        value: &AttributeValue<cityjson_lib::cityjson::resources::storage::OwnedStringStorage>,
    ) -> Option<MetadataValue> {
        match value {
            AttributeValue::Bool(value) => Some(MetadataValue::Bool(*value)),
            AttributeValue::Unsigned(value) => i32::try_from(*value).ok().map(MetadataValue::Int),
            AttributeValue::Integer(value) => i32::try_from(*value).ok().map(MetadataValue::Int),
            AttributeValue::Float(value) => {
                if value.is_finite()
                    && *value >= f64::from(f32::MIN)
                    && *value <= f64::from(f32::MAX)
                {
                    #[allow(clippy::cast_possible_truncation)]
                    Some(MetadataValue::Float(*value as f32))
                } else {
                    None
                }
            }
            AttributeValue::String(value) => Some(MetadataValue::String(value.clone())),
            _ => None,
        }
    }

    fn find_projection_axis(positions: &[[f32; 3]]) -> usize {
        if let Some(normal) = Self::compute_polygon_normal(positions) {
            return Self::dominant_axis(normal);
        }

        let mut min = [f32::INFINITY; 3];
        let mut max = [f32::NEG_INFINITY; 3];
        for pos in positions {
            for (axis, coordinate) in pos.iter().enumerate() {
                min[axis] = min[axis].min(*coordinate);
                max[axis] = max[axis].max(*coordinate);
            }
        }

        (0..3)
            .min_by(|&lhs, &rhs| {
                (max[lhs] - min[lhs])
                    .partial_cmp(&(max[rhs] - min[rhs]))
                    .unwrap()
            })
            .unwrap_or(2)
    }

    fn compute_polygon_normal(positions: &[[f32; 3]]) -> Option<[f32; 3]> {
        if positions.len() < 3 {
            return None;
        }

        let mut normal = [0.0_f32; 3];
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
        if length > f32::EPSILON {
            Some([normal[0] / length, normal[1] / length, normal[2] / length])
        } else {
            None
        }
    }

    fn dominant_axis(normal: [f32; 3]) -> usize {
        (0..3)
            .max_by(|&lhs, &rhs| normal[lhs].abs().partial_cmp(&normal[rhs].abs()).unwrap())
            .unwrap_or(2)
    }

    fn compute_face_normal(p0: [f32; 3], p1: [f32; 3], p2: [f32; 3]) -> [f32; 3] {
        let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let normal = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let length = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if length > f32::EPSILON {
            [normal[0] / length, normal[1] / length, normal[2] / length]
        } else {
            [0.0, 0.0, 1.0]
        }
    }
}

impl RawPrimitiveMesh {
    fn add_triangle<F>(
        &mut self,
        feature_id: u32,
        source_positions: [[f64; 3]; 3],
        local_positions: [[f32; 3]; 3],
        compute_face_normal: F,
        smooth_normals: bool,
    ) where
        F: Fn([f32; 3], [f32; 3], [f32; 3]) -> [f32; 3],
    {
        let face_normal =
            compute_face_normal(local_positions[0], local_positions[1], local_positions[2]);
        self.triangles.push(SourceTriangle {
            feature_id,
            source_positions,
            local_positions,
            face_normal,
        });

        if !smooth_normals {
            return;
        }

        for position in local_positions {
            let key = SmoothVertexKey::new(position, feature_id);
            let entry = self.source_vertex_normals.entry(key).or_insert([0.0; 3]);
            entry[0] += face_normal[0];
            entry[1] += face_normal[1];
            entry[2] += face_normal[2];
        }
    }

    fn build(
        &self,
        coordinate_transform: &CoordinateTransform,
        smooth_normals: bool,
        clip_volume: Option<&ClipVolume>,
    ) -> Result<BuiltPrimitiveMesh> {
        let mut built = BuiltPrimitiveMesh::default();
        let source_vertex_normals = smooth_normals.then(|| self.normalized_source_vertex_normals());
        let mut vertex_map = HashMap::new();

        for triangle in &self.triangles {
            let clipped_triangles = if let Some(clip_volume) = clip_volume {
                Self::clip_triangle_to_volume(triangle.source_positions, clip_volume)?
            } else {
                vec![[
                    ClipVertex {
                        source_position: triangle.source_positions[0],
                        clip_position: triangle.source_positions[0],
                        barycentric: [1.0, 0.0, 0.0],
                    },
                    ClipVertex {
                        source_position: triangle.source_positions[1],
                        clip_position: triangle.source_positions[1],
                        barycentric: [0.0, 1.0, 0.0],
                    },
                    ClipVertex {
                        source_position: triangle.source_positions[2],
                        clip_position: triangle.source_positions[2],
                        barycentric: [0.0, 0.0, 1.0],
                    },
                ]]
            };

            let source_keys = triangle
                .local_positions
                .map(|position| SmoothVertexKey::new(position, triangle.feature_id));

            for clipped_triangle in clipped_triangles {
                let mut positions = [[0.0; 3]; 3];
                let mut normals = [[0.0; 3]; 3];

                for corner in 0..3 {
                    let clipped_corner = clipped_triangle[corner];
                    positions[corner] = if clip_volume.is_some() {
                        coordinate_transform.transform_position(clipped_corner.source_position)?
                    } else {
                        triangle.local_positions[corner]
                    };
                    normals[corner] = if let Some(source_vertex_normals) = &source_vertex_normals {
                        let interpolated = interpolate_normal(
                            clipped_corner.barycentric,
                            source_keys,
                            source_vertex_normals,
                        );
                        normalize_vector(interpolated)
                    } else {
                        triangle.face_normal
                    };
                }

                if is_degenerate_triangle(positions) {
                    continue;
                }

                if smooth_normals {
                    for corner in 0..3 {
                        let vertex_index = built.vertex_index_for_position(
                            positions[corner],
                            normals[corner],
                            triangle.feature_id,
                            &mut vertex_map,
                        );
                        built.indices.push(
                            u32::try_from(vertex_index)
                                .context("GLB vertex count exceeds u32 range")?,
                        );
                    }
                } else {
                    let base_index = u32::try_from(built.vertices.len())
                        .context("GLB vertex count exceeds u32 range")?;
                    for corner in 0..3 {
                        built.vertices.push(Vertex {
                            position: positions[corner],
                            normal: normals[corner],
                            feature_id: triangle.feature_id,
                        });
                    }
                    built
                        .indices
                        .extend_from_slice(&[base_index, base_index + 1, base_index + 2]);
                }
            }
        }

        Ok(built)
    }

    fn normalized_source_vertex_normals(&self) -> HashMap<SmoothVertexKey, [f32; 3]> {
        self.source_vertex_normals
            .iter()
            .map(|(key, normal)| (*key, normalize_vector(*normal)))
            .collect()
    }

    fn clip_triangle_to_volume(
        triangle: [[f64; 3]; 3],
        clip_volume: &ClipVolume,
    ) -> Result<Vec<[ClipVertex; 3]>> {
        let mut polygon = vec![
            ClipVertex {
                source_position: triangle[0],
                clip_position: clip_volume.clip_position(triangle[0])?,
                barycentric: [1.0, 0.0, 0.0],
            },
            ClipVertex {
                source_position: triangle[1],
                clip_position: clip_volume.clip_position(triangle[1])?,
                barycentric: [0.0, 1.0, 0.0],
            },
            ClipVertex {
                source_position: triangle[2],
                clip_position: clip_volume.clip_position(triangle[2])?,
                barycentric: [0.0, 0.0, 1.0],
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
            let clipped_triangle = [polygon[0], polygon[index], polygon[index + 1]];
            if is_degenerate_source_triangle(clipped_triangle.map(|vertex| vertex.source_position))
            {
                continue;
            }
            triangles.push(clipped_triangle);
        }

        Ok(triangles)
    }
}

impl BuiltPrimitiveMesh {
    fn vertex_index_for_position(
        &mut self,
        position: [f32; 3],
        normal: [f32; 3],
        feature_id: u32,
        vertex_map: &mut HashMap<SmoothVertexKey, usize>,
    ) -> usize {
        let key = SmoothVertexKey::new(position, feature_id);
        if let Some(index) = vertex_map.get(&key).copied() {
            return index;
        }

        let index = self.vertices.len();
        self.vertices.push(Vertex {
            position,
            normal,
            feature_id,
        });
        vertex_map.insert(key, index);
        index
    }

    fn optimize(&mut self) -> Result<()> {
        if self.vertices.is_empty() || self.indices.is_empty() {
            return Ok(());
        }

        let (vertex_count, remap) = generate_vertex_remap(&self.vertices, Some(&self.indices));
        let remapped_indices = remap_index_buffer(Some(&self.indices), vertex_count, &remap);
        let remapped_vertices = remap_vertex_buffer(&self.vertices, vertex_count, &remap);

        let mut optimized_indices =
            optimize_vertex_cache(&remapped_indices, remapped_vertices.len());
        optimize_overdraw_in_place_decoder(
            &mut optimized_indices,
            &remapped_vertices,
            OVERDRAW_THRESHOLD,
        );
        let optimized_vertices = optimize_vertex_fetch(&mut optimized_indices, &remapped_vertices);

        self.vertices = optimized_vertices;
        self.indices = optimized_indices;

        if self.vertices.len() > u32::MAX as usize {
            anyhow::bail!("GLB vertex count exceeds u32 index range");
        }

        Ok(())
    }
}

fn interpolate_normal(
    barycentric: [f32; 3],
    keys: [SmoothVertexKey; 3],
    source_vertex_normals: &HashMap<SmoothVertexKey, [f32; 3]>,
) -> [f32; 3] {
    let normals = keys.map(|key| {
        source_vertex_normals
            .get(&key)
            .copied()
            .unwrap_or([0.0, 0.0, 1.0])
    });
    [
        barycentric[0] * normals[0][0]
            + barycentric[1] * normals[1][0]
            + barycentric[2] * normals[2][0],
        barycentric[0] * normals[0][1]
            + barycentric[1] * normals[1][1]
            + barycentric[2] * normals[2][1],
        barycentric[0] * normals[0][2]
            + barycentric[1] * normals[1][2]
            + barycentric[2] * normals[2][2],
    ]
}

fn normalize_vector(vector: [f32; 3]) -> [f32; 3] {
    let length = (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
    if length > f32::EPSILON {
        [vector[0] / length, vector[1] / length, vector[2] / length]
    } else {
        [0.0, 0.0, 1.0]
    }
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

fn is_degenerate_source_triangle(points: [[f64; 3]; 3]) -> bool {
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
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let area_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
    area_sq <= 1.0e-18
}

fn is_degenerate_triangle(points: [[f32; 3]; 3]) -> bool {
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
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let area_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
    area_sq <= 1.0e-12
}

impl ProcessedScene {
    fn from_collector(collector: MeshCollector) -> Result<Self> {
        let MeshCollector {
            features,
            primitives,
            coordinate_transform,
            smooth_normals,
            clip_volume,
            ..
        } = collector;
        let mut primitives: BTreeMap<String, BuiltPrimitiveMesh> = primitives
            .into_iter()
            .map(|(feature_type, primitive)| {
                primitive
                    .build(&coordinate_transform, smooth_normals, clip_volume.as_ref())
                    .map(|built| (feature_type, built))
            })
            .collect::<Result<_>>()?;
        // Clipping (in `RawPrimitiveMesh::build`) can leave a per-feature_type
        // primitive with no surviving triangles when every triangle of that
        // type fell outside the clip volume. Drop those here so the encoder
        // never sees an empty primitive — `Bounds::from_vertices` would
        // otherwise return None and bail with "primitive bounds missing for
        // non-empty mesh", causing the whole tile to be pruned.
        primitives.retain(|_, primitive| {
            !(primitive.vertices.is_empty() && primitive.indices.is_empty())
        });
        let mut bounds = Bounds::empty();
        let mut has_vertices = false;
        for primitive in primitives.values() {
            for vertex in &primitive.vertices {
                bounds.add_point(vertex.position);
                has_vertices = true;
            }
        }

        let center = if has_vertices {
            bounds.center()
        } else {
            [0.0; 3]
        };
        for primitive in primitives.values_mut() {
            for vertex in &mut primitive.vertices {
                vertex.position[0] -= center[0];
                vertex.position[1] -= center[1];
                vertex.position[2] -= center[2];
            }
            primitive.optimize()?;
        }

        let mut final_bounds = Bounds::empty();
        let mut has_final_vertices = false;
        let mut processed_primitives = Vec::with_capacity(primitives.len());
        for (feature_type, primitive) in primitives {
            for vertex in &primitive.vertices {
                final_bounds.add_point(vertex.position);
                has_final_vertices = true;
            }
            processed_primitives.push(ProcessedPrimitiveMesh {
                feature_type,
                vertices: primitive.vertices,
                indices: primitive.indices,
            });
        }

        let processed_features = reorder_features_by_type(features, &mut processed_primitives)?;

        Ok(Self {
            primitives: processed_primitives,
            features: processed_features,
            center,
            node_translation_base: coordinate_transform.node_translation_base,
            bounds: has_final_vertices.then_some(final_bounds),
        })
    }

    fn vertex_count(&self) -> usize {
        self.primitives
            .iter()
            .map(|primitive| primitive.vertices.len())
            .sum()
    }

    fn index_count(&self) -> usize {
        self.primitives
            .iter()
            .map(|primitive| primitive.indices.len())
            .sum()
    }

    fn write_glb<P: AsRef<Path>>(self, output_path: P, options: &ExportOptions) -> Result<()> {
        if self
            .primitives
            .iter()
            .all(|primitive| primitive.vertices.is_empty())
        {
            debug!(
                "No geometry to write, creating empty GLB file at {}",
                output_path.as_ref().display()
            );
            if let Some(parent) = output_path.as_ref().parent() {
                std::fs::create_dir_all(parent).with_context(|| {
                    format!(
                        "Failed to create parent directory for {}",
                        output_path.as_ref().display()
                    )
                })?;
            }
            File::create(output_path.as_ref()).context("Create empty GLB file")?;
            return Ok(());
        }

        debug!(
            "Writing GLB file with {} vertices and {} indices to {}",
            self.vertex_count(),
            self.index_count(),
            output_path.as_ref().display()
        );

        let vertex_count = self.vertex_count();
        let index_count = self.index_count();
        let bounds = self
            .bounds
            .ok_or_else(|| anyhow::anyhow!("geometry bounds missing for non-empty mesh"))?;
        let node_translation = self.node_translation();
        let encoded = self.encode_glb(options)?;
        encoded.write(output_path.as_ref())?;

        debug!("GLB Summary: {}", output_path.as_ref().display());
        debug!("  Vertices: {vertex_count}");
        debug!("  Indices: {index_count}");
        debug!(
            "  Local coordinate range: X [{:.2}, {:.2}], Y [{:.2}, {:.2}], Z [{:.2}, {:.2}]",
            bounds.min[0],
            bounds.max[0],
            bounds.min[1],
            bounds.max[1],
            bounds.min[2],
            bounds.max[2]
        );
        debug!(
            "  Node translation: [{:.2}, {:.2}, {:.2}]",
            node_translation[0], node_translation[1], node_translation[2]
        );

        Ok(())
    }

    fn encode_glb(self, options: &ExportOptions) -> Result<EncodedGlb> {
        if options.quantize_geometry {
            QuantizedScene::from_processed(self)?.encode_glb(options)
        } else {
            self.encode_raw_glb(options)
        }
    }

    fn encode_raw_glb(self, options: &ExportOptions) -> Result<EncodedGlb> {
        let Self {
            primitives,
            features,
            node_translation_base,
            center,
            ..
        } = self;
        let mut buffer_builder = BufferBuilder::default();
        let primitive_encodings = encode_primitives_raw(
            &primitives,
            &mut buffer_builder,
            options.meshopt_compression,
        )?;
        let materials = build_materials(&primitive_encodings, options)?;
        let structural_metadata = StructuralMetadataExtension::from_features(
            &features,
            &mut buffer_builder,
            &options.metadata_class_name,
            options.meshopt_compression,
        )?;

        Ok(build_encoded_glb(
            buffer_builder,
            primitive_encodings,
            materials,
            build_node_matrix(
                1.0,
                [
                    node_translation_base[0] + center[0],
                    node_translation_base[1] + center[1],
                    node_translation_base[2] + center[2],
                ],
            ),
            false,
            options.meshopt_compression,
            structural_metadata.as_ref(),
        ))
    }

    fn node_translation(&self) -> [f32; 3] {
        [
            self.node_translation_base[0] + self.center[0],
            self.node_translation_base[1] + self.center[1],
            self.node_translation_base[2] + self.center[2],
        ]
    }
}

impl ProcessedPrimitiveMesh {
    fn select_index_buffer(&self) -> Result<IndexBuffer> {
        if self.vertices.len() <= (u16::MAX as usize) + 1 {
            let mut indices = Vec::with_capacity(self.indices.len());
            for index in &self.indices {
                indices.push(
                    u16::try_from(*index)
                        .context("GLB index exceeds u16 range after mesh optimization")?,
                );
            }
            Ok(IndexBuffer::U16(indices))
        } else {
            Ok(IndexBuffer::U32(self.indices.clone()))
        }
    }
}

impl QuantizedScene {
    fn from_processed(scene: ProcessedScene) -> Result<Self> {
        let bounds = scene
            .bounds
            .ok_or_else(|| anyhow::anyhow!("geometry bounds missing for non-empty mesh"))?;
        let position_scale = [
            bounds.min[0].abs(),
            bounds.max[0].abs(),
            bounds.min[1].abs(),
            bounds.max[1].abs(),
            bounds.min[2].abs(),
            bounds.max[2].abs(),
        ]
        .into_iter()
        .fold(0.0_f32, f32::max)
        .max(f32::EPSILON);

        let center = scene.node_translation();
        let primitives = scene
            .primitives
            .into_iter()
            .map(|primitive| primitive.quantize(position_scale))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            primitives,
            position_scale,
            center,
            features: scene.features,
        })
    }

    fn encode_glb(self, options: &ExportOptions) -> Result<EncodedGlb> {
        let Self {
            primitives,
            position_scale,
            center,
            features,
        } = self;
        let mut buffer_builder = BufferBuilder::default();
        let primitive_encodings = encode_quantized_primitives(
            &primitives,
            &mut buffer_builder,
            options.meshopt_compression,
        )?;
        let materials = build_materials(&primitive_encodings, options)?;
        let structural_metadata = StructuralMetadataExtension::from_features(
            &features,
            &mut buffer_builder,
            &options.metadata_class_name,
            options.meshopt_compression,
        )?;

        Ok(build_encoded_glb(
            buffer_builder,
            primitive_encodings,
            materials,
            build_node_matrix(position_scale, center),
            true,
            options.meshopt_compression,
            structural_metadata.as_ref(),
        ))
    }
}

impl ProcessedPrimitiveMesh {
    fn quantize(&self, position_scale: f32) -> Result<QuantizedPrimitiveMesh> {
        Bounds::from_vertices(&self.vertices)
            .ok_or_else(|| anyhow::anyhow!("primitive bounds missing for non-empty mesh"))?;

        let positions = self
            .vertices
            .iter()
            .map(|vertex| QuantizedPosition {
                position: [
                    quantize_position_component(vertex.position[0], position_scale),
                    quantize_position_component(vertex.position[1], position_scale),
                    quantize_position_component(vertex.position[2], position_scale),
                    0,
                ],
            })
            .collect::<Vec<_>>();
        let position_bounds = quantized_position_bounds(&positions)
            .ok_or_else(|| anyhow::anyhow!("quantized primitive has no positions"))?;

        Ok(QuantizedPrimitiveMesh {
            feature_type: self.feature_type.clone(),
            positions,
            normals: self
                .vertices
                .iter()
                .map(|vertex| QuantizedNormal {
                    normal: [
                        quantize_normal_component(vertex.normal[0]),
                        quantize_normal_component(vertex.normal[1]),
                        quantize_normal_component(vertex.normal[2]),
                        0,
                    ],
                })
                .collect(),
            feature_ids: self
                .vertices
                .iter()
                .map(|vertex| vertex.feature_id)
                .collect(),
            indices: self.select_index_buffer()?,
            position_bounds,
        })
    }
}

fn quantized_position_bounds(positions: &[QuantizedPosition]) -> Option<QuantizedBounds> {
    let first = positions.first()?;
    let mut min = [first.position[0], first.position[1], first.position[2]];
    let mut max = min;

    for position in &positions[1..] {
        for axis in 0..3 {
            min[axis] = min[axis].min(position.position[axis]);
            max[axis] = max[axis].max(position.position[axis]);
        }
    }

    Some(QuantizedBounds { min, max })
}

fn build_encoded_glb(
    buffer_builder: BufferBuilder,
    primitive_encodings: Vec<PrimitiveEncoding>,
    materials: Vec<json::Material>,
    node_matrix: [f32; 16],
    quantization_enabled: bool,
    meshopt_compression: bool,
    structural_metadata: Option<&StructuralMetadataExtension>,
) -> EncodedGlb {
    let has_structural_metadata = structural_metadata.is_some();
    let mesh = json::Mesh {
        primitives: primitive_encodings
            .into_iter()
            .map(|mut encoding| {
                encoding.primitive.extensions = Some(mesh_feature_extensions(
                    encoding.feature_count,
                    has_structural_metadata,
                ));
                encoding.primitive
            })
            .collect(),
        weights: None,
        extensions: None,
        extras: Option::default(),
        name: None,
    };

    let node = json::Node {
        mesh: Some(json::Index::new(0)),
        camera: None,
        children: None,
        skin: None,
        matrix: Some(node_matrix),
        rotation: None,
        scale: None,
        translation: None,
        weights: None,
        extensions: None,
        extras: Option::default(),
        name: None,
    };

    let scene = json::Scene {
        nodes: vec![json::Index::new(0)],
        extensions: None,
        extras: Option::default(),
        name: None,
    };

    let mut extensions_used = Vec::new();
    let mut extensions_required = Vec::new();
    if quantization_enabled {
        extensions_used.push(QUANTIZATION_EXTENSION.to_string());
        extensions_required.push(QUANTIZATION_EXTENSION.to_string());
    }
    if meshopt_compression {
        extensions_used.push(MESHOPT_EXTENSION.to_string());
        extensions_required.push(MESHOPT_EXTENSION.to_string());
    }
    if !mesh.primitives.is_empty() {
        extensions_used.push(MESH_FEATURES_EXTENSION.to_string());
    }
    if structural_metadata.is_some() {
        extensions_used.push(STRUCTURAL_METADATA_EXTENSION.to_string());
    }

    let BufferBuilder {
        bytes,
        mut buffer_views,
        accessors,
        meshopt_views,
        fallback_buffer_length,
    } = buffer_builder;
    let mut buffers = vec![json::Buffer {
        byte_length: json::validation::USize64(bytes.len() as u64),
        uri: None,
        name: Some("buffer0".into()),
        extensions: None,
        extras: Option::default(),
    }];
    if meshopt_compression {
        buffers.push(json::Buffer {
            byte_length: json::validation::USize64(fallback_buffer_length as u64),
            uri: None,
            name: Some("fallback".into()),
            extensions: Some(meshopt_fallback_buffer_extension()),
            extras: Option::default(),
        });
    }

    for (buffer_view, meshopt_view) in buffer_views.iter_mut().zip(meshopt_views.iter()) {
        let Some(meshopt_view) = meshopt_view else {
            continue;
        };
        buffer_view.extensions = Some(meshopt_buffer_view_extension(meshopt_view));
    }

    let root_extensions = structural_metadata.map(structural_metadata_root_extension);

    EncodedGlb {
        root: json::Root {
            accessors,
            buffers,
            buffer_views,
            materials,
            meshes: vec![mesh],
            nodes: vec![node],
            scenes: vec![scene],
            scene: Some(json::Index::new(0)),
            extensions: root_extensions,
            extensions_used,
            extensions_required,
            asset: json::Asset {
                version: GLTF_VERSION.into(),
                generator: Some("cityjson-convert".into()),
                copyright: None,
                ..Default::default()
            },
            ..Default::default()
        },
        bin_buffer: bytes,
    }
}

fn build_node_matrix(scale: f32, translation: [f32; 3]) -> [f32; 16] {
    [
        scale,
        0.0,
        0.0,
        0.0, //
        0.0,
        0.0,
        -scale,
        0.0, //
        0.0,
        scale,
        0.0,
        0.0, //
        translation[0],
        translation[2],
        -translation[1],
        1.0,
    ]
}

fn reorder_features_by_type(
    mut features: Vec<FeatureRecord>,
    primitives: &mut [ProcessedPrimitiveMesh],
) -> Result<Vec<FeatureRecord>> {
    let mut indexed = features
        .drain(..)
        .enumerate()
        .collect::<Vec<(usize, FeatureRecord)>>();
    indexed.sort_by(|(_, lhs), (_, rhs)| {
        lhs.feature_type
            .cmp(&rhs.feature_type)
            .then_with(|| lhs.object_id.cmp(&rhs.object_id))
    });

    let mut remap = vec![0u32; indexed.len()];
    let mut reordered = Vec::with_capacity(indexed.len());
    for (new_index, (old_index, feature)) in indexed.into_iter().enumerate() {
        remap[old_index] = u32::try_from(new_index).expect("feature index within u32 range");
        reordered.push(feature);
    }

    for primitive in primitives {
        for vertex in &mut primitive.vertices {
            let Some(mapped) = remap.get(vertex.feature_id as usize) else {
                anyhow::bail!("feature ID remap out of range");
            };
            vertex.feature_id = *mapped;
        }
    }

    Ok(reordered)
}

fn encode_primitives_raw(
    primitives: &[ProcessedPrimitiveMesh],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<Vec<PrimitiveEncoding>> {
    let mut primitive_encodings = Vec::with_capacity(primitives.len());
    for (material_index, primitive) in primitives.iter().enumerate() {
        let bounds = Bounds::from_vertices(&primitive.vertices)
            .ok_or_else(|| anyhow::anyhow!("primitive bounds missing for non-empty mesh"))?;
        let vertex_accessors = buffer_builder.push_float_vertices(
            &primitive.vertices,
            &bounds,
            meshopt_compression,
            true,
        )?;
        let index_buffer = primitive.select_index_buffer()?;
        let index_accessor = buffer_builder.push_indices(
            &index_buffer,
            primitive.vertices.len(),
            meshopt_compression,
        )?;
        primitive_encodings.push(PrimitiveEncoding {
            feature_type: primitive.feature_type.clone(),
            primitive: build_primitive(
                vertex_accessors,
                index_accessor,
                u32::try_from(material_index).expect("material index within u32 range"),
            ),
            feature_count: primitive
                .vertices
                .iter()
                .map(|vertex| vertex.feature_id)
                .collect::<BTreeSet<_>>()
                .len(),
        });
    }
    Ok(primitive_encodings)
}

fn encode_quantized_primitives(
    primitives: &[QuantizedPrimitiveMesh],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<Vec<PrimitiveEncoding>> {
    let mut primitive_encodings = Vec::with_capacity(primitives.len());
    for (material_index, primitive) in primitives.iter().enumerate() {
        let vertex_accessors = buffer_builder.push_quantized_vertices(
            &primitive.positions,
            &primitive.normals,
            &primitive.feature_ids,
            &primitive.position_bounds,
            meshopt_compression,
        )?;
        let index_accessor = buffer_builder.push_indices(
            &primitive.indices,
            primitive.positions.len(),
            meshopt_compression,
        )?;
        primitive_encodings.push(PrimitiveEncoding {
            feature_type: primitive.feature_type.clone(),
            primitive: build_primitive(
                vertex_accessors,
                index_accessor,
                u32::try_from(material_index).expect("material index within u32 range"),
            ),
            feature_count: primitive
                .feature_ids
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
                .len(),
        });
    }
    Ok(primitive_encodings)
}

fn build_primitive(
    vertex_accessors: VertexAccessors,
    index_accessor: json::Index<json::Accessor>,
    material_index: u32,
) -> json::mesh::Primitive {
    let mut attributes = BTreeMap::new();
    attributes.insert(
        json::validation::Checked::Valid(json::mesh::Semantic::Positions),
        vertex_accessors.positions,
    );
    attributes.insert(
        json::validation::Checked::Valid(json::mesh::Semantic::Normals),
        vertex_accessors.normals,
    );
    if let Some(feature_ids) = vertex_accessors.feature_ids {
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Extras(
                "FEATURE_ID_0".to_string(),
            )),
            feature_ids,
        );
    }

    json::mesh::Primitive {
        attributes,
        indices: Some(index_accessor),
        material: Some(json::Index::new(material_index)),
        mode: json::validation::Checked::Valid(json::mesh::Mode::Triangles),
        targets: None,
        extensions: None,
        extras: Option::default(),
    }
}

fn build_materials(
    primitive_encodings: &[PrimitiveEncoding],
    options: &ExportOptions,
) -> Result<Vec<json::Material>> {
    primitive_encodings
        .iter()
        .map(|encoding| {
            create_default_material(resolve_feature_type_color(&encoding.feature_type, options))
        })
        .collect()
}

fn resolve_feature_type_color<'a>(feature_type: &str, options: &'a ExportOptions) -> &'a str {
    if let Some(color) = options.feature_type_colors.get(feature_type) {
        return color.as_str();
    }

    default_color_for_feature_type(feature_type).unwrap_or(&options.native_glb_color)
}

fn default_color_for_feature_type(feature_type: &str) -> Option<&'static str> {
    match feature_type {
        "Building" => Some("#d8c3a5"),
        "BuildingPart" => Some("#e6d5b8"),
        "BuildingInstallation" => Some("#b89b7a"),
        "TINRelief" => Some("#9fb37c"),
        "Road" => Some("#8b8b8b"),
        "Railway" => Some("#666666"),
        "TransportSquare" => Some("#a8a8a8"),
        "WaterBody" => Some("#7db8da"),
        "PlantCover" => Some("#8dbb6b"),
        "SolitaryVegetationObject" => Some("#689d45"),
        "LandUse" => Some("#c0cf8c"),
        "CityFurniture" => Some("#c78d5b"),
        "Bridge" => Some("#b39a86"),
        "BridgePart" => Some("#c2aa96"),
        "BridgeInstallation" => Some("#927562"),
        "BridgeConstructiveElement" => Some("#8e705d"),
        "Tunnel" => Some("#968d84"),
        "TunnelPart" => Some("#a39a90"),
        "TunnelInstallation" => Some("#7b726a"),
        "GenericCityObject" => Some("#b48ead"),
        "OtherConstruction" => Some("#b5947d"),
        _ => None,
    }
}

fn build_structural_metadata_columns(
    features: &[FeatureRecord],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<BTreeMap<String, StructuralMetadataColumn>> {
    let mut column_types = BTreeMap::<String, &'static str>::new();
    for feature in features {
        for (name, value) in &feature.attributes {
            let kind = match value {
                MetadataValue::Bool(_) => "bool",
                MetadataValue::Int(_) => "int",
                MetadataValue::Float(_) => "float",
                MetadataValue::String(_) => "string",
            };
            column_types
                .entry(name.clone())
                .and_modify(|existing| {
                    if *existing != kind {
                        *existing = if *existing == "string" || kind == "string" {
                            "string"
                        } else if *existing == "float" || kind == "float" {
                            "float"
                        } else {
                            "int"
                        };
                    }
                })
                .or_insert(kind);
        }
    }

    let mut columns = BTreeMap::new();
    for (name, kind) in column_types {
        let column = match kind {
            "bool" => build_bool_metadata_column(&name, features, buffer_builder),
            "int" => {
                build_int_metadata_column(&name, features, buffer_builder, meshopt_compression)?
            }
            "float" => {
                build_float_metadata_column(&name, features, buffer_builder, meshopt_compression)?
            }
            "string" => {
                build_string_metadata_column(&name, features, buffer_builder, meshopt_compression)?
            }
            _ => continue,
        };
        columns.insert(name, column);
    }

    Ok(columns)
}

fn build_bool_metadata_column(
    name: &str,
    features: &[FeatureRecord],
    buffer_builder: &mut BufferBuilder,
) -> StructuralMetadataColumn {
    let values = features
        .iter()
        .map(|feature| match feature.attributes.get(name) {
            Some(MetadataValue::Bool(value)) => i8::from(*value),
            _ => i8::MAX,
        })
        .collect::<Vec<_>>();
    let view = buffer_builder.push_scalar_buffer_view(&values, json::buffer::Target::ArrayBuffer);
    StructuralMetadataColumn {
        property: json_value!({
            "type": "SCALAR",
            "componentType": "INT8",
            "noData": i8::MAX,
        }),
        property_table_entry: json_value!({
            "values": view.value(),
        }),
    }
}

fn build_int_metadata_column(
    name: &str,
    features: &[FeatureRecord],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<StructuralMetadataColumn> {
    let values = features
        .iter()
        .map(|feature| match feature.attributes.get(name) {
            Some(MetadataValue::Bool(value)) => i32::from(*value),
            Some(MetadataValue::Int(value)) => *value,
            #[allow(clippy::cast_possible_truncation)]
            Some(MetadataValue::Float(value)) => *value as i32,
            _ => i32::MAX,
        })
        .collect::<Vec<_>>();
    let view = buffer_builder.push_metadata_scalar_buffer_view(&values, meshopt_compression)?;
    Ok(StructuralMetadataColumn {
        property: json_value!({
            "type": "SCALAR",
            "componentType": "INT32",
            "noData": i32::MAX,
        }),
        property_table_entry: json_value!({
            "values": view.value(),
        }),
    })
}

fn build_float_metadata_column(
    name: &str,
    features: &[FeatureRecord],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<StructuralMetadataColumn> {
    let values = features
        .iter()
        .map(|feature| match feature.attributes.get(name) {
            Some(MetadataValue::Bool(value)) => f32::from(u8::from(*value)),
            #[allow(clippy::cast_precision_loss)]
            Some(MetadataValue::Int(value)) => *value as f32,
            Some(MetadataValue::Float(value)) => *value,
            _ => f32::MAX,
        })
        .collect::<Vec<_>>();
    let view = buffer_builder.push_metadata_scalar_buffer_view(&values, meshopt_compression)?;
    Ok(StructuralMetadataColumn {
        property: json_value!({
            "type": "SCALAR",
            "componentType": "FLOAT32",
            "noData": f32::MAX,
        }),
        property_table_entry: json_value!({
            "values": view.value(),
        }),
    })
}

fn build_string_metadata_column(
    name: &str,
    features: &[FeatureRecord],
    buffer_builder: &mut BufferBuilder,
    meshopt_compression: bool,
) -> Result<StructuralMetadataColumn> {
    let mut values = Vec::<u8>::new();
    let mut offsets = Vec::<u32>::with_capacity(features.len() + 1);
    offsets.push(0);
    for feature in features {
        if let Some(MetadataValue::String(value)) = feature.attributes.get(name) {
            values.extend_from_slice(value.as_bytes());
        }
        offsets.push(u32::try_from(values.len()).expect("string column offset within u32 range"));
    }

    let values_view =
        buffer_builder.push_byte_buffer_view(&values, json::buffer::Target::ArrayBuffer);
    let offsets_view =
        buffer_builder.push_metadata_scalar_buffer_view(&offsets, meshopt_compression)?;
    Ok(StructuralMetadataColumn {
        property: json_value!({
            "type": "STRING",
        }),
        property_table_entry: json_value!({
            "values": values_view.value(),
            "stringOffsets": offsets_view.value(),
            "stringOffsetType": "UINT32",
        }),
    })
}

impl StructuralMetadataExtension {
    fn from_features(
        features: &[FeatureRecord],
        buffer_builder: &mut BufferBuilder,
        class_name: &str,
        meshopt_compression: bool,
    ) -> Result<Option<Self>> {
        if features.is_empty() {
            return Ok(None);
        }

        let columns =
            build_structural_metadata_columns(features, buffer_builder, meshopt_compression)?;
        if columns.is_empty() {
            return Ok(None);
        }

        Ok(Some(Self {
            class_name: class_name.to_string(),
            columns,
            feature_count: features.len(),
        }))
    }
}

impl BufferBuilder {
    fn push_quantized_vertices(
        &mut self,
        positions: &[QuantizedPosition],
        normals: &[QuantizedNormal],
        feature_ids: &[u32],
        bounds: &QuantizedBounds,
        meshopt_compression: bool,
    ) -> Result<VertexAccessors> {
        let position_view = if meshopt_compression {
            let position_bytes = encode_vertex_buffer(positions)
                .context("failed to meshopt-encode quantized position stream")?;
            self.push_meshopt_buffer_view(
                &position_bytes,
                positions.len() * QUANTIZED_POSITION_STRIDE,
                Some(QUANTIZED_POSITION_STRIDE),
                positions.len(),
                json::buffer::Target::ArrayBuffer,
                "ATTRIBUTES",
                None,
                std::mem::align_of::<i16>(),
            )
        } else {
            self.push_buffer_view(
                positions,
                Some(QUANTIZED_POSITION_STRIDE),
                json::buffer::Target::ArrayBuffer,
            )
        };
        let positions = self.push_accessor(json::Accessor {
            buffer_view: Some(position_view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(positions.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::I16,
            )),
            normalized: true,
            min: Some(json::Value::Array(
                bounds.min.iter().copied().map(json::Value::from).collect(),
            )),
            max: Some(json::Value::Array(
                bounds.max.iter().copied().map(json::Value::from).collect(),
            )),
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: None,
            extras: Option::default(),
            name: None,
            sparse: None,
        });

        let normal_view = if meshopt_compression {
            let normal_bytes = encode_vertex_buffer(normals)
                .context("failed to meshopt-encode quantized normal stream")?;
            self.push_meshopt_buffer_view(
                &normal_bytes,
                normals.len() * QUANTIZED_NORMAL_STRIDE,
                Some(QUANTIZED_NORMAL_STRIDE),
                normals.len(),
                json::buffer::Target::ArrayBuffer,
                "ATTRIBUTES",
                None,
                std::mem::align_of::<i8>(),
            )
        } else {
            self.push_buffer_view(
                normals,
                Some(QUANTIZED_NORMAL_STRIDE),
                json::buffer::Target::ArrayBuffer,
            )
        };
        let normals = self.push_accessor(json::Accessor {
            buffer_view: Some(normal_view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(normals.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::I8,
            )),
            normalized: true,
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: None,
            extras: Option::default(),
            min: None,
            max: None,
            name: None,
            sparse: None,
        });

        let feature_ids = Some(self.push_feature_ids(feature_ids, meshopt_compression)?);

        Ok(VertexAccessors {
            positions,
            normals,
            feature_ids,
        })
    }

    fn push_float_vertices(
        &mut self,
        vertices: &[Vertex],
        bounds: &Bounds,
        meshopt_compression: bool,
        include_feature_ids: bool,
    ) -> Result<VertexAccessors> {
        let positions: Vec<FloatPosition> = vertices
            .iter()
            .map(|vertex| FloatPosition {
                position: vertex.position,
            })
            .collect();
        let normals: Vec<FloatNormal> = vertices
            .iter()
            .map(|vertex| FloatNormal {
                normal: vertex.normal,
            })
            .collect();

        let position_stride = std::mem::size_of::<FloatPosition>();
        let normal_stride = std::mem::size_of::<FloatNormal>();

        let position_view = if meshopt_compression {
            let position_bytes = encode_vertex_buffer(&positions)
                .context("failed to meshopt-encode float position stream")?;
            self.push_meshopt_buffer_view(
                &position_bytes,
                positions.len() * position_stride,
                Some(position_stride),
                positions.len(),
                json::buffer::Target::ArrayBuffer,
                "ATTRIBUTES",
                None,
                std::mem::align_of::<f32>(),
            )
        } else {
            self.push_buffer_view(
                &positions,
                Some(position_stride),
                json::buffer::Target::ArrayBuffer,
            )
        };
        let positions = self.push_accessor(json::Accessor {
            buffer_view: Some(position_view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(positions.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::F32,
            )),
            normalized: false,
            min: Some(json::Value::Array(
                bounds.min.iter().copied().map(json::Value::from).collect(),
            )),
            max: Some(json::Value::Array(
                bounds.max.iter().copied().map(json::Value::from).collect(),
            )),
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: None,
            extras: Option::default(),
            name: None,
            sparse: None,
        });

        let normal_view = if meshopt_compression {
            let normal_bytes = encode_vertex_buffer(&normals)
                .context("failed to meshopt-encode float normal stream")?;
            self.push_meshopt_buffer_view(
                &normal_bytes,
                normals.len() * normal_stride,
                Some(normal_stride),
                normals.len(),
                json::buffer::Target::ArrayBuffer,
                "ATTRIBUTES",
                None,
                std::mem::align_of::<f32>(),
            )
        } else {
            self.push_buffer_view(
                &normals,
                Some(normal_stride),
                json::buffer::Target::ArrayBuffer,
            )
        };
        let normals = self.push_accessor(json::Accessor {
            buffer_view: Some(normal_view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(normals.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::F32,
            )),
            normalized: false,
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: None,
            extras: Option::default(),
            min: None,
            max: None,
            name: None,
            sparse: None,
        });

        let feature_ids = if include_feature_ids {
            let feature_ids = vertices
                .iter()
                .map(|vertex| vertex.feature_id)
                .collect::<Vec<_>>();
            Some(self.push_feature_ids(&feature_ids, meshopt_compression)?)
        } else {
            None
        };

        Ok(VertexAccessors {
            positions,
            normals,
            feature_ids,
        })
    }

    fn push_feature_ids(
        &mut self,
        feature_ids: &[u32],
        meshopt_compression: bool,
    ) -> Result<json::Index<json::Accessor>> {
        let feature_ids = FeatureIdBuffer::from_feature_ids(feature_ids)?;
        let view = match &feature_ids {
            FeatureIdBuffer::U8(feature_id_values) => {
                if meshopt_compression {
                    let padded_feature_ids = feature_id_values
                        .iter()
                        .copied()
                        .map(|feature_id| PaddedFeatureIdU8 {
                            feature_id,
                            _padding: [0; 3],
                        })
                        .collect::<Vec<_>>();
                    let encoded = encode_vertex_buffer(&padded_feature_ids)
                        .context("failed to meshopt-encode feature ID stream")?;
                    self.push_meshopt_buffer_view(
                        &encoded,
                        padded_feature_ids.len() * FeatureIdBuffer::meshopt_byte_stride(),
                        Some(FeatureIdBuffer::meshopt_byte_stride()),
                        feature_id_values.len(),
                        json::buffer::Target::ArrayBuffer,
                        "ATTRIBUTES",
                        None,
                        std::mem::align_of::<PaddedFeatureIdU8>(),
                    )
                } else {
                    self.push_buffer_view(
                        feature_id_values,
                        Some(feature_ids.byte_stride()),
                        json::buffer::Target::ArrayBuffer,
                    )
                }
            }
            FeatureIdBuffer::U16(feature_id_values) => {
                if meshopt_compression {
                    let padded_feature_ids = feature_id_values
                        .iter()
                        .copied()
                        .map(|feature_id| PaddedFeatureIdU16 {
                            feature_id,
                            _padding: 0,
                        })
                        .collect::<Vec<_>>();
                    let encoded = encode_vertex_buffer(&padded_feature_ids)
                        .context("failed to meshopt-encode feature ID stream")?;
                    self.push_meshopt_buffer_view(
                        &encoded,
                        padded_feature_ids.len() * FeatureIdBuffer::meshopt_byte_stride(),
                        Some(FeatureIdBuffer::meshopt_byte_stride()),
                        feature_id_values.len(),
                        json::buffer::Target::ArrayBuffer,
                        "ATTRIBUTES",
                        None,
                        std::mem::align_of::<PaddedFeatureIdU16>(),
                    )
                } else {
                    self.push_buffer_view(
                        feature_id_values,
                        Some(feature_ids.byte_stride()),
                        json::buffer::Target::ArrayBuffer,
                    )
                }
            }
        };

        Ok(self.push_accessor(json::Accessor {
            buffer_view: Some(view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(feature_ids.count() as u64),
            component_type: feature_ids.component_type(),
            normalized: false,
            min: Some(json::Value::from(vec![0])),
            max: Some(json::Value::from(vec![feature_ids.max_value()])),
            type_: json::validation::Checked::Valid(json::accessor::Type::Scalar),
            extensions: None,
            extras: Option::default(),
            name: None,
            sparse: None,
        }))
    }

    fn push_indices(
        &mut self,
        index_buffer: &IndexBuffer,
        vertex_count: usize,
        meshopt_compression: bool,
    ) -> Result<json::Index<json::Accessor>> {
        let view = if meshopt_compression {
            let encoded_indices = encode_index_buffer(&index_buffer.as_u32_vec(), vertex_count)
                .context("failed to meshopt-encode index stream")?;
            self.push_meshopt_buffer_view(
                &encoded_indices,
                index_buffer.byte_length(),
                None,
                index_buffer.count(),
                json::buffer::Target::ElementArrayBuffer,
                "TRIANGLES",
                None,
                index_buffer.byte_stride(),
            )
        } else {
            match index_buffer {
                IndexBuffer::U16(indices) => {
                    self.push_buffer_view(indices, None, json::buffer::Target::ElementArrayBuffer)
                }
                IndexBuffer::U32(indices) => {
                    self.push_buffer_view(indices, None, json::buffer::Target::ElementArrayBuffer)
                }
            }
        };
        Ok(self.push_accessor(json::Accessor {
            buffer_view: Some(view),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(index_buffer.count() as u64),
            component_type: index_buffer.component_type(),
            normalized: false,
            min: Some(json::Value::from(vec![0])),
            max: Some(json::Value::from(vec![index_buffer.max_value()])),
            type_: json::validation::Checked::Valid(json::accessor::Type::Scalar),
            extensions: None,
            extras: Option::default(),
            name: None,
            sparse: None,
        }))
    }

    fn push_buffer_view<T>(
        &mut self,
        data: &[T],
        byte_stride: Option<usize>,
        target: json::buffer::Target,
    ) -> json::Index<json::buffer::View> {
        let byte_offset = self.bytes.len();
        let byte_length = std::mem::size_of_val(data);
        let data_bytes =
            unsafe { std::slice::from_raw_parts(data.as_ptr().cast::<u8>(), byte_length) };
        self.bytes.extend_from_slice(data_bytes);

        let index = self.buffer_views.len();
        self.buffer_views.push(json::buffer::View {
            buffer: json::Index::new(0),
            byte_length: json::validation::USize64(byte_length as u64),
            byte_offset: Some(json::validation::USize64(byte_offset as u64)),
            byte_stride: byte_stride.map(json::buffer::Stride),
            target: Some(json::validation::Checked::Valid(target)),
            extensions: None,
            extras: Option::default(),
            name: None,
        });
        self.meshopt_views.push(None);

        json::Index::new(u32::try_from(index).expect("buffer view index exceeds u32 range"))
    }

    fn push_scalar_buffer_view<T>(
        &mut self,
        data: &[T],
        target: json::buffer::Target,
    ) -> json::Index<json::buffer::View> {
        self.push_buffer_view(data, None, target)
    }

    fn push_metadata_scalar_buffer_view<T>(
        &mut self,
        data: &[T],
        meshopt_compression: bool,
    ) -> Result<json::Index<json::buffer::View>> {
        if meshopt_compression && std::mem::size_of::<T>() % 4 == 0 {
            let encoded = encode_vertex_buffer(data)
                .context("failed to meshopt-encode structural metadata column")?;
            Ok(self.push_meshopt_buffer_view(
                &encoded,
                std::mem::size_of_val(data),
                Some(std::mem::size_of::<T>()),
                data.len(),
                json::buffer::Target::ArrayBuffer,
                "ATTRIBUTES",
                None,
                std::mem::align_of::<T>(),
            ))
        } else {
            Ok(self.push_scalar_buffer_view(data, json::buffer::Target::ArrayBuffer))
        }
    }

    fn push_byte_buffer_view(
        &mut self,
        data: &[u8],
        target: json::buffer::Target,
    ) -> json::Index<json::buffer::View> {
        self.push_buffer_view(data, None, target)
    }

    fn push_meshopt_buffer_view(
        &mut self,
        compressed_data: &[u8],
        fallback_byte_length: usize,
        byte_stride: Option<usize>,
        count: usize,
        target: json::buffer::Target,
        mode: &'static str,
        filter: Option<&'static str>,
        fallback_alignment: usize,
    ) -> json::Index<json::buffer::View> {
        let compressed_byte_offset = self.bytes.len();
        let compressed_byte_length = compressed_data.len();
        self.bytes.extend_from_slice(compressed_data);

        let fallback_byte_offset = align_length(self.fallback_buffer_length, fallback_alignment);
        self.fallback_buffer_length = fallback_byte_offset + fallback_byte_length;

        let index = self.buffer_views.len();
        self.buffer_views.push(json::buffer::View {
            buffer: json::Index::new(1),
            byte_length: json::validation::USize64(fallback_byte_length as u64),
            byte_offset: Some(json::validation::USize64(fallback_byte_offset as u64)),
            byte_stride: byte_stride.map(json::buffer::Stride),
            target: Some(json::validation::Checked::Valid(target)),
            extensions: None,
            extras: Option::default(),
            name: None,
        });
        self.meshopt_views.push(Some(MeshoptBufferView {
            buffer: 0,
            byte_offset: compressed_byte_offset as u64,
            byte_length: compressed_byte_length as u64,
            byte_stride: byte_stride.unwrap_or(fallback_alignment) as u64,
            count: count as u64,
            mode,
            filter,
        }));

        json::Index::new(u32::try_from(index).expect("buffer view index exceeds u32 range"))
    }

    fn push_accessor(&mut self, accessor: json::Accessor) -> json::Index<json::Accessor> {
        let index = self.accessors.len();
        self.accessors.push(accessor);
        json::Index::new(u32::try_from(index).expect("accessor index exceeds u32 range"))
    }
}

impl EncodedGlb {
    fn write<P: AsRef<Path>>(self, output_path: P) -> Result<()> {
        let mut json_bytes =
            serde_json::to_vec(&self.root).context("failed to encode glTF JSON chunk")?;
        let json_padding = (4 - (json_bytes.len() % 4)) % 4;
        json_bytes.extend(std::iter::repeat_n(b' ', json_padding));

        let mut bin_buffer = self.bin_buffer;
        let bin_padding = (4 - (bin_buffer.len() % 4)) % 4;
        bin_buffer.extend(std::iter::repeat_n(0u8, bin_padding));

        let total_length = 12 + 8 + json_bytes.len() + 8 + bin_buffer.len();
        let total_length_u32 =
            u32::try_from(total_length).context("GLB total length exceeds u32 range")?;
        let json_length_u32 =
            u32::try_from(json_bytes.len()).context("GLB JSON chunk exceeds u32 range")?;
        let bin_length_u32 =
            u32::try_from(bin_buffer.len()).context("GLB BIN chunk exceeds u32 range")?;

        let mut glb_bytes = Vec::with_capacity(total_length);
        glb_bytes.extend_from_slice(b"glTF");
        glb_bytes.extend_from_slice(&2u32.to_le_bytes());
        glb_bytes.extend_from_slice(&total_length_u32.to_le_bytes());

        glb_bytes.extend_from_slice(&json_length_u32.to_le_bytes());
        glb_bytes.extend_from_slice(b"JSON");
        glb_bytes.extend_from_slice(&json_bytes);

        glb_bytes.extend_from_slice(&bin_length_u32.to_le_bytes());
        glb_bytes.extend_from_slice(b"BIN\0");
        glb_bytes.extend_from_slice(&bin_buffer);

        if let Some(parent) = output_path.as_ref().parent() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Failed to create parent directory for {}",
                    output_path.as_ref().display()
                )
            })?;
        }

        let mut file = File::create(output_path.as_ref())?;
        file.write_all(&glb_bytes)?;

        Ok(())
    }
}

fn structural_metadata_json(structural_metadata: &StructuralMetadataExtension) -> JsonValue {
    let mut class_properties = JsonMap::new();
    let mut property_table_properties = JsonMap::new();
    for (name, column) in &structural_metadata.columns {
        class_properties.insert(name.clone(), column.property.clone());
        property_table_properties.insert(name.clone(), column.property_table_entry.clone());
    }

    json_value!({
        "schema": {
            "id": "schema_0",
            "classes": {
                structural_metadata.class_name.clone(): {
                    "name": structural_metadata.class_name,
                    "properties": JsonValue::Object(class_properties),
                }
            }
        },
        "propertyTables": [{
            "class": structural_metadata.class_name,
            "count": structural_metadata.feature_count,
            "properties": JsonValue::Object(property_table_properties),
        }]
    })
}

fn structural_metadata_root_extension(
    structural_metadata: &StructuralMetadataExtension,
) -> json::extensions::root::Root {
    let mut others = JsonMap::new();
    others.insert(
        STRUCTURAL_METADATA_EXTENSION.to_string(),
        structural_metadata_json(structural_metadata),
    );
    json::extensions::root::Root { others }
}

fn meshopt_fallback_buffer_extension() -> json::extensions::buffer::Buffer {
    let mut others = JsonMap::new();
    others.insert(
        MESHOPT_EXTENSION.to_string(),
        json_value!({ "fallback": true }),
    );
    json::extensions::buffer::Buffer { others }
}

fn meshopt_buffer_view_extension(
    meshopt_view: &MeshoptBufferView,
) -> json::extensions::buffer::View {
    let mut extension = JsonMap::new();
    extension.insert("buffer".into(), JsonValue::from(meshopt_view.buffer));
    extension.insert(
        "byteOffset".into(),
        JsonValue::from(meshopt_view.byte_offset),
    );
    extension.insert(
        "byteLength".into(),
        JsonValue::from(meshopt_view.byte_length),
    );
    extension.insert(
        "byteStride".into(),
        JsonValue::from(meshopt_view.byte_stride),
    );
    extension.insert("count".into(), JsonValue::from(meshopt_view.count));
    extension.insert("mode".into(), JsonValue::from(meshopt_view.mode));
    if let Some(filter) = meshopt_view.filter {
        extension.insert("filter".into(), JsonValue::from(filter));
    }

    let mut others = JsonMap::new();
    others.insert(MESHOPT_EXTENSION.to_string(), JsonValue::Object(extension));
    json::extensions::buffer::View { others }
}

fn mesh_feature_extensions(
    feature_count: usize,
    has_structural_metadata: bool,
) -> json::extensions::mesh::Primitive {
    let mut feature_id = json_value!({
        "featureCount": feature_count,
        "attribute": 0,
    });
    if has_structural_metadata {
        feature_id
            .as_object_mut()
            .expect("feature ID JSON should be an object")
            .insert("propertyTable".to_string(), json_value!(0));
    }

    let mut others = JsonMap::new();
    others.insert(
        MESH_FEATURES_EXTENSION.to_string(),
        json_value!({
            "featureIds": [feature_id]
        }),
    );
    json::extensions::mesh::Primitive { others }
}

fn align_length(length: usize, alignment: usize) -> usize {
    if alignment <= 1 {
        return length;
    }
    let remainder = length % alignment;
    if remainder == 0 {
        length
    } else {
        length + (alignment - remainder)
    }
}
