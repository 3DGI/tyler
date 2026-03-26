use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};
use earcutr::earcut;
use gltf::json as json;

use crate::cli::TilesVersion;
use crate::material::MaterialConfig;
use crate::parser::{CityJSONFeatureVertices, CityObjectType, Geometry, Transform, World};
use crate::proj::Proj;
use crate::spatial_structs::{QuadTree, QuadTreeNodeId};

const GLTF_VERSION: &str = "2.0";

/// Create PBR material with specified metallic and roughness parameters.
/// Base color is white so per-vertex COLOR_0 provides the actual coloring.
fn create_material(base_color: &str, metallic: f32, roughness: f32) -> Result<json::Material, anyhow::Error> {
    let base_color_rgba = crate::material::hex_to_rgba(base_color)?;
    
    Ok(json::Material {
        name: None,
        extensions: Default::default(),
        extras: Default::default(),
        pbr_metallic_roughness: json::material::PbrMetallicRoughness {
            base_color_factor: json::material::PbrBaseColorFactor(base_color_rgba),
            metallic_factor: json::material::StrengthFactor(metallic),
            roughness_factor: json::material::StrengthFactor(roughness),
            base_color_texture: None,
            metallic_roughness_texture: None,
            extensions: Default::default(),
            extras: Default::default(),
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

pub fn write_tile_glb<P: AsRef<Path>>(
    world: &World,
    quadtree: &QuadTree,
    qtree_node_id: QuadTreeNodeId,
    output_path: P,
    material_config: &MaterialConfig,
    tiles_version: TilesVersion,
) -> Result<()> {
    let qtree_node = quadtree
        .node(&qtree_node_id)
        .context("Tile not present in quadtree")?;

    let epsg_code = world
        .crs
        .to_epsg()
        .map_err(|e| anyhow::anyhow!("Failed to read EPSG code from metadata: {}", e))?;
    let crs_from = format!("EPSG:{}", epsg_code);

    // Transform coordinates to ECEF (EPSG:4978) to match root transform coordinate system
    // Root transform is in ECEF, so GLB content must also be in ECEF for correct positioning
    let transformer_to_ecef =
        Proj::new_known_crs(&crs_from, "EPSG:4978", None).context("Create CRS to ECEF transformer")?;
    
    // Also need transformer for vertical geoid correction (local to tile)
    let transformer_crs_to_ell =
        Proj::new_known_crs(&crs_from, "EPSG:4979", None).context("Create CRS to ellipsoidal transformer")?;

    // Calculate root center in input CRS, then transform to ECEF
    // GLB coordinates will be relative to root center in ECEF to match root transform
    let root_bbox = quadtree.bbox(&world.grid);
    let root_center_input_crs = [
        (root_bbox[0] + root_bbox[3]) * 0.5,
        (root_bbox[1] + root_bbox[4]) * 0.5,
        (root_bbox[2] + root_bbox[5]) * 0.5,
    ];
    // Transform root center to ECEF for GLB content coordinates
    let root_center_ecef = transformer_to_ecef
        .convert((root_center_input_crs[0], root_center_input_crs[1], root_center_input_crs[2]))
        .context("Transform root center to ECEF")?;
    

    // Use tile center for vertical geoid correction (local to tile)
    let tile_bbox = qtree_node.bbox(&world.grid);
    let tile_center_original = [
        (tile_bbox[0] + tile_bbox[3]) * 0.5,
        (tile_bbox[1] + tile_bbox[4]) * 0.5,
        tile_bbox[2],
    ];

    let vertical_geoid_n = transformer_crs_to_ell
        .convert((
            tile_center_original[0],
            tile_center_original[1],
            tile_center_original[2],
        ))
        .map(|(_, _, h_ell)| h_ell - tile_center_original[2])
        .unwrap_or(0.0);

    // For 3D Tiles with root transform, GLB coordinates must be in ECEF and relative to ROOT center in ECEF
    // This ensures coordinate system consistency: root transform (ECEF) + GLB content (ECEF) = correct positioning
    let mut builder = MeshBuilder::new(
        transformer_to_ecef,
        root_center_ecef,
        vertical_geoid_n,
        material_config.metallic_factor,
        material_config.roughness_factor,
    );

    // Deduplicate by feature id: a building can be referenced by multiple cells (e.g. bbox
    // intersection), so we add each feature at most once per tile to avoid duplicate geometry
    // and the same BuildingId appearing on multiple meshes.
    let mut seen_fids: HashSet<usize> = HashSet::new();
    for cellid in qtree_node.cells() {
        let cell = world.grid.cell(cellid);
        for fid in cell.feature_ids.iter() {
            if !seen_fids.insert(*fid) {
                continue; // already added this feature to this tile
            }
            let feature = &world.features[*fid];
            let cf = CityJSONFeatureVertices::from_file(&feature.path_jsonl)
                .map_err(|e| anyhow::anyhow!("Failed to read {:?}: {}", feature.path_jsonl, e))?;
            builder.add_feature(&cf, &world.transform, &material_config.color_map)?;
        }
    }

    builder.write_glb(output_path, tiles_version)
}

struct MeshBuilder {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    colors: Vec<[f32; 4]>,
    batch_ids: Vec<u32>,
    indices: Vec<u32>,
    next_batch_index: u32,
    /// CityJSON object id per batch index (batch_id_to_cityobject_id[i] = id for batch index i)
    batch_id_to_cityobject_id: Vec<String>,
    /// CityObjectType per batch index (e.g. "Building", "WaterBody")
    batch_id_to_cityobject_type: Vec<String>,
    transformer_to_ecef: Proj,
    root_center_ecef: (f64, f64, f64),
    vertical_bias: f64,
    metallic_factor: f32,
    roughness_factor: f32,
}

impl MeshBuilder {
    fn new(
        transformer_to_ecef: Proj,
        root_center_ecef: (f64, f64, f64),
        vertical_bias: f64,
        metallic_factor: f32,
        roughness_factor: f32,
    ) -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            colors: Vec::new(),
            batch_ids: Vec::new(),
            indices: Vec::new(),
            next_batch_index: 0,
            batch_id_to_cityobject_id: Vec::new(),
            batch_id_to_cityobject_type: Vec::new(),
            transformer_to_ecef,
            root_center_ecef,
            vertical_bias,
            metallic_factor,
            roughness_factor,
        }
    }

    /// Add one feature (one CityJSON feature = one building). Uses a single batch index and the
    /// Building id (feature id or Building CityObject key), not BuildingPart ids.
    fn add_feature(
        &mut self,
        feature: &CityJSONFeatureVertices,
        transform: &Transform,
        color_map: &HashMap<CityObjectType, [f32; 4]>,
    ) -> Result<()> {
        let building_id = feature
            .id
            .clone()
            .or_else(|| {
                feature
                    .cityobjects
                    .iter()
                    .find(|(_, co)| co.cotype == CityObjectType::Building)
                    .map(|(k, _)| k.clone())
            })
            .or_else(|| {
                feature
                    .cityobjects
                    .iter()
                    .find(|(_, co)| co.cotype == CityObjectType::BuildingPart)
                    .and_then(|(k, _)| k.strip_suffix("-0").map(|s| s.to_string()))
            })
            .unwrap_or_else(|| "unknown".to_string());

        let batch_index = self.next_batch_index;
        self.next_batch_index += 1;
        self.batch_id_to_cityobject_id.push(building_id);

        // Resolve primary CityObjectType: Building > BuildingPart > first type
        let primary_type = feature
            .cityobjects
            .values()
            .find(|co| co.cotype == CityObjectType::Building)
            .or_else(|| {
                feature
                    .cityobjects
                    .values()
                    .find(|co| co.cotype == CityObjectType::BuildingPart)
            })
            .unwrap_or_else(|| feature.cityobjects.values().next().unwrap());
        self.batch_id_to_cityobject_type
            .push(format!("{:?}", primary_type.cotype));

        let default_color = [1.0_f32, 0.753, 0.796, 1.0]; // #FFC0CB pink
        let mut vertex_cache: HashMap<usize, u32> = HashMap::new();
        for (_, co) in feature.cityobjects.iter() {
            let color = color_map.get(&co.cotype).copied().unwrap_or(default_color);
            if let Some(geoms) = &co.geometry {
                for geometry in geoms {
                    match geometry {
                        Geometry::MultiSurface { boundaries } => {
                            for surface in boundaries {
                                self.add_surface(surface, &feature.vertices, transform, &mut vertex_cache, batch_index, color)?;
                            }
                        }
                        Geometry::Solid { boundaries } => {
                            for shell in boundaries {
                                for surface in shell {
                                    self.add_surface(surface, &feature.vertices, transform, &mut vertex_cache, batch_index, color)?;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn add_surface(
        &mut self,
        surface: &Vec<Vec<usize>>,
        vertices_qc: &[[i64; 3]],
        transform: &Transform,
        cache: &mut HashMap<usize, u32>,
        batch_index: u32,
        color: [f32; 4],
    ) -> Result<()> {
        if surface.is_empty() {
            return Ok(());
        }
        let exterior = &surface[0];
        if exterior.len() < 3 {
            return Ok(());
        }

        let mut local_positions: Vec<[f32; 3]> = Vec::new();
        let mut glb_indices: Vec<u32> = Vec::new();
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
                let position = self.compute_local_position(vertex_id, vertices_qc, transform)?;
                let glb_index = self.vertex_index(vertex_id, position, cache, batch_index, color);
                local_positions.push(position);
                glb_indices.push(glb_index);
                vertex_count += 1;
            }
        }

        if glb_indices.len() < 3 {
            return Ok(());
        }

        // Fast path: single ring with 3 or 4 vertices (triangles and quads, e.g. tree
        // crown sides and trunk sides). Triangulate in 3D so vertical faces are not
        // collapsed by the drop-axis projection.
        if surface.len() == 1 && hole_indices.is_empty() {
            let n = exterior.len();
            if n == 3 {
                self.emit_triangles(vec![glb_indices[0], glb_indices[1], glb_indices[2]]);
                return Ok(());
            }
            if n == 4 {
                self.emit_triangles(vec![
                    glb_indices[0],
                    glb_indices[1],
                    glb_indices[2],
                    glb_indices[0],
                    glb_indices[2],
                    glb_indices[3],
                ]);
                return Ok(());
            }
        }

        let mut min = [f32::MAX; 3];
        let mut max = [f32::MIN; 3];
        for pos in &local_positions {
            for axis in 0..3 {
                if pos[axis] < min[axis] {
                    min[axis] = pos[axis];
                }
                if pos[axis] > max[axis] {
                    max[axis] = pos[axis];
                }
            }
        }

        let mut ranges = [0.0f32; 3];
        for axis in 0..3 {
            ranges[axis] = max[axis] - min[axis];
        }
        let drop_axis = ranges
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(idx, _)| idx)
            .unwrap_or(2);

        let mut flat_coords: Vec<f64> = Vec::with_capacity(local_positions.len() * 2);
        for pos in &local_positions {
            match drop_axis {
                0 => {
                    flat_coords.push(pos[1] as f64);
                    flat_coords.push(pos[2] as f64);
                }
                1 => {
                    flat_coords.push(pos[0] as f64);
                    flat_coords.push(pos[2] as f64);
                }
                _ => {
                    flat_coords.push(pos[0] as f64);
                    flat_coords.push(pos[1] as f64);
                }
            }
        }

        let triangulated = earcut(&flat_coords, &hole_indices, 2);
        if triangulated.len() < 3 {
            return Ok(());
        }

        let mut face_indices = Vec::with_capacity(triangulated.len());
        for idx in triangulated {
            face_indices.push(glb_indices[idx]);
        }

        self.emit_triangles(face_indices);
        Ok(())
    }

    fn vertex_index(
        &mut self,
        idx: usize,
        position: [f32; 3],
        cache: &mut HashMap<usize, u32>,
        batch_index: u32,
        color: [f32; 4],
    ) -> u32 {
        if let Some(&existing) = cache.get(&idx) {
            return existing;
        }
        self.positions.push(position);
        self.normals.push([0.0, 0.0, 0.0]);
        self.colors.push(color);
        self.batch_ids.push(batch_index);
        let index = (self.positions.len() - 1) as u32;
        cache.insert(idx, index);
        index
    }

    fn emit_triangles(&mut self, face_indices: Vec<u32>) {
        for tri in face_indices.chunks_exact(3) {
            let i0 = tri[0] as usize;
            let i1 = tri[1] as usize;
            let i2 = tri[2] as usize;

            let v0 = self.positions[i0];
            let v1 = self.positions[i1];
            let v2 = self.positions[i2];

            let u = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let v = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let normal = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];

            for &i in tri {
                let n = &mut self.normals[i as usize];
                n[0] += normal[0];
                n[1] += normal[1];
                n[2] += normal[2];
            }

            self.indices.extend_from_slice(tri);
        }
    }

    fn compute_local_position(
        &self,
        idx: usize,
        vertices_qc: &[[i64; 3]],
        transform: &Transform,
    ) -> Result<[f32; 3], anyhow::Error> {
        let [x_qc, y_qc, z_qc] = vertices_qc[idx];
        // 1. Dequantize coordinates to input CRS
        let x_input = (x_qc as f64 * transform.scale[0]) + transform.translate[0];
        let y_input = (y_qc as f64 * transform.scale[1]) + transform.translate[1];
        let z_input = (z_qc as f64 * transform.scale[2]) + transform.translate[2] + self.vertical_bias;

        // 2. Transform coordinates from input CRS to ECEF (EPSG:4978)
        // This ensures coordinate system consistency with root transform (which is in ECEF)
        let (x_ecef, y_ecef, z_ecef) = self.transformer_to_ecef
            .convert((x_input, y_input, z_input))
            .context("Transform vertex coordinates to ECEF")?;

        // 3. Make coordinates relative to root center in ECEF
        // Root transform will translate these relative coordinates to the correct ECEF position
        let x_local = (x_ecef - self.root_center_ecef.0) as f32;
        let y_local = (y_ecef - self.root_center_ecef.1) as f32;
        let z_local = (z_ecef - self.root_center_ecef.2) as f32;
        

        // Return ECEF coordinates relative to root center
        // Y-up transformation in glTF node will convert from ECEF Z-up to glTF Y-up standard
        Ok([x_local, y_local, z_local])
    }

    fn write_glb<P: AsRef<Path>>(&mut self, output_path: P, tiles_version: TilesVersion) -> Result<()> {
        self.normalize_normals();

        if self.positions.is_empty() {
            // Create parent directories if they don't exist
            if let Some(parent) = output_path.as_ref().parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("Failed to create parent directory for {:?}", output_path.as_ref()))?;
            }
            File::create(output_path.as_ref()).context("Create empty tile file")?;
            return Ok(());
        }

        let mut bin_buffer: Vec<u8> = Vec::new();

        let positions_offset = 0;
        for p in &self.positions {
            for component in p {
                bin_buffer.extend_from_slice(&component.to_le_bytes());
            }
        }

        let normals_offset = bin_buffer.len();
        for n in &self.normals {
            for component in n {
                bin_buffer.extend_from_slice(&component.to_le_bytes());
            }
        }

        // Per-vertex RGBA colors (COLOR_0)
        let colors_offset = bin_buffer.len();
        for c in &self.colors {
            for component in c {
                bin_buffer.extend_from_slice(&component.to_le_bytes());
            }
        }

        // Use U16 for batch/feature IDs so validators accept it (glTF 2.0 mesh attributes cannot use UNSIGNED_INT).
        let feature_count = self.next_batch_index as usize;
        if feature_count > 65535 {
            bail!(
                "Tile has {} features; glTF mesh attribute uses UNSIGNED_SHORT (max 65535)",
                feature_count
            );
        }
        let batch_ids_offset = bin_buffer.len();
        for &bid in &self.batch_ids {
            bin_buffer.extend_from_slice(&(bid as u16).to_le_bytes());
        }

        // Align to 4 bytes before u32 indices (batch_ids are u16, so when vertex
        // count is odd the buffer length is not a multiple of 4).
        let pad = (4 - (bin_buffer.len() % 4)) % 4;
        bin_buffer.extend(std::iter::repeat(0u8).take(pad));

        let indices_offset = bin_buffer.len();
        for index in &self.indices {
            bin_buffer.extend_from_slice(&index.to_le_bytes());
        }

        // EXT_structural_metadata string property tables (1.1 only)
        let mut string_data_offset = 0;
        let mut string_offsets_offset = 0;
        let mut string_offsets_len = 0;
        let mut id_string_data_len = 0usize;
        let mut type_string_data_offset = 0;
        let mut type_string_offsets_offset = 0;
        let mut type_string_offsets_len = 0;
        let mut type_string_data_len = 0usize;

        if tiles_version == TilesVersion::V1_1 {
            string_data_offset = bin_buffer.len();
            let mut string_offsets: Vec<u32> = Vec::with_capacity(self.batch_id_to_cityobject_id.len() + 1);
            let mut offset_acc: u32 = 0;
            for id in &self.batch_id_to_cityobject_id {
                string_offsets.push(offset_acc);
                let bytes = id.as_bytes();
                bin_buffer.extend_from_slice(bytes);
                bin_buffer.push(0); // null terminator
                offset_acc += (bytes.len() + 1) as u32;
            }
            string_offsets.push(offset_acc);
            id_string_data_len = bin_buffer.len() - string_data_offset;
            // Align to 4 bytes before u32 string offset array
            let pad = (4 - (bin_buffer.len() % 4)) % 4;
            bin_buffer.extend(std::iter::repeat(0u8).take(pad));
            string_offsets_offset = bin_buffer.len();
            string_offsets_len = string_offsets.len();
            for &o in &string_offsets {
                bin_buffer.extend_from_slice(&o.to_le_bytes());
            }

            type_string_data_offset = bin_buffer.len();
            let mut type_string_offsets: Vec<u32> =
                Vec::with_capacity(self.batch_id_to_cityobject_type.len() + 1);
            let mut type_offset_acc: u32 = 0;
            for type_name in &self.batch_id_to_cityobject_type {
                type_string_offsets.push(type_offset_acc);
                let bytes = type_name.as_bytes();
                bin_buffer.extend_from_slice(bytes);
                bin_buffer.push(0); // null terminator
                type_offset_acc += (bytes.len() + 1) as u32;
            }
            type_string_offsets.push(type_offset_acc);
            type_string_data_len = bin_buffer.len() - type_string_data_offset;
            // Align to 4 bytes before u32 type string offset array
            let pad = (4 - (bin_buffer.len() % 4)) % 4;
            bin_buffer.extend(std::iter::repeat(0u8).take(pad));
            type_string_offsets_offset = bin_buffer.len();
            type_string_offsets_len = type_string_offsets.len();
            for &o in &type_string_offsets {
                bin_buffer.extend_from_slice(&o.to_le_bytes());
            }
        }

        let accessor_positions = json::Accessor {
            buffer_view: Some(json::Index::new(0)),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(self.positions.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::F32,
            )),
            normalized: false,
            min: Some(json::Value::Array(
                (0..3)
                    .map(|axis| {
                        let min = self.positions.iter().map(|v| v[axis]).fold(f32::INFINITY, f32::min);
                        json::Value::from(min)
                    })
                    .collect(),
            )),
            max: Some(json::Value::Array(
                (0..3)
                    .map(|axis| {
                        let max = self.positions.iter().map(|v| v[axis]).fold(f32::NEG_INFINITY, f32::max);
                        json::Value::from(max)
                    })
                    .collect(),
            )),
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
            sparse: None,
        };

        let accessor_normals = json::Accessor {
            buffer_view: Some(json::Index::new(1)),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(self.normals.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::F32,
            )),
            normalized: false,
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec3),
            extensions: Default::default(),
            extras: Default::default(),
            min: None,
            max: None,
            name: None,
            sparse: None,
        };

        let accessor_colors = json::Accessor {
            buffer_view: Some(json::Index::new(2)),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(self.colors.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::F32,
            )),
            normalized: false,
            type_: json::validation::Checked::Valid(json::accessor::Type::Vec4),
            extensions: Default::default(),
            extras: Default::default(),
            min: None,
            max: None,
            name: None,
            sparse: None,
        };

        let accessor_batch_ids = json::Accessor {
            buffer_view: Some(json::Index::new(3)),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(self.batch_ids.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U16,
            )),
            normalized: false,
            min: Some(json::Value::from(vec![*self.batch_ids.iter().min().unwrap_or(&0)])),
            max: Some(json::Value::from(vec![*self.batch_ids.iter().max().unwrap_or(&0)])),
            type_: json::validation::Checked::Valid(json::accessor::Type::Scalar),
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
            sparse: None,
        };

        let accessor_indices = json::Accessor {
            buffer_view: Some(json::Index::new(4)),
            byte_offset: Some(json::validation::USize64(0)),
            count: json::validation::USize64(self.indices.len() as u64),
            component_type: json::validation::Checked::Valid(json::accessor::GenericComponentType(
                json::accessor::ComponentType::U32,
            )),
            normalized: false,
            min: Some(json::Value::from(vec![0])),
            max: Some(json::Value::from(vec![self.indices.iter().copied().max().unwrap_or(0)])),
            type_: json::validation::Checked::Valid(json::accessor::Type::Scalar),
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
            sparse: None,
        };

        // BV 0-4: common to both 1.0 and 1.1
        let mut buffer_views = vec![
            // BV 0: POSITION
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.positions.len() * 12) as u64),
                byte_offset: Some(json::validation::USize64(positions_offset as u64)),
                byte_stride: None,
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
            // BV 1: NORMAL
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.normals.len() * 12) as u64),
                byte_offset: Some(json::validation::USize64(normals_offset as u64)),
                byte_stride: None,
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
            // BV 2: COLOR_0 (Vec4 F32 = 16 bytes/vertex)
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.colors.len() * 16) as u64),
                byte_offset: Some(json::validation::USize64(colors_offset as u64)),
                byte_stride: None,
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
            // BV 3: _BATCHID (1.0) or _FEATURE_ID_0 (1.1)
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.batch_ids.len() * 2) as u64),
                byte_offset: Some(json::validation::USize64(batch_ids_offset as u64)),
                byte_stride: None,
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
            // BV 4: indices
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.indices.len() * 4) as u64),
                byte_offset: Some(json::validation::USize64(indices_offset as u64)),
                byte_stride: None,
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ElementArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
        ];

        // BV 5-8: EXT_structural_metadata string property tables (1.1 only)
        if tiles_version == TilesVersion::V1_1 {
            buffer_views.extend([
                // BV 5: BuildingId string data
                json::buffer::View {
                    buffer: json::Index::new(0),
                    byte_length: json::validation::USize64(id_string_data_len as u64),
                    byte_offset: Some(json::validation::USize64(string_data_offset as u64)),
                    byte_stride: None,
                    target: None,
                    extensions: Default::default(),
                    extras: Default::default(),
                    name: None,
                },
                // BV 6: BuildingId string offsets
                json::buffer::View {
                    buffer: json::Index::new(0),
                    byte_length: json::validation::USize64((string_offsets_len * 4) as u64),
                    byte_offset: Some(json::validation::USize64(string_offsets_offset as u64)),
                    byte_stride: None,
                    target: None,
                    extensions: Default::default(),
                    extras: Default::default(),
                    name: None,
                },
                // BV 7: CityObjectType string data
                json::buffer::View {
                    buffer: json::Index::new(0),
                    byte_length: json::validation::USize64(type_string_data_len as u64),
                    byte_offset: Some(json::validation::USize64(type_string_data_offset as u64)),
                    byte_stride: None,
                    target: None,
                    extensions: Default::default(),
                    extras: Default::default(),
                    name: None,
                },
                // BV 8: CityObjectType string offsets
                json::buffer::View {
                    buffer: json::Index::new(0),
                    byte_length: json::validation::USize64((type_string_offsets_len * 4) as u64),
                    byte_offset: Some(json::validation::USize64(type_string_offsets_offset as u64)),
                    byte_stride: None,
                    target: None,
                    extensions: Default::default(),
                    extras: Default::default(),
                    name: None,
                },
            ]);
        }

        let mut attributes = std::collections::BTreeMap::new();
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Positions),
            json::Index::new(0),
        );
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Normals),
            json::Index::new(1),
        );
        // Per-vertex RGBA color
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Colors(0)),
            json::Index::new(2),
        );
        // 1.1: _FEATURE_ID_0 for EXT_mesh_features; 1.0: _BATCHID for batch table
        let batch_attr_name = match tiles_version {
            TilesVersion::V1_1 => "FEATURE_ID_0",
            TilesVersion::V1_0 => "BATCHID",
        };
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Extras(batch_attr_name.into())),
            json::Index::new(3),
        );

        // White base material — vertex colors provide the actual per-type coloring
        let material = create_material("#FFFFFF", self.metallic_factor, self.roughness_factor)?;

        let feature_count = self.next_batch_index;

        // EXT_mesh_features primitive extension (1.1 only)
        let primitive_extensions = match tiles_version {
            TilesVersion::V1_1 => {
                let mut ext_others = serde_json::Map::new();
                ext_others.insert(
                    "EXT_mesh_features".to_string(),
                    serde_json::json!({
                        "featureIds": [{
                            "attribute": 3,
                            "featureCount": feature_count,
                            "propertyTable": 0
                        }]
                    }),
                );
                Some(json::extensions::mesh::Primitive {
                    others: ext_others,
                    ..Default::default()
                })
            }
            TilesVersion::V1_0 => None,
        };

        let primitive = json::mesh::Primitive {
            attributes,
            indices: Some(json::Index::new(4)),
            material: Some(json::Index::new(0)),
            mode: json::validation::Checked::Valid(json::mesh::Mode::Triangles),
            targets: None,
            extensions: primitive_extensions,
            extras: Default::default(),
        };

        let mesh = json::Mesh {
            primitives: vec![primitive],
            weights: None,
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
        };

        // Apply Y-up transformation matrix to convert from ECEF (Z-up) to glTF standard (Y-up)
        // GLB content coordinates are in ECEF (Z-up), and this matrix converts them to glTF Y-up format
        // This matches pg2b3dm's approach - the Y-up matrix is needed in GLB node
        // Matrix format: [1,0,0,0, 0,0,-1,0, 0,1,0,0, 0,0,0,1] (column-major in glTF JSON)
        // Transformation: X'=X, Y'=-Z (ECEF Z becomes glTF -Y), Z'=Y (ECEF Y becomes glTF Z)
        let y_up_matrix = [
            1.0, 0.0, 0.0, 0.0,   // Column 0: [1, 0, 0, 0] - X axis
            0.0, 0.0, -1.0, 0.0,  // Column 1: [0, 0, -1, 0] - Y axis becomes -Z
            0.0, 1.0, 0.0, 0.0,   // Column 2: [0, 1, 0, 0] - Z axis becomes Y
            0.0, 0.0, 0.0, 1.0,   // Column 3: [0, 0, 0, 1] - Translation/scale
        ];

        let node = json::Node {
            mesh: Some(json::Index::new(0)),
            camera: None,
            children: None,
            skin: None,
            matrix: Some(y_up_matrix),
            rotation: None,
            scale: None,
            translation: None,
            weights: None,
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
        };

        let scene = json::Scene {
            nodes: vec![json::Index::new(0)],
            extensions: Default::default(),
            extras: Default::default(),
            name: None,
        };

        // EXT_structural_metadata + extensions_used (1.1 only)
        let (root_extensions, extensions_used) = match tiles_version {
            TilesVersion::V1_1 => {
                let n_features = self.batch_id_to_cityobject_id.len();
                let structural_metadata_ext = serde_json::json!({
                    "schema": {
                        "classes": {
                            "Feature": {
                                "properties": {
                                    "BuildingId": {
                                        "type": "STRING",
                                        "stringOffsetType": "UINT32"
                                    },
                                    "CityObjectType": {
                                        "type": "STRING",
                                        "stringOffsetType": "UINT32"
                                    }
                                }
                            }
                        }
                    },
                    "propertyTables": [{
                        "class": "Feature",
                        "count": n_features,
                        "properties": {
                            "BuildingId": {
                                "values": 5,
                                "stringOffsets": 6
                            },
                            "CityObjectType": {
                                "values": 7,
                                "stringOffsets": 8
                            }
                        }
                    }]
                });
                let mut root_ext_others = serde_json::Map::new();
                root_ext_others.insert("EXT_structural_metadata".to_string(), structural_metadata_ext);
                (
                    Some(json::extensions::root::Root {
                        others: root_ext_others,
                        ..Default::default()
                    }),
                    vec!["EXT_mesh_features".into(), "EXT_structural_metadata".into()],
                )
            }
            TilesVersion::V1_0 => (None, vec![]),
        };

        let root = json::Root {
            accessors: vec![accessor_positions, accessor_normals, accessor_colors, accessor_batch_ids, accessor_indices],
            extensions_used,
            extensions: root_extensions,
            buffers: vec![json::Buffer {
                byte_length: json::validation::USize64(bin_buffer.len() as u64),
                uri: None,
                name: Some("buffer0".into()),
                extensions: Default::default(),
                extras: Default::default(),
            }],
            buffer_views,
            materials: vec![material],
            meshes: vec![mesh],
            nodes: vec![node],
            scenes: vec![scene],
            scene: Some(json::Index::new(0)),
            asset: json::Asset {
                version: GLTF_VERSION.into(),
                generator: Some("tyler".into()),
                copyright: None,
                ..Default::default()
            },
            ..Default::default()
        };

        let mut json_bytes = json::serialize::to_string(&root)?.into_bytes();
        let json_padding = (4 - (json_bytes.len() % 4)) % 4;
        json_bytes.extend(std::iter::repeat(b' ').take(json_padding));

        let bin_padding = (4 - (bin_buffer.len() % 4)) % 4;
        bin_buffer.extend(std::iter::repeat(0).take(bin_padding));

        let total_length = 12 + 8 + json_bytes.len() + 8 + bin_buffer.len();
        let mut glb_bytes = Vec::with_capacity(total_length);
        glb_bytes.extend_from_slice(b"glTF");
        glb_bytes.extend_from_slice(&2u32.to_le_bytes());
        glb_bytes.extend_from_slice(&(total_length as u32).to_le_bytes());

        glb_bytes.extend_from_slice(&(json_bytes.len() as u32).to_le_bytes());
        glb_bytes.extend_from_slice(b"JSON");
        glb_bytes.extend_from_slice(&json_bytes);

        glb_bytes.extend_from_slice(&(bin_buffer.len() as u32).to_le_bytes());
        glb_bytes.extend_from_slice(b"BIN\0");
        glb_bytes.extend_from_slice(&bin_buffer);

        // For 1.0: wrap GLB in B3DM container with feature table + batch table
        let final_bytes = match tiles_version {
            TilesVersion::V1_1 => glb_bytes,
            TilesVersion::V1_0 => wrap_b3dm(
                &glb_bytes,
                &self.batch_id_to_cityobject_id,
                &self.batch_id_to_cityobject_type,
            ),
        };

        // Create parent directories if they don't exist
        if let Some(parent) = output_path.as_ref().parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create parent directory for {:?}", output_path.as_ref()))?;
        }

        let mut file = File::create(output_path)?;
        file.write_all(&final_bytes)?;

        Ok(())
    }

    fn normalize_normals(&mut self) {
        for normal in self.normals.iter_mut() {
            let length = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
            if length > f32::EPSILON {
                normal[0] /= length;
                normal[1] /= length;
                normal[2] /= length;
            } else {
                // Zero-length normal (degenerate triangles). Use default up vector [0, 1, 0]
                // which is appropriate for glTF Y-up coordinate system
                normal[0] = 0.0;
                normal[1] = 1.0;
                normal[2] = 0.0;
            }
        }
    }
}

/// Wrap a GLB binary in a B3DM container with feature table and batch table.
/// B3DM format: 28-byte header + feature table JSON + batch table JSON + GLB body.
fn wrap_b3dm(glb_bytes: &[u8], building_ids: &[String], city_object_types: &[String]) -> Vec<u8> {
    let batch_length = building_ids.len();

    // Feature table JSON: {"BATCH_LENGTH": n}
    let mut ft_json_bytes = serde_json::to_vec(&serde_json::json!({ "BATCH_LENGTH": batch_length })).unwrap();
    // Pad feature table JSON to 8-byte alignment (header is 28 bytes)
    let ft_padding = (8 - ((28 + ft_json_bytes.len()) % 8)) % 8;
    ft_json_bytes.extend(std::iter::repeat(b' ').take(ft_padding));

    // Batch table JSON: {"BuildingId": [...], "CityObjectType": [...]}
    let mut bt_json_bytes = serde_json::to_vec(&serde_json::json!({
        "BuildingId": building_ids,
        "CityObjectType": city_object_types,
    })).unwrap();
    // Pad batch table JSON to 8-byte alignment
    let bt_padding = (8 - ((28 + ft_json_bytes.len() + bt_json_bytes.len()) % 8)) % 8;
    bt_json_bytes.extend(std::iter::repeat(b' ').take(bt_padding));

    let total = 28 + ft_json_bytes.len() + bt_json_bytes.len() + glb_bytes.len();
    let mut b3dm = Vec::with_capacity(total);

    // 28-byte header
    b3dm.extend_from_slice(b"b3dm");                                       // magic
    b3dm.extend_from_slice(&1u32.to_le_bytes());                           // version
    b3dm.extend_from_slice(&(total as u32).to_le_bytes());                 // byteLength
    b3dm.extend_from_slice(&(ft_json_bytes.len() as u32).to_le_bytes());   // featureTableJSONByteLength
    b3dm.extend_from_slice(&0u32.to_le_bytes());                           // featureTableBinaryByteLength
    b3dm.extend_from_slice(&(bt_json_bytes.len() as u32).to_le_bytes());   // batchTableJSONByteLength
    b3dm.extend_from_slice(&0u32.to_le_bytes());                           // batchTableBinaryByteLength

    // Body
    b3dm.extend_from_slice(&ft_json_bytes);
    b3dm.extend_from_slice(&bt_json_bytes);
    b3dm.extend_from_slice(glb_bytes);

    b3dm
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a MeshBuilder with synthetic data and write to a GLB file.
    /// Returns (raw_glb_bytes, raw_file_bytes) — for V1_0, raw_file_bytes is
    /// the B3DM wrapper and raw_glb_bytes is the embedded GLB extracted from it.
    /// For V1_1, both are identical.
    fn build_test_glb(
        n_vertices: usize,
        n_features: usize,
        tiles_version: TilesVersion,
    ) -> (Vec<u8>, Vec<u8>) {
        assert!(n_vertices % 3 == 0, "n_vertices must be a multiple of 3");
        assert!(n_features >= 1, "need at least 1 feature");

        let transformer = Proj::new_known_crs("EPSG:7415", "EPSG:4978", None)
            .expect("PROJ transformer");
        let mut builder = MeshBuilder::new(transformer, (0.0, 0.0, 0.0), 0.0, 0.0, 1.0);

        // Populate synthetic vertex data
        for i in 0..n_vertices {
            let f = i as f32;
            builder.positions.push([f, f + 1.0, f + 2.0]);
            builder.normals.push([0.0, 1.0, 0.0]);
            builder.colors.push([1.0, 0.0, 0.0, 1.0]);
            builder.batch_ids.push((i % n_features) as u32);
        }
        for i in 0..n_vertices as u32 {
            builder.indices.push(i);
        }
        builder.next_batch_index = n_features as u32;

        // Varying string lengths to exercise alignment edge cases
        let id_patterns = ["A", "AB", "ABC"];
        let type_patterns = ["Building", "WaterBody", "Bridge"];
        for i in 0..n_features {
            builder
                .batch_id_to_cityobject_id
                .push(id_patterns[i % id_patterns.len()].to_string());
            builder
                .batch_id_to_cityobject_type
                .push(type_patterns[i % type_patterns.len()].to_string());
        }

        let tmp = tempfile::NamedTempFile::new().expect("tempfile");
        let path = tmp.path().to_path_buf();
        builder
            .write_glb(&path, tiles_version)
            .expect("write_glb failed");

        let file_bytes = std::fs::read(&path).expect("read GLB file");

        let glb_bytes = match tiles_version {
            TilesVersion::V1_1 => file_bytes.clone(),
            TilesVersion::V1_0 => {
                // Extract embedded GLB from B3DM: skip 28-byte header + feature/batch table JSON
                let ft_json_len =
                    u32::from_le_bytes(file_bytes[12..16].try_into().unwrap()) as usize;
                let bt_json_len =
                    u32::from_le_bytes(file_bytes[20..24].try_into().unwrap()) as usize;
                let glb_start = 28 + ft_json_len + bt_json_len;
                file_bytes[glb_start..].to_vec()
            }
        };

        (glb_bytes, file_bytes)
    }

    #[test]
    fn glb_no_invalid_byte_stride() {
        for &n in &[3, 5, 6, 7, 100, 101] {
            for &version in &[TilesVersion::V1_0, TilesVersion::V1_1] {
                let (glb_bytes, _) = build_test_glb(n, 2, version);
                let gltf = gltf::Gltf::from_slice(&glb_bytes)
                    .unwrap_or_else(|e| panic!("Failed to parse GLB (n={n}, v={version:?}): {e}"));

                for bv in gltf.document.views() {
                    if let Some(stride) = bv.stride() {
                        assert!(
                            stride >= 4,
                            "byteStride {} < 4 on BV {} (n={n}, v={version:?})",
                            stride,
                            bv.index()
                        );
                        assert!(
                            stride % 4 == 0,
                            "byteStride {} not multiple of 4 on BV {} (n={n}, v={version:?})",
                            stride,
                            bv.index()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn glb_buffer_view_offsets_aligned() {
        for &n in &[3, 5, 7, 9, 101] {
            for &version in &[TilesVersion::V1_0, TilesVersion::V1_1] {
                let (glb_bytes, _) = build_test_glb(n, 2, version);
                let gltf = gltf::Gltf::from_slice(&glb_bytes)
                    .unwrap_or_else(|e| panic!("Failed to parse GLB (n={n}, v={version:?}): {e}"));

                for accessor in gltf.document.accessors() {
                    let bv = accessor.view().unwrap();
                    let comp_size = match accessor.data_type() {
                        gltf::accessor::DataType::U8 | gltf::accessor::DataType::I8 => 1,
                        gltf::accessor::DataType::U16 | gltf::accessor::DataType::I16 => 2,
                        gltf::accessor::DataType::U32 | gltf::accessor::DataType::F32 => 4,
                    };
                    let total_offset = accessor.offset() + bv.offset();
                    assert!(
                        total_offset % comp_size == 0,
                        "Accessor {} (BV {}): offset {} not aligned to component size {} (n={n}, v={version:?})",
                        accessor.index(),
                        bv.index(),
                        total_offset,
                        comp_size
                    );
                }
            }
        }
    }

    #[test]
    fn glb_buffer_view_bounds_valid() {
        for &n in &[3, 6, 99, 100] {
            for &version in &[TilesVersion::V1_0, TilesVersion::V1_1] {
                let n_features = if n >= 6 { 3 } else { 1 };
                let (glb_bytes, _) = build_test_glb(n, n_features, version);
                let gltf = gltf::Gltf::from_slice(&glb_bytes)
                    .unwrap_or_else(|e| panic!("Failed to parse GLB (n={n}, v={version:?}): {e}"));

                let buffer_len = gltf.document.buffers().next().unwrap().length();
                for bv in gltf.document.views() {
                    assert!(
                        bv.offset() + bv.length() <= buffer_len,
                        "BV {} exceeds buffer: offset {} + length {} > {} (n={n}, v={version:?})",
                        bv.index(),
                        bv.offset(),
                        bv.length(),
                        buffer_len
                    );
                }
            }
        }
    }

    #[test]
    fn glb_valid_structure() {
        for &version in &[TilesVersion::V1_0, TilesVersion::V1_1] {
            let (glb_bytes, _) = build_test_glb(6, 2, version);
            let gltf = gltf::Gltf::from_slice(&glb_bytes)
                .unwrap_or_else(|e| panic!("Failed to parse GLB (v={version:?}): {e}"));

            assert_eq!(gltf.document.meshes().count(), 1);
            assert_eq!(gltf.document.nodes().count(), 1);
            assert_eq!(gltf.document.scenes().count(), 1);
            assert_eq!(gltf.document.buffers().count(), 1);
            assert_eq!(gltf.document.accessors().count(), 5);

            let expected_bvs = match version {
                TilesVersion::V1_1 => 9,
                TilesVersion::V1_0 => 5,
            };
            assert_eq!(
                gltf.document.views().count(),
                expected_bvs,
                "Wrong BV count for {version:?}"
            );

            // Verify asset info by parsing JSON chunk directly
            let json_chunk_len =
                u32::from_le_bytes(glb_bytes[12..16].try_into().unwrap()) as usize;
            let json_bytes = &glb_bytes[20..20 + json_chunk_len];
            let root: serde_json::Value = serde_json::from_slice(json_bytes).unwrap();
            assert_eq!(root["asset"]["version"], "2.0");
            assert_eq!(root["asset"]["generator"], "tyler");

            // Verify accessor types and counts
            let accessors: Vec<_> = gltf.document.accessors().collect();
            // Accessor 0: POSITION (Vec3/F32, count=6)
            assert_eq!(accessors[0].count(), 6);
            assert_eq!(accessors[0].data_type(), gltf::accessor::DataType::F32);
            assert_eq!(accessors[0].dimensions(), gltf::accessor::Dimensions::Vec3);
            // Accessor 1: NORMAL (Vec3/F32)
            assert_eq!(accessors[1].data_type(), gltf::accessor::DataType::F32);
            assert_eq!(accessors[1].dimensions(), gltf::accessor::Dimensions::Vec3);
            // Accessor 2: COLOR_0 (Vec4/F32)
            assert_eq!(accessors[2].data_type(), gltf::accessor::DataType::F32);
            assert_eq!(accessors[2].dimensions(), gltf::accessor::Dimensions::Vec4);
            // Accessor 3: BATCHID/FEATURE_ID (Scalar/U16)
            assert_eq!(accessors[3].data_type(), gltf::accessor::DataType::U16);
            assert_eq!(
                accessors[3].dimensions(),
                gltf::accessor::Dimensions::Scalar
            );
            // Accessor 4: indices (Scalar/U32)
            assert_eq!(accessors[4].data_type(), gltf::accessor::DataType::U32);
            assert_eq!(
                accessors[4].dimensions(),
                gltf::accessor::Dimensions::Scalar
            );
            assert_eq!(accessors[4].count(), 6);
        }
    }

    #[test]
    fn glb_chunk_alignment() {
        for &n in &[3, 5] {
            let (glb_bytes, _) = build_test_glb(n, 1, TilesVersion::V1_1);

            // GLB header
            assert_eq!(&glb_bytes[0..4], b"glTF");
            let version = u32::from_le_bytes(glb_bytes[4..8].try_into().unwrap());
            assert_eq!(version, 2);
            let total_len = u32::from_le_bytes(glb_bytes[8..12].try_into().unwrap()) as usize;
            assert_eq!(total_len, glb_bytes.len());

            // JSON chunk alignment
            let json_chunk_len =
                u32::from_le_bytes(glb_bytes[12..16].try_into().unwrap()) as usize;
            assert!(
                json_chunk_len % 4 == 0,
                "JSON chunk length {} not 4-byte aligned (n={n})",
                json_chunk_len
            );

            // BIN chunk alignment
            let bin_chunk_start = 12 + 8 + json_chunk_len;
            let bin_chunk_len =
                u32::from_le_bytes(glb_bytes[bin_chunk_start..bin_chunk_start + 4].try_into().unwrap())
                    as usize;
            assert!(
                bin_chunk_len % 4 == 0,
                "BIN chunk length {} not 4-byte aligned (n={n})",
                bin_chunk_len
            );
        }
    }

    #[test]
    fn b3dm_wrapping_preserves_valid_glb() {
        let (glb_bytes, file_bytes) = build_test_glb(6, 2, TilesVersion::V1_0);

        // Verify B3DM header
        assert_eq!(&file_bytes[0..4], b"b3dm");
        let b3dm_version = u32::from_le_bytes(file_bytes[4..8].try_into().unwrap());
        assert_eq!(b3dm_version, 1);
        let b3dm_total = u32::from_le_bytes(file_bytes[8..12].try_into().unwrap()) as usize;
        assert_eq!(b3dm_total, file_bytes.len());

        // Verify the embedded GLB is valid and passes alignment checks
        let gltf = gltf::Gltf::from_slice(&glb_bytes)
            .expect("B3DM embedded GLB should be valid glTF 2.0");

        for bv in gltf.document.views() {
            if let Some(stride) = bv.stride() {
                assert!(stride >= 4, "byteStride {} < 4 in B3DM GLB", stride);
                assert!(stride % 4 == 0, "byteStride {} not multiple of 4 in B3DM GLB", stride);
            }
        }

        for accessor in gltf.document.accessors() {
            let bv = accessor.view().unwrap();
            let comp_size = match accessor.data_type() {
                gltf::accessor::DataType::U8 | gltf::accessor::DataType::I8 => 1,
                gltf::accessor::DataType::U16 | gltf::accessor::DataType::I16 => 2,
                gltf::accessor::DataType::U32 | gltf::accessor::DataType::F32 => 4,
            };
            let total_offset = accessor.offset() + bv.offset();
            assert!(
                total_offset % comp_size == 0,
                "B3DM GLB: accessor {} offset {} not aligned to {}",
                accessor.index(),
                total_offset,
                comp_size
            );
        }
    }
}
