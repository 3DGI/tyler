use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};
use earcutr::earcut;
use gltf::json as json;

use crate::parser::{CityJSONFeatureVertices, Geometry, Transform, World};
use crate::proj::Proj;
use crate::spatial_structs::{QuadTree, QuadTreeNodeId};

const GLTF_VERSION: &str = "2.0";

/// Parse hex color string (#RRGGBB) to RGBA f32 array [R, G, B, A]
fn hex_to_rgba(hex: &str) -> Result<[f32; 4], anyhow::Error> {
    if hex.len() != 7 || !hex.starts_with('#') {
        bail!("Invalid hex color format: expected #RRGGBB");
    }
    let hex_digits = &hex[1..];
    let r = u8::from_str_radix(&hex_digits[0..2], 16)?;
    let g = u8::from_str_radix(&hex_digits[2..4], 16)?;
    let b = u8::from_str_radix(&hex_digits[4..6], 16)?;
    Ok([
        r as f32 / 255.0,
        g as f32 / 255.0,
        b as f32 / 255.0,
        1.0, // Alpha, always opaque
    ])
}

/// Create default PBR material matching pg2b3dm's structure
fn create_default_material(base_color: &str) -> Result<json::Material, anyhow::Error> {
    let base_color_rgba = hex_to_rgba(base_color)?;
    
    // Metallic roughness factor from pg2b3dm default: #008000
    // Green channel = 128/255 = 0.501960... (roughness)
    // Red channel = 0/255 = 0.0 (metallic)
    let roughness_factor = 128.0 / 255.0;
    let metallic_factor = 0.0;
    
    Ok(json::Material {
        name: None,
        extensions: Default::default(),
        extras: Default::default(),
        pbr_metallic_roughness: json::material::PbrMetallicRoughness {
            base_color_factor: json::material::PbrBaseColorFactor(base_color_rgba),
            metallic_factor: json::material::StrengthFactor(metallic_factor),
            roughness_factor: json::material::StrengthFactor(roughness_factor),
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

/// Write a GLB file for a single tile.
/// 
/// # Arguments
/// * `world` - The world containing features and grid
/// * `quadtree` - The quadtree structure
/// * `qtree_node_id` - The node ID of the tile to write
/// * `output_path` - Where to write the GLB file
/// * `default_color` - Default color for the mesh
pub fn write_tile_glb<P: AsRef<Path>>(
    world: &World,
    quadtree: &QuadTree,
    qtree_node_id: QuadTreeNodeId,
    output_path: P,
    default_color: &str,
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
    
    // Optimized path: estimate buffer sizes based on tile contents (Priority 3)
    // Estimate: avg ~50 vertices per feature, ~100 triangles (300 indices)
    let estimated_features: usize = qtree_node.cells().iter()
        .map(|cellid| world.grid.cell(cellid).feature_ids.len())
        .sum();
    let estimated_vertices = estimated_features * 50;
    let estimated_indices = estimated_features * 300;
    
    let mut builder = MeshBuilder::new_with_capacity(
        transformer_to_ecef,
        root_center_ecef,
        vertical_geoid_n,
        estimated_vertices,
        estimated_indices,
    );

    for cellid in qtree_node.cells() {
        let cell = world.grid.cell(cellid);
        for fid in cell.feature_ids.iter() {
            let feature = &world.features[*fid];
            // Check FCB cache first (for FCB-sourced features), otherwise read from file (JSONL path)
            if let Some(ref cache) = world.fcb_feature_cache {
                // FCB path - use cached feature (reference, no clone)
                let cf = cache.get(&feature.path_jsonl)
                    .ok_or_else(|| anyhow::anyhow!("Feature not found in FCB cache: {:?}", feature.path_jsonl))?;
                
                // Optimized path: batch PROJ transformation (Priority 2)
                builder.add_feature_batch(cf, &world.transform)?;
            } else {
                // JSONL path - read from file (existing behavior)
                let cf = CityJSONFeatureVertices::from_file(&feature.path_jsonl)
                    .map_err(|e| anyhow::anyhow!("Failed to read {:?}: {}", feature.path_jsonl, e))?;
                
                // Optimized path: batch PROJ transformation (Priority 2)
                builder.add_feature_batch(&cf, &world.transform)?;
            }
        }
    }

    builder.write_glb(output_path, default_color)
}

struct MeshBuilder {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    indices: Vec<u32>,
    transformer_to_ecef: Proj,
    root_center_ecef: (f64, f64, f64),
    vertical_bias: f64,
    // Reusable buffers for batch operations (tip 9 - reuse buffers)
    batch_input_buffer: Vec<(f64, f64, f64)>,
    batch_output_buffer: Vec<(f64, f64, f64)>,
}

impl MeshBuilder {
    /// Create MeshBuilder with pre-allocated capacity (Priority 3 optimization)
    /// 
    /// # Arguments
    /// * `estimated_vertices` - Expected number of vertices (positions and normals)
    /// * `estimated_indices` - Expected number of indices (triangles * 3)
    #[allow(dead_code)]
    fn new_with_capacity(
        transformer_to_ecef: Proj,
        root_center_ecef: (f64, f64, f64),
        vertical_bias: f64,
        estimated_vertices: usize,
        estimated_indices: usize,
    ) -> Self {
        Self {
            positions: Vec::with_capacity(estimated_vertices),
            normals: Vec::with_capacity(estimated_vertices),
            indices: Vec::with_capacity(estimated_indices),
            transformer_to_ecef,
            root_center_ecef,
            vertical_bias,
            // Pre-allocate batch buffers for typical feature vertex count
            batch_input_buffer: Vec::with_capacity(256),
            batch_output_buffer: Vec::with_capacity(256),
        }
    }

    /// Add feature with batch PROJ transformation (Priority 2 optimization)
    /// Steps:
    /// 1. Collect unique vertex indices used by all geometries
    /// 2. Dequantize all vertices in one pass
    /// 3. Batch transform through PROJ (5-15% gain from reduced FFI overhead)
    /// 4. Build surfaces using pre-transformed vertex cache
    fn add_feature_batch(&mut self, feature: &CityJSONFeatureVertices, transform: &Transform) -> Result<()> {
        // Step 1: Collect unique vertex indices (avoid processing unused vertices)
        let mut used_indices: Vec<usize> = Vec::with_capacity(feature.vertices.len());
        let mut index_set: std::collections::HashSet<usize> = std::collections::HashSet::with_capacity(feature.vertices.len());
        
        for (_id, co) in feature.cityobjects.iter() {
            if let Some(geoms) = &co.geometry {
                for geometry in geoms {
                    match geometry {
                        Geometry::MultiSurface { boundaries } => {
                            for surface in boundaries {
                                for ring in surface {
                                    for &idx in ring {
                                        if index_set.insert(idx) {
                                            used_indices.push(idx);
                                        }
                                    }
                                }
                            }
                        }
                        Geometry::Solid { boundaries } => {
                            for shell in boundaries {
                                for surface in shell {
                                    for ring in surface {
                                        for &idx in ring {
                                            if index_set.insert(idx) {
                                                used_indices.push(idx);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        drop(index_set); // Free memory early
        
        if used_indices.is_empty() {
            return Ok(());
        }

        // Step 2: Dequantize vertices into batch input buffer (reuse buffer - tip 9)
        self.batch_input_buffer.clear();
        self.batch_input_buffer.reserve(used_indices.len());
        
        for &idx in &used_indices {
            let [x_qc, y_qc, z_qc] = feature.vertices[idx];
            let x_input = (x_qc as f64 * transform.scale[0]) + transform.translate[0];
            let y_input = (y_qc as f64 * transform.scale[1]) + transform.translate[1];
            let z_input = (z_qc as f64 * transform.scale[2]) + transform.translate[2] + self.vertical_bias;
            self.batch_input_buffer.push((x_input, y_input, z_input));
        }

        // Step 3: Batch transform through PROJ (Priority 2 - reduces FFI call overhead)
        self.transformer_to_ecef
            .convert_batch_into(&self.batch_input_buffer, &mut self.batch_output_buffer)
            .context("Batch transform vertices to ECEF")?;

        // Step 4: Build transformed vertex cache (vertex_idx -> local position)
        // Pre-allocate with exact capacity (tip 1)
        let mut transformed_cache: HashMap<usize, [f32; 3]> = HashMap::with_capacity(used_indices.len());
        
        for (i, &idx) in used_indices.iter().enumerate() {
            let (x_ecef, y_ecef, z_ecef) = self.batch_output_buffer[i];
            let local_pos = [
                (x_ecef - self.root_center_ecef.0) as f32,
                (y_ecef - self.root_center_ecef.1) as f32,
                (z_ecef - self.root_center_ecef.2) as f32,
            ];
            transformed_cache.insert(idx, local_pos);
        }

        // Step 5: Build surfaces using pre-transformed vertices
        let mut vertex_cache: HashMap<usize, u32> = HashMap::with_capacity(used_indices.len());
        
        for (_id, co) in feature.cityobjects.iter() {
            if let Some(geoms) = &co.geometry {
                for geometry in geoms {
                    match geometry {
                        Geometry::MultiSurface { boundaries } => {
                            for surface in boundaries {
                                self.add_surface_with_cache(surface, &transformed_cache, &mut vertex_cache)?;
                            }
                        }
                        Geometry::Solid { boundaries } => {
                            for shell in boundaries {
                                for surface in shell {
                                    self.add_surface_with_cache(surface, &transformed_cache, &mut vertex_cache)?;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Add surface using pre-transformed vertex cache (for batch optimization)
    fn add_surface_with_cache(
        &mut self,
        surface: &[Vec<usize>],
        transformed_cache: &HashMap<usize, [f32; 3]>,
        vertex_cache: &mut HashMap<usize, u32>,
    ) -> Result<()> {
        if surface.is_empty() {
            return Ok(());
        }
        let exterior = &surface[0];
        if exterior.len() < 3 {
            return Ok(());
        }

        // Pre-allocate local buffers with estimated capacity (tip 1)
        let estimated_verts = surface.iter().map(|r| r.len()).sum::<usize>();
        let mut local_positions: Vec<[f32; 3]> = Vec::with_capacity(estimated_verts);
        let mut glb_indices: Vec<u32> = Vec::with_capacity(estimated_verts);
        let mut hole_indices: Vec<usize> = Vec::with_capacity(surface.len().saturating_sub(1));
        let mut vertex_count = 0usize;

        for (ring_idx, ring) in surface.iter().enumerate() {
            if ring.len() < 3 {
                continue;
            }
            if ring_idx > 0 {
                hole_indices.push(vertex_count);
            }
            for &idx in ring {
                // Get pre-transformed position from cache
                let pos = transformed_cache.get(&idx)
                    .ok_or_else(|| anyhow::anyhow!("Vertex index {} not found in transformed cache", idx))?;
                
                // Check if this original vertex index is already in the GLB mesh
                let glb_idx = if let Some(&existing_idx) = vertex_cache.get(&idx) {
                    existing_idx
                } else {
                    let new_idx = self.positions.len() as u32;
                    self.positions.push(*pos);
                    self.normals.push([0.0, 0.0, 0.0]);
                    vertex_cache.insert(idx, new_idx);
                    new_idx
                };

                local_positions.push(*pos);
                glb_indices.push(glb_idx);
                vertex_count += 1;
            }
        }

        if glb_indices.len() < 3 {
            return Ok(());
        }

        // Determine which axis to drop for 2D triangulation (same logic as add_surface)
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

        // Build flat_coords for earcut, dropping the flattest axis
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

        // Triangulate and emit
        let triangulated = earcut(&flat_coords, &hole_indices, 2)?;
        if triangulated.len() < 3 {
            return Ok(());
        }

        // Pre-allocate face indices (tip 1)
        let mut face_indices = Vec::with_capacity(triangulated.len());
        for idx in triangulated {
            face_indices.push(glb_indices[idx]);
        }

        self.emit_triangles(face_indices);
        Ok(())
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

    fn write_glb<P: AsRef<Path>>(&mut self, output_path: P, default_color: &str) -> Result<()> {
        self.normalize_normals();

        if self.positions.is_empty() {
            // Create parent directories if they don't exist
            if let Some(parent) = output_path.as_ref().parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("Failed to create parent directory for {:?}", output_path.as_ref()))?;
            }
            File::create(output_path.as_ref()).context("Create empty GLB file")?;
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

        let indices_offset = bin_buffer.len();
        for index in &self.indices {
            bin_buffer.extend_from_slice(&index.to_le_bytes());
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

        let accessor_indices = json::Accessor {
            buffer_view: Some(json::Index::new(2)),
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

        let buffer_views = vec![
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.positions.len() * 12) as u64),
                byte_offset: Some(json::validation::USize64(positions_offset as u64)),
                byte_stride: Some(json::buffer::Stride(12)),
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
            json::buffer::View {
                buffer: json::Index::new(0),
                byte_length: json::validation::USize64((self.normals.len() * 12) as u64),
                byte_offset: Some(json::validation::USize64(normals_offset as u64)),
                byte_stride: Some(json::buffer::Stride(12)),
                target: Some(json::validation::Checked::Valid(json::buffer::Target::ArrayBuffer)),
                extensions: Default::default(),
                extras: Default::default(),
                name: None,
            },
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

        let mut attributes = std::collections::BTreeMap::new();
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Positions),
            json::Index::new(0),
        );
        attributes.insert(
            json::validation::Checked::Valid(json::mesh::Semantic::Normals),
            json::Index::new(1),
        );

        // Create default PBR material matching pg2b3dm's structure
        let material = create_default_material(default_color)?;
        
        let primitive = json::mesh::Primitive {
            attributes,
            indices: Some(json::Index::new(2)),
            material: Some(json::Index::new(0)),
            mode: json::validation::Checked::Valid(json::mesh::Mode::Triangles),
            targets: None,
            extensions: Default::default(),
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

        let root = json::Root {
            accessors: vec![accessor_positions, accessor_normals, accessor_indices],
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

        // Create parent directories if they don't exist
        if let Some(parent) = output_path.as_ref().parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create parent directory for {:?}", output_path.as_ref()))?;
        }
        
        let mut file = File::create(output_path)?;
        file.write_all(&glb_bytes)?;
        
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

