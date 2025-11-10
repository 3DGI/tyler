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
    let mut builder = MeshBuilder::new(
        transformer_to_ecef,
        root_center_ecef,
        vertical_geoid_n
    );

    for cellid in qtree_node.cells() {
        let cell = world.grid.cell(cellid);
        for fid in cell.feature_ids.iter() {
            let feature = &world.features[*fid];
            let cf = CityJSONFeatureVertices::from_file(&feature.path_jsonl)
                .map_err(|e| anyhow::anyhow!("Failed to read {:?}: {}", feature.path_jsonl, e))?;
            builder.add_feature(&cf, &world.transform)?;
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
}

impl MeshBuilder {
    fn new(
        transformer_to_ecef: Proj,
        root_center_ecef: (f64, f64, f64),
        vertical_bias: f64,
    ) -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            transformer_to_ecef,
            root_center_ecef,
            vertical_bias,
        }
    }

    fn add_feature(&mut self, feature: &CityJSONFeatureVertices, transform: &Transform) -> Result<()> {
        let mut vertex_cache: HashMap<usize, u32> = HashMap::new();

        for (_id, co) in feature.cityobjects.iter() {
            if let Some(geoms) = &co.geometry {
                for geometry in geoms {
                    match geometry {
                        Geometry::MultiSurface { boundaries } => {
                            for surface in boundaries {
                                self.add_surface(surface, &feature.vertices, transform, &mut vertex_cache)?;
                            }
                        }
                        Geometry::Solid { boundaries } => {
                            for shell in boundaries {
                                for surface in shell {
                                    self.add_surface(surface, &feature.vertices, transform, &mut vertex_cache)?;
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
                let position = self.compute_local_position(vertex_id, vertices_qc, transform)?;
                let glb_index = self.vertex_index(vertex_id, position, cache);
                local_positions.push(position);
                glb_indices.push(glb_index);
                vertex_count += 1;
            }
        }

        if glb_indices.len() < 3 {
            return Ok(());
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
    ) -> u32 {
        if let Some(&existing) = cache.get(&idx) {
            return existing;
        }
        self.positions.push(position);
        self.normals.push([0.0, 0.0, 0.0]);
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

