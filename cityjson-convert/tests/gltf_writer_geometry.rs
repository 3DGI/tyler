use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use cityjson_convert::{convert_to_glb, ExportOptions, GeographicClipRegion, GeometryPlacement};
use cityjson_lib::json;
use serde_json::Value;

fn stable_output_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("output")
        .join(format!("{name}.glb"))
}

fn read_glb_json(bytes: &[u8]) -> Value {
    assert!(
        bytes.len() >= 20,
        "glb file should contain a header and JSON chunk"
    );
    assert_eq!(&bytes[0..4], b"glTF");
    assert_eq!(u32::from_le_bytes(bytes[4..8].try_into().unwrap()), 2);

    let declared_length = u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize;
    assert_eq!(declared_length, bytes.len());

    let json_length = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
    assert_eq!(&bytes[16..20], b"JSON");

    serde_json::from_slice(&bytes[20..20 + json_length]).expect("GLB JSON chunk should parse")
}

#[allow(clippy::cast_possible_truncation)]
fn read_glb_bin(bytes: &[u8]) -> &[u8] {
    let json_length = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
    let bin_header_offset = 20 + json_length;
    let bin_length = u32::from_le_bytes(
        bytes[bin_header_offset..bin_header_offset + 4]
            .try_into()
            .unwrap(),
    ) as usize;
    assert_eq!(
        &bytes[bin_header_offset + 4..bin_header_offset + 8],
        b"BIN\0"
    );
    &bytes[bin_header_offset + 8..bin_header_offset + 8 + bin_length]
}

#[allow(clippy::cast_possible_truncation)]
fn read_f32_vec3_accessor(root: &Value, bin: &[u8], accessor_index: usize) -> Vec<[f32; 3]> {
    let accessor = &root["accessors"][accessor_index];
    assert_eq!(accessor["componentType"].as_u64().unwrap(), 5126);
    assert_eq!(accessor["type"].as_str().unwrap(), "VEC3");
    let buffer_view_index = accessor["bufferView"].as_u64().unwrap() as usize;
    let buffer_view = &root["bufferViews"][buffer_view_index];
    let count = accessor["count"].as_u64().unwrap() as usize;
    let view_offset = buffer_view
        .get("byteOffset")
        .and_then(Value::as_u64)
        .unwrap_or(0) as usize;
    let accessor_offset = accessor
        .get("byteOffset")
        .and_then(Value::as_u64)
        .unwrap_or(0) as usize;
    let stride = buffer_view
        .get("byteStride")
        .and_then(Value::as_u64)
        .unwrap_or(12) as usize;
    let start = view_offset + accessor_offset;

    (0..count)
        .map(|index| {
            let offset = start + index * stride;
            [
                f32::from_le_bytes(bin[offset..offset + 4].try_into().unwrap()),
                f32::from_le_bytes(bin[offset + 4..offset + 8].try_into().unwrap()),
                f32::from_le_bytes(bin[offset + 8..offset + 12].try_into().unwrap()),
            ]
        })
        .collect()
}

fn normalize(vector: [f32; 3]) -> [f32; 3] {
    let length = (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
    [vector[0] / length, vector[1] / length, vector[2] / length]
}

fn ecef_from_lon_lat(lon_degrees: f64, lat_degrees: f64) -> [f64; 3] {
    let radius = 6_378_137.0_f64;
    let lon = lon_degrees.to_radians();
    let lat = lat_degrees.to_radians();
    [
        radius * lat.cos() * lon.cos(),
        radius * lat.cos() * lon.sin(),
        radius * lat.sin(),
    ]
}

fn lon_lat_from_ecef(position: [f64; 3]) -> [f64; 2] {
    let horizontal = position[0].hypot(position[1]);
    [
        position[1].atan2(position[0]).to_degrees(),
        position[2].atan2(horizontal).to_degrees(),
    ]
}

#[allow(clippy::uninlined_format_args)]
fn assert_vec3_approx_eq(actual: [f32; 3], expected: [f32; 3]) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() < 1.0e-5,
            "expected {:?}, got {:?}",
            expected,
            actual
        );
    }
}

#[allow(clippy::cast_possible_truncation)]
fn node_matrix_local_translation(node_matrix: &[Value]) -> [f32; 3] {
    [
        node_matrix[12].as_f64().unwrap() as f32,
        -node_matrix[14].as_f64().unwrap() as f32,
        node_matrix[13].as_f64().unwrap() as f32,
    ]
}

fn assert_positions_contain(positions: &[[f32; 3]], expected: [f32; 3]) {
    assert!(
        positions.iter().any(|actual| {
            actual
                .iter()
                .zip(expected)
                .all(|(actual, expected)| (*actual - expected).abs() < 1.0e-4)
        }),
        "expected positions to contain {expected:?}, got {positions:?}"
    );
}

fn has_extension(root: &Value, name: &str) -> bool {
    root["extensionsUsed"]
        .as_array()
        .is_some_and(|extensions| extensions.iter().any(|value| value.as_str() == Some(name)))
}

fn accessor_component_types(root: &Value) -> Vec<u64> {
    root["accessors"]
        .as_array()
        .expect("glTF should contain accessors")
        .iter()
        .map(|accessor| accessor["componentType"].as_u64().unwrap())
        .collect()
}

fn assert_rgba_eq(actual: &Value, expected: [f64; 4]) {
    let actual = actual
        .as_array()
        .expect("material baseColorFactor should be an array");
    assert_eq!(actual.len(), 4);
    for (component, expected) in actual.iter().zip(expected) {
        let actual = component.as_f64().unwrap();
        assert!(
            (actual - expected).abs() < 1.0e-6,
            "expected rgba component {expected}, got {actual}"
        );
    }
}

fn buffer_view_is_meshopt_compressed(root: &Value, buffer_view_index: u64) -> bool {
    root["bufferViews"][usize::try_from(buffer_view_index).unwrap()]["extensions"]
        ["EXT_meshopt_compression"]
        .is_object()
}

#[allow(clippy::too_many_lines, clippy::float_cmp)]
#[test]
fn convert_to_glb_writes_expected_geometry_layout() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let output_path = stable_output_path("ams-up-geometry");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    let accessors = root["accessors"]
        .as_array()
        .expect("glTF should contain accessors");
    assert_eq!(accessors.len(), 4);

    let positions = &accessors[0];
    let normals = &accessors[1];
    let feature_ids = &accessors[2];
    let indices = &accessors[3];

    let position_count = positions["count"].as_u64().unwrap();
    let normal_count = normals["count"].as_u64().unwrap();
    let index_count = indices["count"].as_u64().unwrap();

    assert!(position_count > 0, "mesh should contain vertices");
    assert_eq!(position_count, normal_count);
    assert_eq!(position_count, feature_ids["count"].as_u64().unwrap());
    assert!(
        index_count > position_count,
        "triangle soup should be deduplicated into indexed geometry"
    );
    assert_eq!(index_count % 3, 0, "index stream should describe triangles");
    assert_eq!(
        positions["componentType"].as_u64().unwrap(),
        5122,
        "positions should be quantized to i16"
    );
    assert_eq!(
        normals["componentType"].as_u64().unwrap(),
        5120,
        "normals should be quantized to i8"
    );
    assert!(positions["normalized"].as_bool().unwrap());
    assert!(normals["normalized"].as_bool().unwrap());
    assert!(matches!(
        feature_ids["componentType"].as_u64().unwrap(),
        5121 | 5123
    ));
    assert!(matches!(
        indices["componentType"].as_u64().unwrap(),
        5123 | 5125
    ));

    let extensions_used = root["extensionsUsed"]
        .as_array()
        .expect("quantized glTF should declare extensionsUsed");
    let extensions_required = root["extensionsRequired"]
        .as_array()
        .expect("quantized glTF should declare extensionsRequired");
    assert!(extensions_used
        .iter()
        .any(|value| value.as_str() == Some("KHR_mesh_quantization")));
    assert!(extensions_used
        .iter()
        .any(|value| value.as_str() == Some("EXT_meshopt_compression")));
    assert!(extensions_used
        .iter()
        .any(|value| value.as_str() == Some("EXT_mesh_features")));
    assert!(extensions_used
        .iter()
        .any(|value| value.as_str() == Some("EXT_structural_metadata")));
    assert!(extensions_required
        .iter()
        .any(|value| value.as_str() == Some("KHR_mesh_quantization")));
    assert!(extensions_required
        .iter()
        .any(|value| value.as_str() == Some("EXT_meshopt_compression")));

    let buffers = root["buffers"]
        .as_array()
        .expect("compressed glTF should contain explicit source and fallback buffers");
    assert_eq!(buffers.len(), 2);

    let min = positions["min"]
        .as_array()
        .expect("positions accessor should have min");
    let max = positions["max"]
        .as_array()
        .expect("positions accessor should have max");
    for axis in 0..3 {
        assert!(
            min[axis].as_i64().is_some(),
            "quantized POSITION min values must be raw integer components"
        );
        assert!(
            max[axis].as_i64().is_some(),
            "quantized POSITION max values must be raw integer components"
        );
        assert!(min[axis].as_f64().unwrap() < 0.0);
        assert!(max[axis].as_f64().unwrap() > 0.0);
    }

    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");
    assert_eq!(node_matrix.len(), 16);
    let dequant_scale = node_matrix[0].as_f64().unwrap();
    assert!(dequant_scale > 0.0);
    assert_eq!(node_matrix[5].as_f64().unwrap(), 0.0);
    assert_eq!(node_matrix[6].as_f64().unwrap(), -dequant_scale);
    assert_eq!(node_matrix[9].as_f64().unwrap(), dequant_scale);
    assert_eq!(node_matrix[10].as_f64().unwrap(), 0.0);
    assert_ne!(node_matrix[12].as_f64().unwrap(), 0.0);
    assert_ne!(node_matrix[13].as_f64().unwrap(), 0.0);
    assert_ne!(node_matrix[14].as_f64().unwrap(), 0.0);

    let mesh = &root["meshes"][0]["primitives"][0];
    if let Some(mode) = mesh.get("mode") {
        assert_eq!(mode.as_u64().unwrap(), 4);
    }
    assert_eq!(mesh["attributes"]["POSITION"].as_u64().unwrap(), 0);
    assert_eq!(mesh["attributes"]["NORMAL"].as_u64().unwrap(), 1);
    assert_eq!(mesh["attributes"]["_FEATURE_ID_0"].as_u64().unwrap(), 2);
    assert_eq!(mesh["indices"].as_u64().unwrap(), 3);
    assert_eq!(
        mesh["extensions"]["EXT_mesh_features"]["featureIds"][0]["propertyTable"]
            .as_u64()
            .unwrap(),
        0
    );
    assert!(
        root["extensions"]["EXT_structural_metadata"]["propertyTables"]
            .as_array()
            .is_some_and(|tables| !tables.is_empty())
    );
}

#[test]
fn convert_to_glb_omits_metadata_table_when_features_have_no_attributes() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building-1":{
                    "type":"Building",
                    "geometry":[{
                        "type":"MultiSurface",
                        "lod":"2.2",
                        "boundaries":[[[0,1,2]]]
                    }]
                }
            },
            "vertices":[
                [0.0,0.0,0.0],
                [1.0,0.0,0.0],
                [0.0,1.0,0.0]
            ]
        }"#,
    )
    .expect("inline CityJSON fixture should parse");
    let output_path = stable_output_path("no-attribute-metadata");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed without feature attributes");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    assert!(has_extension(&root, "EXT_mesh_features"));
    assert!(!has_extension(&root, "EXT_structural_metadata"));
    assert!(
        root.get("extensions").is_none(),
        "no attribute columns should mean no top-level structural metadata"
    );

    let feature_id =
        &root["meshes"][0]["primitives"][0]["extensions"]["EXT_mesh_features"]["featureIds"][0];
    assert_eq!(feature_id["attribute"].as_u64().unwrap(), 0);
    assert_eq!(feature_id["featureCount"].as_u64().unwrap(), 1);
    assert!(
        feature_id.get("propertyTable").is_none(),
        "mesh feature IDs must not reference a missing metadata property table"
    );
}

#[test]
fn convert_to_glb_can_reproject_geometry_to_ecef() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let output_path = stable_output_path("ams-up-ecef");
    let options = ExportOptions {
        geometry_placement: GeometryPlacement::EcefRelative {
            source_crs: "EPSG:7415".to_string(),
            origin: [0.0; 3],
        },
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options).expect("ECEF GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");

    let translation = [
        node_matrix[12].as_f64().unwrap(),
        node_matrix[13].as_f64().unwrap(),
        node_matrix[14].as_f64().unwrap(),
    ];
    assert!(
        translation
            .iter()
            .any(|component| component.abs() > 1_000_000.0),
        "ECEF translation should be in Earth-centered coordinates: {translation:?}"
    );
}

#[allow(clippy::float_cmp)]
#[test]
fn convert_to_glb_can_place_geometry_in_enu() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building-1":{
                    "type":"Building",
                    "geometry":[
                        {
                            "type":"MultiSurface",
                            "lod":"2.2",
                            "boundaries":[
                                [[0,1,2]]
                            ]
                        }
                    ]
                }
            },
            "vertices":[
                [1000.0,2000.0,3000.0],
                [1010.0,2000.0,3000.0],
                [1000.0,2020.0,3030.0]
            ]
        }"#,
    )
    .expect("inline CityJSON fixture should parse");
    let output_path = stable_output_path("enu-placement");
    let options = ExportOptions {
        geometry_placement: GeometryPlacement::Enu {
            source_crs: "EPSG:4978".to_string(),
            ecef_origin: [1000.0, 2000.0, 3000.0],
            east: [1.0, 0.0, 0.0],
            north: [0.0, 1.0, 0.0],
            up: [0.0, 0.0, 1.0],
        },
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options).expect("ENU GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let bin = read_glb_bin(&glb_bytes);
    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");

    assert_eq!(node_matrix[0].as_f64().unwrap(), 1.0);
    assert_eq!(node_matrix[5].as_f64().unwrap(), 0.0);
    assert_eq!(node_matrix[6].as_f64().unwrap(), -1.0);
    assert_eq!(node_matrix[9].as_f64().unwrap(), 1.0);
    assert_eq!(node_matrix[10].as_f64().unwrap(), 0.0);

    let center = node_matrix_local_translation(node_matrix);
    assert_vec3_approx_eq(center, [5.0, 10.0, 15.0]);
    assert!(
        center.iter().all(|component| component.abs() < 100.0),
        "ENU node translation should stay local-scale: {center:?}"
    );

    let centered_positions = read_f32_vec3_accessor(&root, bin, 0);
    let local_positions = centered_positions
        .into_iter()
        .map(|position| {
            [
                position[0] + center[0],
                position[1] + center[1],
                position[2] + center[2],
            ]
        })
        .collect::<Vec<_>>();

    assert_positions_contain(&local_positions, [0.0, 0.0, 0.0]);
    assert_positions_contain(&local_positions, [10.0, 0.0, 0.0]);
    assert_positions_contain(&local_positions, [0.0, 20.0, 30.0]);
}

#[allow(clippy::float_cmp)]
#[test]
fn convert_to_glb_can_disable_quantization_and_compression() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let output_path = stable_output_path("ams-up-raw");
    let options = ExportOptions {
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options).expect("raw GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    let component_types = accessor_component_types(&root);
    assert_eq!(component_types[0], 5126);
    assert_eq!(component_types[1], 5126);
    assert!(matches!(component_types[2], 5121 | 5123));
    assert!(matches!(component_types[3], 5123 | 5125));
    assert!(!has_extension(&root, "KHR_mesh_quantization"));
    assert!(!has_extension(&root, "EXT_meshopt_compression"));
    assert!(has_extension(&root, "EXT_mesh_features"));
    assert!(has_extension(&root, "EXT_structural_metadata"));

    let buffers = root["buffers"]
        .as_array()
        .expect("glTF should contain buffers");
    assert_eq!(
        buffers.len(),
        1,
        "uncompressed output should contain a single GLB buffer"
    );

    let buffer_views = root["bufferViews"]
        .as_array()
        .expect("glTF should contain bufferViews");
    for (index, buffer_view) in buffer_views.iter().enumerate() {
        assert_eq!(
            buffer_view["buffer"].as_u64().unwrap(),
            0,
            "bufferView {index} should point at the GLB BIN buffer"
        );
        assert!(
            buffer_view.get("extensions").is_none(),
            "bufferView {index} should not carry meshopt extension metadata"
        );
    }

    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");
    assert_eq!(node_matrix[0].as_f64().unwrap(), 1.0);
    assert_eq!(node_matrix[6].as_f64().unwrap(), -1.0);
    assert_eq!(node_matrix[9].as_f64().unwrap(), 1.0);
}

#[test]
fn convert_to_glb_can_disable_only_meshopt_compression() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let output_path = stable_output_path("ams-up-quantized-uncompressed");
    let options = ExportOptions {
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options)
        .expect("quantized uncompressed GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let accessors = root["accessors"]
        .as_array()
        .expect("glTF should contain accessors");

    assert_eq!(accessors[0]["componentType"].as_u64().unwrap(), 5122);
    assert_eq!(accessors[1]["componentType"].as_u64().unwrap(), 5120);
    assert!(matches!(
        accessors[2]["componentType"].as_u64().unwrap(),
        5121 | 5123
    ));
    assert!(has_extension(&root, "KHR_mesh_quantization"));
    assert!(!has_extension(&root, "EXT_meshopt_compression"));
    assert!(has_extension(&root, "EXT_mesh_features"));
    assert!(has_extension(&root, "EXT_structural_metadata"));
    assert_eq!(
        root["buffers"]
            .as_array()
            .expect("glTF should contain buffers")
            .len(),
        1,
        "disabling meshopt compression should remove the fallback buffer"
    );

    let buffer_views = root["bufferViews"]
        .as_array()
        .expect("glTF should contain bufferViews");
    for (index, buffer_view) in buffer_views.iter().enumerate() {
        assert_eq!(
            buffer_view["buffer"].as_u64().unwrap(),
            0,
            "bufferView {index} should point at the GLB BIN buffer"
        );
        assert!(
            buffer_view.get("extensions").is_none(),
            "bufferView {index} should not carry meshopt extension metadata"
        );
    }
}

#[test]
fn convert_to_glb_can_generate_smooth_normals() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let hard_output_path = stable_output_path("ams-up-hard-normals");
    let smooth_output_path = stable_output_path("ams-up-smooth-normals");

    let hard_options = ExportOptions {
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };
    let smooth_options = ExportOptions {
        smooth_normals: true,
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &hard_output_path, &hard_options)
        .expect("hard-normal GLB conversion should succeed");
    convert_to_glb(&model, &smooth_output_path, &smooth_options)
        .expect("smooth-normal GLB conversion should succeed");

    let hard_glb_bytes = fs::read(&hard_output_path).expect("hard-normal test GLB should exist");
    let smooth_glb_bytes =
        fs::read(&smooth_output_path).expect("smooth-normal test GLB should exist");
    let hard_root = read_glb_json(&hard_glb_bytes);
    let smooth_root = read_glb_json(&smooth_glb_bytes);

    let hard_accessors = hard_root["accessors"]
        .as_array()
        .expect("hard-normal glTF should contain accessors");
    let smooth_accessors = smooth_root["accessors"]
        .as_array()
        .expect("smooth-normal glTF should contain accessors");

    let hard_vertex_count = hard_accessors[0]["count"].as_u64().unwrap();
    let smooth_vertex_count = smooth_accessors[0]["count"].as_u64().unwrap();
    let hard_index_count = hard_accessors[3]["count"].as_u64().unwrap();
    let smooth_index_count = smooth_accessors[3]["count"].as_u64().unwrap();

    assert!(
        smooth_vertex_count < hard_vertex_count,
        "smooth normals should deduplicate shared vertices"
    );
    assert_eq!(
        smooth_index_count, hard_index_count,
        "smoothing should preserve triangle topology"
    );
}

#[test]
fn convert_to_glb_can_clip_to_bbox_with_preclip_smooth_normals() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building-1":{
                    "type":"Building",
                    "geometry":[
                        {
                            "type":"MultiSurface",
                            "lod":"2.2",
                            "boundaries":[
                                [[0,1,2]],
                                [[1,3,2]]
                            ]
                        }
                    ]
                }
            },
            "vertices":[
                [0.0,0.0,0.0],
                [1.0,0.0,0.0],
                [0.0,1.0,0.0],
                [1.0,1.0,1.0]
            ]
        }"#,
    )
    .expect("inline CityJSON fixture should parse");
    let output_path = stable_output_path("clipped-preclip-smooth-normals");
    let options = ExportOptions {
        clip_bbox: Some([-1.0, -1.0, -1.0, 0.5, 2.0, 2.0]),
        smooth_normals: true,
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options)
        .expect("clipped smooth-normal GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("clipped test GLB should exist");
    let root = read_glb_json(&glb_bytes);
    let bin = read_glb_bin(&glb_bytes);
    let positions = read_f32_vec3_accessor(&root, bin, 0);
    let normals = read_f32_vec3_accessor(&root, bin, 1);

    assert_eq!(
        positions.len(),
        5,
        "clipped mesh should emit five unique vertices"
    );
    assert_eq!(positions.len(), normals.len());

    let clipped_boundary_index = positions
        .iter()
        .position(|position| {
            (position[0] - 0.25).abs() < 1.0e-6
                && (position[1] + 0.5).abs() < 1.0e-6
                && (position[2] + 0.25).abs() < 1.0e-6
        })
        .expect("clipped boundary vertex should be present");

    let first_face_normal = [0.0, 0.0, 1.0];
    let second_face_normal = normalize([-1.0, -1.0, 1.0]);
    let shared_normal = normalize([
        first_face_normal[0] + second_face_normal[0],
        first_face_normal[1] + second_face_normal[1],
        first_face_normal[2] + second_face_normal[2],
    ]);
    let expected_boundary_normal = normalize([
        0.5 * first_face_normal[0] + 0.5 * shared_normal[0],
        0.5 * first_face_normal[1] + 0.5 * shared_normal[1],
        0.5 * first_face_normal[2] + 0.5 * shared_normal[2],
    ]);

    assert_vec3_approx_eq(normals[clipped_boundary_index], expected_boundary_normal);
}

#[test]
fn convert_to_glb_can_clip_to_geographic_region() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "surface-1":{
                    "type":"BuildingPart",
                    "geometry":[
                        {
                            "type":"MultiSurface",
                            "lod":"2.2",
                            "boundaries":[
                                [[0,1,2]],
                                [[1,3,2]]
                            ]
                        }
                    ]
                }
            },
            "vertices":[
                [0.0,0.0,0.0],
                [2.0,0.0,0.0],
                [0.0,2.0,0.0],
                [2.0,2.0,0.0]
            ]
        }"#,
    )
    .expect("inline CityJSON fixture should parse");
    let output_path = stable_output_path("clipped-geographic-region");
    let options = ExportOptions {
        clip_geographic_region: Some(GeographicClipRegion {
            source_crs: "EPSG:4979".to_string(),
            west: 0.0,
            south: 0.0,
            east: 1.0,
            north: 1.0,
        }),
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options)
        .expect("geographically clipped GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("clipped test GLB should exist");
    let root = read_glb_json(&glb_bytes);
    let bin = read_glb_bin(&glb_bytes);
    let positions = read_f32_vec3_accessor(&root, bin, 0);
    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");
    let center = node_matrix_local_translation(node_matrix);
    let mut touches_east = false;
    let mut touches_north = false;

    for position in &positions {
        let world_x = position[0] + center[0];
        let world_y = position[1] + center[1];
        assert!((0.0..=1.0).contains(&world_x), "x outside clip: {world_x}");
        assert!((0.0..=1.0).contains(&world_y), "y outside clip: {world_y}");
        touches_east |= (world_x - 1.0).abs() < 1.0e-6;
        touches_north |= (world_y - 1.0).abs() < 1.0e-6;
    }

    assert!(
        touches_east,
        "clipped mesh should include the east boundary"
    );
    assert!(
        touches_north,
        "clipped mesh should include the north boundary"
    );
}

#[test]
fn convert_to_glb_solves_geographic_clip_intersections_in_source_space() {
    let [a_x, a_y, a_z] = ecef_from_lon_lat(0.0, 0.0);
    let [b_x, b_y, b_z] = ecef_from_lon_lat(10.0, 0.0);
    let [c_x, c_y, c_z] = ecef_from_lon_lat(0.0, 1.0);
    let model = json::from_slice(
        format!(
            r#"{{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{{
                "surface-1":{{
                    "type":"TINRelief",
                    "geometry":[
                        {{
                            "type":"MultiSurface",
                            "lod":"1.0",
                            "boundaries":[
                                [[0,1,2]]
                            ]
                        }}
                    ]
                }}
            }},
            "vertices":[
                [{a_x},{a_y},{a_z}],
                [{b_x},{b_y},{b_z}],
                [{c_x},{c_y},{c_z}]
            ]
        }}"#
        )
        .as_bytes(),
    )
    .expect("inline ECEF CityJSON fixture should parse");
    let output_path = stable_output_path("clipped-geographic-ecef");
    let options = ExportOptions {
        clip_geographic_region: Some(GeographicClipRegion {
            source_crs: "EPSG:4978".to_string(),
            west: -1.0,
            south: -1.0,
            east: 2.0,
            north: 2.0,
        }),
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options)
        .expect("geographically clipped ECEF GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("clipped test GLB should exist");
    let root = read_glb_json(&glb_bytes);
    let bin = read_glb_bin(&glb_bytes);
    let positions = read_f32_vec3_accessor(&root, bin, 0);
    let node_matrix = root["nodes"][0]["matrix"]
        .as_array()
        .expect("root node should carry the local-to-world transform");
    let center = node_matrix_local_translation(node_matrix);
    let mut touches_east = false;

    for position in &positions {
        let source_position = [
            f64::from(position[0] + center[0]),
            f64::from(position[1] + center[1]),
            f64::from(position[2] + center[2]),
        ];
        let [lon, lat] = lon_lat_from_ecef(source_position);
        assert!(lon <= 2.0001, "longitude outside clip: {lon}");
        assert!(lat <= 2.0001, "latitude outside clip: {lat}");
        touches_east |= (lon - 2.0).abs() < 0.0001;
    }

    assert!(
        touches_east,
        "clipped mesh should include a source-space intersection on the east boundary"
    );
}

#[test]
fn convert_to_glb_merges_multiple_features_into_one_centered_mesh() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture feature stream should parse");
    let output_path = stable_output_path("multi-feature-types");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    let primitives = root["meshes"][0]["primitives"]
        .as_array()
        .expect("glTF should contain primitives");
    assert_eq!(
        primitives.len(),
        2,
        "writer should group geometry into one primitive per feature type"
    );

    let materials = root["materials"]
        .as_array()
        .expect("glTF should contain materials");
    assert_eq!(
        materials.len(),
        2,
        "writer should create one material per primitive type group"
    );
    assert_ne!(
        materials[0]["pbrMetallicRoughness"]["baseColorFactor"],
        materials[1]["pbrMetallicRoughness"]["baseColorFactor"],
        "different feature types should receive different default colors"
    );

    for primitive in primitives {
        assert!(primitive["attributes"].get("_FEATURE_ID_0").is_some());
        assert_eq!(
            primitive["extensions"]["EXT_mesh_features"]["featureIds"][0]["featureCount"]
                .as_u64()
                .unwrap(),
            1
        );
    }

    assert_eq!(
        root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["count"]
            .as_u64()
            .unwrap(),
        2
    );
}

#[test]
fn convert_to_glb_writes_meshopt_compression_extension() {
    let model = json::merge_feature_stream_slice(include_bytes!("data/ams_up_holes.city.jsonl"))
        .expect("fixture feature stream should parse");
    let output_path = stable_output_path("ams-up-meshopt");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    let buffers = root["buffers"]
        .as_array()
        .expect("compressed glTF should declare buffers");
    assert_eq!(
        buffers.len(),
        2,
        "compressed output should carry a source buffer and a fallback placeholder"
    );
    assert!(
        buffers[1]["extensions"]["EXT_meshopt_compression"]["fallback"]
            .as_bool()
            .unwrap(),
        "fallback buffer should be explicitly tagged"
    );

    let buffer_views = root["bufferViews"]
        .as_array()
        .expect("compressed glTF should declare bufferViews");
    let meshopt_views = buffer_views
        .iter()
        .filter(|buffer_view| buffer_view.get("extensions").is_some())
        .collect::<Vec<_>>();
    assert!(
        meshopt_views.len() >= 4,
        "compressed output should include at least the geometry streams"
    );

    for (index, buffer_view) in meshopt_views.iter().enumerate() {
        assert_eq!(
            buffer_view["buffer"].as_u64().unwrap(),
            1,
            "bufferView {index} should reference the fallback buffer layout"
        );
        let extension = &buffer_view["extensions"]["EXT_meshopt_compression"];
        assert_eq!(
            extension["buffer"].as_u64().unwrap(),
            0,
            "bufferView {index} should source compressed bytes from buffer 0"
        );
        assert!(extension["byteLength"].as_u64().unwrap() > 0);
        assert!(extension["count"].as_u64().unwrap() > 0);
    }

    let accessors = root["accessors"]
        .as_array()
        .expect("compressed glTF should declare accessors");
    let geometry_buffer_views = [
        accessors[0]["bufferView"].as_u64().unwrap(),
        accessors[1]["bufferView"].as_u64().unwrap(),
        accessors[2]["bufferView"].as_u64().unwrap(),
        accessors[3]["bufferView"].as_u64().unwrap(),
    ];
    assert_eq!(
        buffer_views[usize::try_from(geometry_buffer_views[0]).unwrap()]["extensions"]
            ["EXT_meshopt_compression"]["mode"]
            .as_str()
            .unwrap(),
        "ATTRIBUTES"
    );
    assert_eq!(
        buffer_views[usize::try_from(geometry_buffer_views[1]).unwrap()]["extensions"]
            ["EXT_meshopt_compression"]["mode"]
            .as_str()
            .unwrap(),
        "ATTRIBUTES"
    );
    assert_eq!(
        buffer_views[usize::try_from(geometry_buffer_views[2]).unwrap()]["extensions"]
            ["EXT_meshopt_compression"]["mode"]
            .as_str()
            .unwrap(),
        "ATTRIBUTES"
    );
    assert_eq!(
        buffer_views[usize::try_from(geometry_buffer_views[3]).unwrap()]["extensions"]
            ["EXT_meshopt_compression"]["mode"]
            .as_str()
            .unwrap(),
        "TRIANGLES"
    );
}

#[test]
fn convert_to_glb_uses_custom_metadata_class_and_feature_colors() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture feature stream should parse");
    let output_path = stable_output_path("multi-feature-types-custom");
    let options = ExportOptions {
        metadata_class_name: "test".to_string(),
        feature_type_colors: BTreeMap::from([
            ("Building".to_string(), "#010203".to_string()),
            ("WaterBody".to_string(), "#ABCDEF".to_string()),
        ]),
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options).expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    assert!(
        root["extensions"]["EXT_structural_metadata"]["schema"]["classes"]
            .get("test")
            .is_some()
    );
    assert_eq!(
        root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["class"]
            .as_str()
            .unwrap(),
        "test"
    );

    let materials = root["materials"]
        .as_array()
        .expect("glTF should contain materials");
    assert_rgba_eq(
        &materials[0]["pbrMetallicRoughness"]["baseColorFactor"],
        [1.0 / 255.0, 2.0 / 255.0, 3.0 / 255.0, 1.0],
    );
    assert_rgba_eq(
        &materials[1]["pbrMetallicRoughness"]["baseColorFactor"],
        [171.0 / 255.0, 205.0 / 255.0, 239.0 / 255.0, 1.0],
    );
}

#[test]
fn convert_to_glb_compresses_metadata_numeric_columns() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture feature stream should parse");
    let output_path = stable_output_path("multi-feature-types-metadata-compressed");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let properties = root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]
        ["properties"]
        .as_object()
        .expect("structural metadata property table should contain properties");

    let levels_values = properties["levels"]["values"].as_u64().unwrap();
    let depth_values = properties["depth"]["values"].as_u64().unwrap();
    let name_offsets = properties["name"]["stringOffsets"].as_u64().unwrap();
    let name_values = properties["name"]["values"].as_u64().unwrap();
    let occupied_values = properties["occupied"]["values"].as_u64().unwrap();

    assert!(buffer_view_is_meshopt_compressed(&root, levels_values));
    assert!(buffer_view_is_meshopt_compressed(&root, depth_values));
    assert!(buffer_view_is_meshopt_compressed(&root, name_offsets));
    assert!(
        !buffer_view_is_meshopt_compressed(&root, name_values),
        "string byte blobs still use raw buffer views"
    );
    assert!(
        !buffer_view_is_meshopt_compressed(&root, occupied_values),
        "INT8 bool columns remain raw because meshopt attribute encoding requires 4-byte elements"
    );
}

#[test]
fn convert_to_glb_writes_all_geometry_it_receives() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_lod_building_part.city.jsonl"))
            .expect("fixture feature stream should parse");
    let output_path = stable_output_path("multi-lod-building-part");
    let options = ExportOptions {
        quantize_geometry: false,
        meshopt_compression: false,
        ..ExportOptions::default()
    };

    convert_to_glb(&model, &output_path, &options).expect("GLB conversion should succeed");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let accessors = root["accessors"]
        .as_array()
        .expect("glTF should contain accessors");

    assert!(
        accessors[3]["count"].as_u64().unwrap() > 36,
        "the writer should serialize every geometry left in the input model"
    );
    assert_eq!(
        root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["count"]
            .as_u64()
            .unwrap(),
        1,
        "metadata rows should stay aligned with the surviving feature geometry"
    );
}
