use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use cityjson_convert::{
    convert_to_cityjson, convert_to_cityjsonseq, convert_to_glb, convert_to_obj,
    CityJsonSeqExportOptions, ExportOptions, GeographicClipRegion, GeometryPlacement,
    JsonExportOptions, ObjExportOptions, VoxelRemeshOptions,
};
use cityjson_lib::json;
use serde_json::Value;

fn stable_output_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("output")
        .join(format!("{name}.glb"))
}

fn stable_output_path_with_suffix(name: &str, suffix: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("output")
        .join(format!("{name}.{suffix}"))
}

fn legacy_feature_dataset_path(name: &str) -> PathBuf {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("output")
        .join(name);
    if path.exists() {
        fs::remove_dir_all(&path).expect("remove previous legacy feature dataset");
    }
    fs::create_dir_all(&path).expect("create legacy feature dataset");

    let items = include_bytes!("data/multi_feature_types.city.jsonl")
        .split(|byte| *byte == b'\n')
        .filter(|item| !item.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(items.len(), 3);
    let mut metadata: Value =
        serde_json::from_slice(items[0]).expect("parse dataset metadata fixture");
    metadata["transform"]["translate"] = serde_json::json!([100.0, 200.0, 0.0]);
    fs::write(
        path.join("metadata.json"),
        serde_json::to_vec(&metadata).expect("serialize dataset metadata"),
    )
    .expect("write dataset metadata");
    fs::write(path.join("building.city.jsonl"), items[1]).expect("write building feature");
    fs::write(path.join("water.city.jsonl"), items[2]).expect("write water feature");
    path
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

fn read_json_lines(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .expect("read json lines")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str(line).expect("json line should parse"))
        .collect()
}

fn read_obj_lines(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .expect("read OBJ output")
        .lines()
        .map(str::to_string)
        .collect()
}

fn read_obj_vertices_and_faces(path: &Path) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    for line in read_obj_lines(path) {
        if let Some(rest) = line.strip_prefix("v ") {
            let coordinates = rest
                .split_whitespace()
                .map(|value| value.parse::<f64>().expect("OBJ coordinate should parse"))
                .collect::<Vec<_>>();
            assert_eq!(coordinates.len(), 3);
            vertices.push([coordinates[0], coordinates[1], coordinates[2]]);
        } else if let Some(rest) = line.strip_prefix("f ") {
            let indices = rest
                .split_whitespace()
                .map(|value| value.parse::<usize>().expect("OBJ face index should parse") - 1)
                .collect::<Vec<_>>();
            assert_eq!(indices.len(), 3);
            faces.push([indices[0], indices[1], indices[2]]);
        }
    }
    (vertices, faces)
}

fn inline_feature(id: &str, vertices: &str, boundaries: &str) -> Vec<u8> {
    format!(
        r#"{{"type":"CityJSONFeature","id":"{id}","transform":{{"scale":[1.0,1.0,1.0],"translate":[0.0,0.0,0.0]}},"CityObjects":{{"{id}":{{"type":"Building","geometry":[{{"type":"MultiSurface","lod":"1","boundaries":{boundaries}}}]}}}},"vertices":{vertices}}}"#
    )
    .into_bytes()
}

fn face_normal(vertices: &[[f64; 3]], face: [usize; 3]) -> [f64; 3] {
    let p0 = vertices[face[0]];
    let p1 = vertices[face[1]];
    let p2 = vertices[face[2]];
    let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn dot_f64(lhs: [f64; 3], rhs: [f64; 3]) -> f64 {
    lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
}

#[test]
fn convert_to_obj_writes_single_triangle_in_source_coordinates() {
    let input = inline_feature(
        "source triangle",
        "[[1000.5,2000.25,3.0],[1001.5,2000.25,3.0],[1000.5,2001.25,4.0]]",
        "[[[0,1,2]]]",
    );
    let model = json::from_feature_slice(&input).expect("inline feature should parse");
    let output_path = stable_output_path_with_suffix("convert_to_obj_single_triangle", "obj");

    convert_to_obj(&model, &output_path, &ObjExportOptions::default())
        .expect("OBJ conversion should succeed");

    let lines = read_obj_lines(&output_path);
    assert_eq!(lines[0], "o source_triangle");
    assert_eq!(lines[1], "v 1000.5 2000.25 3");
    assert_eq!(lines[2], "v 1001.5 2000.25 3");
    assert_eq!(lines[3], "v 1000.5 2001.25 4");
    assert_eq!(lines[4], "f 1 2 3");
}

#[test]
fn convert_to_obj_groups_cityobjects_and_uses_global_indices() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture should merge");
    let output_path = stable_output_path_with_suffix("convert_to_obj_multi_object", "obj");

    convert_to_obj(&model, &output_path, &ObjExportOptions::default())
        .expect("OBJ conversion should succeed");

    let lines = read_obj_lines(&output_path);
    let object_lines = lines
        .iter()
        .filter(|line| line.starts_with("o "))
        .collect::<Vec<_>>();
    assert!(object_lines
        .iter()
        .any(|line| *line == "o building-feature"));
    assert!(object_lines.iter().any(|line| *line == "o water-feature"));
    assert!(lines.iter().any(|line| line == "f 1 2 3"));
    assert!(lines
        .iter()
        .filter_map(|line| line.strip_prefix("f "))
        .any(|face| face != "1 2 3"));
}

#[test]
fn convert_to_obj_triangulates_polygons_and_preserves_unclipped_geometry_by_default() {
    let input = inline_feature(
        "quad",
        "[[0.0,0.0,0.0],[2.0,0.0,0.0],[2.0,2.0,0.0],[0.0,2.0,0.0]]",
        "[[[0,1,2,3]]]",
    );
    let model = json::from_feature_slice(&input).expect("inline feature should parse");
    let output_path = stable_output_path_with_suffix("convert_to_obj_quad", "obj");

    convert_to_obj(&model, &output_path, &ObjExportOptions::default())
        .expect("OBJ conversion should succeed");

    let lines = read_obj_lines(&output_path);
    assert_eq!(
        lines.iter().filter(|line| line.starts_with("f ")).count(),
        2
    );
    assert!(lines.iter().any(|line| line == "v 2 2 0"));
}

#[test]
fn convert_to_obj_preserves_winding_for_y_dominant_polygons() {
    let input = inline_feature(
        "y-normal-quad",
        "[[0.0,0.0,0.0],[0.0,0.0,2.0],[2.0,0.0,2.0],[2.0,0.0,0.0]]",
        "[[[0,1,2,3]]]",
    );
    let model = json::from_feature_slice(&input).expect("inline feature should parse");
    let output_path = stable_output_path_with_suffix("convert_to_obj_y_winding", "obj");

    convert_to_obj(&model, &output_path, &ObjExportOptions::default())
        .expect("OBJ conversion should succeed");

    let (vertices, faces) = read_obj_vertices_and_faces(&output_path);
    assert_eq!(faces.len(), 2);
    for face in faces {
        assert!(
            dot_f64(face_normal(&vertices, face), [0.0, 1.0, 0.0]) > 0.0,
            "OBJ face {face:?} should preserve the source +Y winding"
        );
    }
}

#[test]
fn convert_to_obj_preserves_reversed_winding_for_y_dominant_polygons() {
    let input = inline_feature(
        "negative-y-normal-quad",
        "[[0.0,0.0,0.0],[2.0,0.0,0.0],[2.0,0.0,2.0],[0.0,0.0,2.0]]",
        "[[[0,1,2,3]]]",
    );
    let model = json::from_feature_slice(&input).expect("inline feature should parse");
    let output_path = stable_output_path_with_suffix("convert_to_obj_negative_y_winding", "obj");

    convert_to_obj(&model, &output_path, &ObjExportOptions::default())
        .expect("OBJ conversion should succeed");

    let (vertices, faces) = read_obj_vertices_and_faces(&output_path);
    assert_eq!(faces.len(), 2);
    for face in faces {
        assert!(
            dot_f64(face_normal(&vertices, face), [0.0, -1.0, 0.0]) > 0.0,
            "OBJ face {face:?} should preserve the source -Y winding"
        );
    }
}

#[test]
fn convert_to_obj_can_clip_to_bbox() {
    let input = inline_feature(
        "clip",
        "[[-1.0,0.0,0.0],[2.0,0.0,0.0],[0.0,2.0,0.0]]",
        "[[[0,1,2]]]",
    );
    let model = json::from_feature_slice(&input).expect("inline feature should parse");
    let output_path = stable_output_path_with_suffix("convert_to_obj_clipped", "obj");
    let options = ObjExportOptions {
        clip_bbox: Some([0.0, 0.0, -1.0, 1.0, 1.0, 1.0]),
    };

    convert_to_obj(&model, &output_path, &options).expect("OBJ conversion should succeed");

    let lines = read_obj_lines(&output_path);
    let vertices = lines
        .iter()
        .filter_map(|line| line.strip_prefix("v "))
        .map(|line| {
            line.split_whitespace()
                .map(|value| value.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert!(!vertices.is_empty());
    assert!(vertices.iter().all(|vertex| {
        vertex[0] >= -1.0e-9
            && vertex[0] <= 1.0 + 1.0e-9
            && vertex[1] >= -1.0e-9
            && vertex[1] <= 1.0 + 1.0e-9
    }));
}

#[test]
fn convert_to_cityjson_writes_parseable_cityjson() {
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture should merge");
    let output_path = stable_output_path_with_suffix("convert_to_cityjson", "city.json");

    convert_to_cityjson(&model, &output_path, &JsonExportOptions::default())
        .expect("CityJSON conversion should succeed");

    let root: Value =
        serde_json::from_slice(&fs::read(&output_path).expect("read CityJSON output"))
            .expect("CityJSON output should parse");
    assert_eq!(root["type"], "CityJSON");
}

#[test]
fn convert_to_cityjsonseq_writes_header_and_features() {
    let base = json::from_slice(
        include_bytes!("data/multi_feature_types.city.jsonl")
            .split(|b| *b == b'\n')
            .next()
            .unwrap(),
    )
    .expect("cityjson base should parse");
    let feature_a = json::from_feature_slice(
        include_bytes!("data/multi_feature_types.city.jsonl")
            .split(|b| *b == b'\n')
            .nth(1)
            .unwrap(),
    )
    .expect("first feature should parse");
    let feature_b = json::from_feature_slice(
        include_bytes!("data/multi_feature_types.city.jsonl")
            .split(|b| *b == b'\n')
            .nth(2)
            .unwrap(),
    )
    .expect("second feature should parse");
    let features = vec![feature_a, feature_b];
    let output_path = stable_output_path_with_suffix("convert_to_cityjsonseq", "city.jsonl");

    convert_to_cityjsonseq(
        &base,
        &features,
        &output_path,
        &CityJsonSeqExportOptions::default(),
    )
    .expect("CityJSONSeq conversion should succeed");

    let items = read_json_lines(&output_path);
    assert_eq!(items[0]["type"], "CityJSON");
    assert_eq!(items[1]["type"], "CityJSONFeature");
    assert_eq!(items[2]["type"], "CityJSONFeature");
}

fn assert_cjconvert_lod_selection() {
    let directory = tempfile::tempdir().unwrap();
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_lod_building_part.city.jsonl"))
            .unwrap();
    let input = directory.path().join("input.city.json");
    convert_to_cityjson(&model, &input, &JsonExportOptions::default()).unwrap();
    let mut document: Value = serde_json::from_slice(&fs::read(&input).unwrap()).unwrap();
    let child = document["CityObjects"]["multi-lod-building-part"]
        .as_object_mut()
        .unwrap();
    let attributes = child.remove("attributes").unwrap();
    child.insert("parents".into(), serde_json::json!(["parent-building"]));
    document["CityObjects"]["parent-building"] = serde_json::json!({
        "type": "Building",
        "attributes": attributes,
        "children": ["multi-lod-building-part"]
    });
    fs::write(&input, serde_json::to_vec(&document).unwrap()).unwrap();
    let mut counts = Vec::new();
    for (name, args) in [
        ("all", Vec::new()),
        ("selected", vec!["--lod", "2.2"]),
        ("proxy", vec!["--lod", "2.2", "--remesh-error", "4"]),
        ("proxy-all", vec!["--remesh-error", "4"]),
    ] {
        let output = directory.path().join(format!("{name}.glb"));
        let result = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
            .arg(&input)
            .arg("--output")
            .arg(&output)
            .args(["--no-quantization", "--no-meshopt-compression"])
            .args(args)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let root = read_glb_json(&fs::read(output).unwrap());
        let properties =
            &root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["properties"];
        assert!(
            properties.get("name").is_some(),
            "parent metadata must survive LoD selection"
        );
        let primitive = &root["meshes"][0]["primitives"][0];
        let accessor = usize::try_from(primitive["indices"].as_u64().unwrap()).unwrap();
        counts.push(root["accessors"][accessor]["count"].as_u64().unwrap());
        if name.starts_with("proxy") {
            assert!(String::from_utf8_lossy(&result.stderr).contains("estimated error"));
        }
    }
    assert!(counts[0] > counts[1]);
    for args in [
        vec!["--lod", "999"],
        vec!["--lod", "2.2", "--remesh-error", "NaN"],
        vec!["--lod", "2.2", "--remesh-error", "0"],
        vec!["--format", "obj", "--remesh-error", "4"],
    ] {
        let output = directory.path().join("invalid.glb");
        let result = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
            .arg(&input)
            .arg("--output")
            .arg(&output)
            .args(args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(!output.exists());
    }
}

#[test]
fn cjconvert_can_write_all_formats_with_optional_lod_selection() {
    let cityjson_output = stable_output_path_with_suffix("cjconvert_cityjson", "city.json");
    let cityjsonseq_output = stable_output_path_with_suffix("cjconvert_cityjsonseq", "city.jsonl");
    let obj_output = stable_output_path_with_suffix("cjconvert_obj", "obj");
    let cityjson_input = stable_output_path_with_suffix("cjconvert_input", "city.json");
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_feature_types.city.jsonl"))
            .expect("fixture should merge");
    convert_to_cityjson(&model, &cityjson_input, &JsonExportOptions::default())
        .expect("write cjconvert input");

    let cityjson_status = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
        .arg(&cityjson_input)
        .arg("--output")
        .arg(&cityjson_output)
        .arg("--format")
        .arg("cityjson")
        .status()
        .expect("run cjconvert cityjson");
    assert!(cityjson_status.success());
    let root: Value =
        serde_json::from_slice(&fs::read(&cityjson_output).expect("read cjconvert CityJSON"))
            .expect("parse cjconvert CityJSON");
    assert_eq!(root["type"], "CityJSON");

    let cityjsonseq_status = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
        .arg(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/multi_feature_types.city.jsonl"),
        )
        .arg("--output")
        .arg(&cityjsonseq_output)
        .arg("--format")
        .arg("cityjsonseq")
        .status()
        .expect("run cjconvert cityjsonseq");
    assert!(cityjsonseq_status.success());
    let items = read_json_lines(&cityjsonseq_output);
    assert_eq!(items[0]["type"], "CityJSON");
    assert_eq!(items[1]["type"], "CityJSONFeature");

    let obj_status = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
        .arg(&cityjson_input)
        .arg("--output")
        .arg(&obj_output)
        .arg("--format")
        .arg("obj")
        .status()
        .expect("run cjconvert obj");
    assert!(obj_status.success());
    let obj_lines = read_obj_lines(&obj_output);
    assert!(obj_lines.iter().any(|line| line.starts_with("o ")));
    assert!(obj_lines.iter().any(|line| line.starts_with("v ")));
    assert!(obj_lines.iter().any(|line| line.starts_with("f ")));
    assert_cjconvert_lod_selection();
    for format in ["cityjson", "cityjsonseq", "obj", "tsv", "gpkg"] {
        assert_cjconvert_format_lod_selection(format);
    }
}

fn assert_cjconvert_format_lod_selection(format: &str) {
    let directory = tempfile::tempdir().unwrap();
    let model =
        json::merge_feature_stream_slice(include_bytes!("data/multi_lod_building_part.city.jsonl"))
            .unwrap();
    let input = directory.path().join("input.city.json");
    convert_to_cityjson(&model, &input, &JsonExportOptions::default()).unwrap();
    let mut document: Value = serde_json::from_slice(&fs::read(&input).unwrap()).unwrap();
    document["metadata"] =
        serde_json::json!({"referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/7415"});
    for geometry in document["CityObjects"]["multi-lod-building-part"]["geometry"]
        .as_array_mut()
        .unwrap()
    {
        let values = if geometry["type"] == "Solid" {
            serde_json::json!([[0, 0, 0, 0, 0, 0]])
        } else {
            serde_json::json!([0])
        };
        geometry["semantics"] = serde_json::json!({
            "surfaces": [{"type": "WallSurface"}], "values": values
        });
    }
    fs::write(&input, serde_json::to_vec(&document).unwrap()).unwrap();
    let input = if format == "cityjsonseq" {
        let mut feature: Value = serde_json::from_slice(&fs::read(&input).unwrap()).unwrap();
        feature["type"] = serde_json::json!("CityJSONFeature");
        feature["id"] = serde_json::json!("multi-lod-building-part");
        let base = serde_json::json!({
            "type": "CityJSON", "version": "2.0", "CityObjects": {}, "vertices": [],
            "transform": {"scale": [1, 1, 1], "translate": [0, 0, 0]}
        });
        let unmatched = serde_json::json!({
            "type": "CityJSONFeature", "id": "unmatched", "vertices": [[0,0,0],[1,0,0],[0,1,0]],
            "CityObjects": {"unmatched": {"type": "BuildingPart", "attributes": {"name": "keep"},
                "geometry": [{"type": "MultiSurface", "lod": "1.2", "boundaries": [[[0,1,2]]]}]}}
        });
        let sequence = directory.path().join("input.city.jsonl");
        fs::write(&sequence, format!("{base}\n{feature}\n{unmatched}\n")).unwrap();
        sequence
    } else {
        input
    };
    let output = directory.path().join(format!("selected.{format}"));
    let result = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
        .arg(&input)
        .args([
            "--format",
            format,
            "--lod",
            "2.2",
            "--tsv-include-semantics",
        ])
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_selected_lod_output(format, &output);
    let result = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
        .arg(&input)
        .args(["--format", format, "--lod", "999", "--output"])
        .arg(directory.path().join("missing"))
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("matched no geometry"));
    assert!(!directory.path().join("missing").exists());
}

fn assert_selected_lod_output(format: &str, output: &Path) {
    match format {
        "cityjson" | "cityjsonseq" => {
            let documents = if format == "cityjsonseq" {
                let items = read_json_lines(output);
                assert_eq!(items.len(), 3);
                assert_eq!(
                    items[2]["CityObjects"]["unmatched"]["attributes"]["name"],
                    "keep"
                );
                items
            } else {
                vec![serde_json::from_slice(&fs::read(output).unwrap()).unwrap()]
            };
            let geometries = documents
                .iter()
                .flat_map(|document| document["CityObjects"].as_object().unwrap().values())
                .filter_map(|object| object["geometry"].as_array())
                .flatten()
                .collect::<Vec<_>>();
            assert_eq!(geometries.len(), 1);
            assert_eq!(geometries[0]["lod"], "2.2");
        }
        "obj" => assert_eq!(
            read_obj_lines(output)
                .iter()
                .filter(|line| line.starts_with("f "))
                .count(),
            12
        ),
        "tsv" => {
            let semantics =
                fs::read_to_string(output.with_file_name("selected_semantics.tsv")).unwrap();
            let mut lines = semantics.lines();
            let lod = lines
                .next()
                .unwrap()
                .split('\t')
                .position(|column| column == "geometry_lod")
                .unwrap();
            let rows = lines.collect::<Vec<_>>();
            assert!(!rows.is_empty());
            assert!(rows
                .iter()
                .all(|row| row.split('\t').nth(lod) == Some("2.2")));
        }
        "gpkg" => {
            let connection = rusqlite::Connection::open(output).unwrap();
            let table: String = connection
                .query_row("SELECT table_name FROM gpkg_geometry_columns", [], |row| {
                    row.get(0)
                })
                .unwrap();
            let lod: String = connection
                .query_row(
                    &format!("SELECT DISTINCT lod FROM \"{table}\""),
                    [],
                    |row| row.get(0),
                )
                .unwrap();
            assert_eq!(lod, "2.2");
            let count: u64 = connection
                .query_row(&format!("SELECT COUNT(*) FROM \"{table}\""), [], |row| {
                    row.get(0)
                })
                .unwrap();
            assert_eq!(count, 1);
        }
        _ => unreachable!(),
    }
}

#[test]
fn cjconvert_supports_legacy_feature_file_datasets() {
    let dataset = legacy_feature_dataset_path("cjconvert_legacy_feature_dataset");
    let cityjson_output = stable_output_path_with_suffix("cjconvert_legacy", "city.json");
    let cityjsonseq_output = stable_output_path_with_suffix("cjconvert_legacy", "city.jsonl");
    let obj_output = stable_output_path_with_suffix("cjconvert_legacy", "obj");
    let glb_output = stable_output_path("cjconvert_legacy");

    for (format, output) in [
        ("cityjson", &cityjson_output),
        ("cityjsonseq", &cityjsonseq_output),
        ("obj", &obj_output),
        ("glb", &glb_output),
    ] {
        let status = Command::new(env!("CARGO_BIN_EXE_cjconvert"))
            .arg(&dataset)
            .arg("--output")
            .arg(output)
            .arg("--format")
            .arg(format)
            .status()
            .unwrap_or_else(|error| panic!("run cjconvert {format}: {error}"));
        assert!(status.success(), "cjconvert {format} should succeed");
    }

    assert!(dataset.join(".cityjson-index.sqlite").is_file());

    let cityjson: Value = serde_json::from_slice(
        &fs::read(&cityjson_output).expect("read legacy dataset CityJSON output"),
    )
    .expect("parse legacy dataset CityJSON output");
    assert_eq!(cityjson["type"], "CityJSON");
    assert!(cityjson["CityObjects"].get("building-feature").is_some());
    assert!(cityjson["CityObjects"].get("water-feature").is_some());

    let sequence = read_json_lines(&cityjsonseq_output);
    assert_eq!(sequence.len(), 3);
    assert_eq!(sequence[0]["type"], "CityJSON");
    assert_eq!(sequence[1]["type"], "CityJSONFeature");
    assert_eq!(sequence[2]["type"], "CityJSONFeature");

    let obj = read_obj_lines(&obj_output);
    assert!(obj.iter().any(|line| line == "o building-feature"));
    assert!(obj.iter().any(|line| line == "o water-feature"));
    assert!(obj.iter().any(|line| line == "v 100 200 0"));

    let glb = fs::read(&glb_output).expect("read legacy dataset GLB output");
    let glb_json = read_glb_json(&glb);
    assert!(glb_json["meshes"]
        .as_array()
        .is_some_and(|meshes| !meshes.is_empty()));
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

fn read_integer_accessor(root: &Value, bin: &[u8], accessor_index: usize) -> Vec<u32> {
    let accessor = &root["accessors"][accessor_index];
    let view =
        &root["bufferViews"][usize::try_from(accessor["bufferView"].as_u64().unwrap()).unwrap()];
    let count = usize::try_from(accessor["count"].as_u64().unwrap()).unwrap();
    let component = accessor["componentType"].as_u64().unwrap();
    let size = match component {
        5121 => 1,
        5123 => 2,
        5125 => 4,
        _ => panic!("unexpected integer component"),
    };
    let offset = usize::try_from(accessor["byteOffset"].as_u64().unwrap_or(0)).unwrap();
    let decoded;
    let (bytes, stride) = if let Some(extension) = view["extensions"].get("EXT_meshopt_compression")
    {
        let start = usize::try_from(extension["byteOffset"].as_u64().unwrap()).unwrap();
        let length = usize::try_from(extension["byteLength"].as_u64().unwrap()).unwrap();
        let encoded = &bin[start..start + length];
        if extension["mode"] == "TRIANGLES" {
            return meshopt::decode_index_buffer::<u32>(encoded, count).unwrap();
        }
        assert_eq!(extension["byteStride"], 4);
        decoded = meshopt::decode_vertex_buffer::<[u8; 4]>(encoded, count)
            .unwrap()
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        (decoded.as_slice(), 4)
    } else {
        let start = usize::try_from(view["byteOffset"].as_u64().unwrap_or(0)).unwrap();
        let length = usize::try_from(view["byteLength"].as_u64().unwrap()).unwrap();
        let stride = view["byteStride"]
            .as_u64()
            .map_or(size, |value| usize::try_from(value).unwrap());
        (&bin[start..start + length], stride)
    };
    (0..count)
        .map(|index| {
            let start = offset + index * stride;
            match size {
                1 => u32::from(bytes[start]),
                2 => u32::from(u16::from_le_bytes(
                    bytes[start..start + 2].try_into().unwrap(),
                )),
                4 => u32::from_le_bytes(bytes[start..start + 4].try_into().unwrap()),
                _ => unreachable!(),
            }
        })
        .collect()
}

fn two_boxes_model(gap: f64) -> cityjson_lib::CityModel {
    let mut vertices = Vec::new();
    let mut objects = serde_json::Map::new();
    for (id, kind, offset) in [("a", "Building", 0.0), ("b", "BuildingPart", 10.0 + gap)] {
        let base = vertices.len();
        for point in [
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 10.0],
            [10.0, 0.0, 10.0],
            [10.0, 10.0, 10.0],
            [0.0, 10.0, 10.0],
        ] {
            vertices.push([point[0] + offset, point[1], point[2]]);
        }
        let boundaries = [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ]
        .map(|ring| vec![ring.map(|index| index + base)]);
        objects.insert(id.into(), serde_json::json!({"type":kind,"attributes":{"name":id},"geometry":[{"type":"MultiSurface","lod":"2.2","boundaries":boundaries}]}));
    }
    let input = serde_json::to_vec(&serde_json::json!({"type":"CityJSON","version":"2.0","transform":{"scale":[1,1,1],"translate":[0,0,0]},"CityObjects":objects,"vertices":vertices})).unwrap();
    json::from_slice(&input).unwrap()
}

fn assert_remeshed_glb(bytes: &[u8], quantize: bool, compress: bool, output_triangles: usize) {
    let root = read_glb_json(bytes);
    let bin = read_glb_bin(bytes);
    assert_eq!(
        root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["count"],
        2
    );
    let mut owners = std::collections::BTreeSet::new();
    let mut triangle_count = 0;
    for primitive in root["meshes"][0]["primitives"].as_array().unwrap() {
        let ids = read_integer_accessor(
            &root,
            bin,
            usize::try_from(primitive["attributes"]["_FEATURE_ID_0"].as_u64().unwrap()).unwrap(),
        );
        let indices = read_integer_accessor(
            &root,
            bin,
            usize::try_from(primitive["indices"].as_u64().unwrap()).unwrap(),
        );
        let feature_count = primitive["extensions"]["EXT_mesh_features"]["featureIds"][0]
            ["featureCount"]
            .as_u64()
            .unwrap();
        for triangle in indices.chunks_exact(3) {
            let owner = ids[usize::try_from(triangle[0]).unwrap()];
            assert!(u64::from(owner) < feature_count);
            assert!(owner < 2);
            assert!(triangle
                .iter()
                .all(|index| ids[usize::try_from(*index).unwrap()] == owner));
            owners.insert(owner);
            triangle_count += 1;
        }
        if !quantize && !compress {
            let normals = read_f32_vec3_accessor(
                &root,
                bin,
                usize::try_from(primitive["attributes"]["NORMAL"].as_u64().unwrap()).unwrap(),
            );
            assert!(normals.iter().flatten().all(|value| value.is_finite()));
            let position_accessor =
                usize::try_from(primitive["attributes"]["POSITION"].as_u64().unwrap()).unwrap();
            let positions = read_f32_vec3_accessor(&root, bin, position_accessor);
            assert!(positions.iter().flatten().all(|value| value.is_finite()));
            for position in &positions {
                for (axis, value) in position.iter().enumerate() {
                    let accessor = &root["accessors"][position_accessor];
                    assert!(f64::from(*value) >= accessor["min"][axis].as_f64().unwrap());
                    assert!(f64::from(*value) <= accessor["max"][axis].as_f64().unwrap());
                }
            }
        }
    }
    assert_eq!(owners.len(), 2);
    assert_eq!(triangle_count, output_triangles);
    assert!(root["materials"]
        .as_array()
        .unwrap()
        .iter()
        .all(|material| !material["doubleSided"].as_bool().unwrap_or(false)));
}

/// Purpose: verify end-to-end categorical transfer. Nearby and separated boxes
/// with different types/attributes must retain metadata, finite normals, valid
/// bounds, wide rectangular grids, and one owner per triangle in raw and
/// compressed/quantized GLBs.
#[test]
fn remeshed_glbs_transfer_one_feature_id_per_triangle() {
    for gap in [0.1, 25.0, 1000.0] {
        let model = two_boxes_model(gap);
        for (quantize, compress, smooth) in [
            (false, false, false),
            (true, false, true),
            (false, true, true),
            (true, true, false),
        ] {
            let directory = tempfile::tempdir().unwrap();
            let output = directory.path().join("proxy.glb");
            let report = convert_to_glb(
                &model,
                &output,
                &ExportOptions {
                    voxel_remesh: Some(VoxelRemeshOptions {
                        geometric_error: 2.0,
                    }),
                    quantize_geometry: quantize,
                    meshopt_compression: compress,
                    smooth_normals: smooth,
                    ..Default::default()
                },
            )
            .unwrap();
            let remesh = report.remeshing.unwrap();
            assert_eq!(remesh.input_triangles, 24);
            assert!(remesh.output_triangles > 0);
            assert!(remesh.remeshing_error_estimate <= remesh.requested_geometric_error / 2.0);
            if gap > 100.0 {
                assert!(remesh.grid_dimensions[0] > 1024);
                assert!(remesh.grid_dimensions[1] < 32);
                assert!(remesh.grid_dimensions[2] < 32);
            }
            assert!(remesh.estimated_geometric_error.is_finite());
            assert!(
                (remesh.estimated_geometric_error
                    - remesh.remeshing_error_estimate
                    - remesh.simplification_error
                    - remesh.encoding_error_estimate)
                    .abs()
                    < 1.0e-10
            );
            let bytes = fs::read(&output).unwrap();
            assert_eq!(report.output_bytes, u64::try_from(bytes.len()).unwrap());
            assert_remeshed_glb(&bytes, quantize, compress, remesh.output_triangles);
        }
    }
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
fn convert_to_glb_omits_empty_string_metadata_columns() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building-1":{
                    "type":"Building",
                    "attributes":{"name":""},
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
    let output_path = stable_output_path("empty-string-metadata");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed with empty string metadata");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);

    assert!(!has_extension(&root, "EXT_structural_metadata"));

    let buffer_views = root["bufferViews"]
        .as_array()
        .expect("glTF should contain bufferViews");
    for (index, buffer_view) in buffer_views.iter().enumerate() {
        assert!(
            buffer_view["byteLength"].as_u64().unwrap() > 0,
            "bufferView {index} should not have zero byteLength"
        );
    }
}

#[test]
fn convert_to_glb_sanitizes_structural_metadata_property_names() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building-1":{
                    "type":"Building",
                    "attributes":{
                        "3df_id":"building-1",
                        "3df_class":"Building"
                    },
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
    let output_path = stable_output_path("sanitized-metadata-properties");

    convert_to_glb(&model, &output_path, &ExportOptions::default())
        .expect("GLB conversion should succeed with numeric-leading metadata names");

    let glb_bytes = fs::read(&output_path).expect("test GLB should be written");
    let root = read_glb_json(&glb_bytes);
    let properties = root["extensions"]["EXT_structural_metadata"]["schema"]["classes"]
        ["cityobject"]["properties"]
        .as_object()
        .expect("metadata schema should contain sanitized properties");

    assert!(properties.contains_key("_3df_id"));
    assert!(properties.contains_key("_3df_class"));
    assert!(!properties.contains_key("3df_id"));
    assert!(!properties.contains_key("3df_class"));

    let table_properties = root["extensions"]["EXT_structural_metadata"]["propertyTables"][0]
        ["properties"]
        .as_object()
        .expect("property table should use the same sanitized names");
    assert!(table_properties.contains_key("_3df_id"));
    assert!(table_properties.contains_key("_3df_class"));
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

    let remesh_options = ExportOptions {
        voxel_remesh: Some(VoxelRemeshOptions {
            geometric_error: 1.0,
        }),
        ..options
    };
    assert!(convert_to_glb(&model, &output_path, &remesh_options)
        .unwrap_err()
        .to_string()
        .contains("clipped geometry"));
    assert_eq!(fs::read(&output_path).unwrap(), glb_bytes);
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

    let remesh_options = ExportOptions {
        voxel_remesh: Some(VoxelRemeshOptions {
            geometric_error: 1.0,
        }),
        ..options
    };
    assert!(convert_to_glb(&model, &output_path, &remesh_options)
        .unwrap_err()
        .to_string()
        .contains("clipped geometry"));
    assert_eq!(fs::read(&output_path).unwrap(), glb_bytes);
}

#[test]
fn convert_to_glb_drops_feature_types_fully_outside_clip_region() {
    // One TINRelief triangle sits inside the clip region; one Road triangle
    // sits entirely outside it. Before the empty-primitive fix, the Road
    // primitive survived the BTreeMap with zero vertices and the encoder
    // bailed with "primitive bounds missing for non-empty mesh", causing
    // the whole tile to be pruned.
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "relief-1":{
                    "type":"TINRelief",
                    "geometry":[
                        {
                            "type":"MultiSurface",
                            "lod":"1.0",
                            "boundaries":[
                                [[0,1,2]]
                            ]
                        }
                    ]
                },
                "road-1":{
                    "type":"Road",
                    "geometry":[
                        {
                            "type":"MultiSurface",
                            "lod":"1.0",
                            "boundaries":[
                                [[3,4,5]]
                            ]
                        }
                    ]
                }
            },
            "vertices":[
                [0.1,0.1,0.0],
                [0.9,0.1,0.0],
                [0.5,0.9,0.0],
                [5.0,5.0,0.0],
                [6.0,5.0,0.0],
                [5.5,6.0,0.0]
            ]
        }"#,
    )
    .expect("inline CityJSON fixture should parse");
    let output_path = stable_output_path("clip-drops-empty-feature-type");
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
        .expect("GLB conversion should succeed when one feature_type is fully clipped away");

    let glb_bytes = fs::read(&output_path).expect("test GLB should exist");
    let root = read_glb_json(&glb_bytes);
    let primitives = root["meshes"][0]["primitives"]
        .as_array()
        .expect("mesh should expose primitives");
    assert_eq!(
        primitives.len(),
        1,
        "primitive for fully-clipped Road should be dropped, leaving only TINRelief"
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

    for (index, primitive) in primitives.iter().enumerate() {
        assert!(primitive["attributes"].get("_FEATURE_ID_0").is_some());
        assert_eq!(
            primitive["extensions"]["EXT_mesh_features"]["featureIds"][0]["featureCount"]
                .as_u64()
                .unwrap(),
            u64::try_from(index + 1).unwrap()
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
