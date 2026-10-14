#![cfg(any())]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use serde_json::Value;

fn unique_test_dir(prefix: &str) -> PathBuf {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time")
        .as_nanos();
    let path = std::env::temp_dir().join(format!("tyler-{prefix}-{unique}"));
    fs::create_dir_all(&path).expect("create test dir");
    path
}

fn repo_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn read_fixture(relative: &str) -> String {
    fs::read_to_string(repo_path(relative)).expect("read fixture")
}

fn write_ndjson_dataset(prefix: &str, metadata: &str, feature_blobs: &[String]) -> PathBuf {
    let dataset = unique_test_dir(prefix);
    let mut contents = String::new();
    contents.push_str(metadata.trim_end());
    contents.push('\n');
    for feature_blob in feature_blobs {
        contents.push_str(feature_blob.trim_end());
        contents.push('\n');
    }
    fs::write(dataset.join("source.city.jsonl"), contents).expect("write ndjson source");
    dataset
}

fn run_tyler(dataset: &Path, output: &Path, args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_tyler"));
    command.arg(dataset).arg("--output").arg(output);
    for arg in args {
        command.arg(arg);
    }
    command.output().expect("run tyler")
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read json file")).expect("parse json file")
}

fn read_glb_json(bytes: &[u8]) -> Value {
    assert!(bytes.len() >= 20, "glb should contain a header");
    assert_eq!(&bytes[0..4], b"glTF");
    assert_eq!(u32::from_le_bytes(bytes[4..8].try_into().unwrap()), 2);

    let declared_length = u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize;
    assert_eq!(declared_length, bytes.len());

    let json_length = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
    assert_eq!(&bytes[16..20], b"JSON");

    serde_json::from_slice(&bytes[20..20 + json_length]).expect("GLB JSON chunk should parse")
}

fn find_first_glb(dir: &Path) -> Option<PathBuf> {
    for entry in fs::read_dir(dir).ok()? {
        let entry = entry.ok()?;
        let path = entry.path();
        if path.is_dir() {
            if let Some(found) = find_first_glb(&path) {
                return Some(found);
            }
        } else if path.extension().is_some_and(|extension| extension == "glb") {
            return Some(path);
        }
    }
    None
}

fn collect_paths_with_extension(dir: &Path, extension: &str, out: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                collect_paths_with_extension(&path, extension, out);
            } else if path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with(&format!(".{extension}")))
            {
                out.push(path);
            }
        }
    }
}

fn region_from_bounding_volume(bounding_volume: &Value) -> Vec<f64> {
    bounding_volume["region"]
        .as_array()
        .expect("region should be an array")
        .iter()
        .map(|value| value.as_f64().expect("region component should be numeric"))
        .collect()
}

fn first_glb_json(output_dir: &Path) -> Value {
    let glb_path = find_first_glb(&output_dir.join("t")).expect("expected at least one GLB");
    read_glb_json(&fs::read(glb_path).expect("read glb"))
}

fn material_count(glb_json: &Value) -> usize {
    glb_json["materials"]
        .as_array()
        .expect("glTF materials should exist")
        .len()
}

fn metadata_feature_count(glb_json: &Value) -> u64 {
    glb_json["extensions"]["EXT_structural_metadata"]["propertyTables"][0]["count"]
        .as_u64()
        .expect("structural metadata feature count should exist")
}

fn first_tile_with_content(tile: &Value) -> Option<&Value> {
    if let Some(children) = tile.get("children").and_then(Value::as_array) {
        return children.iter().find_map(first_tile_with_content);
    }
    tile.get("content").is_some().then_some(tile)
}

fn leaf_content_tile(tile: &Value) -> &Value {
    first_tile_with_content(tile).expect("expected at least one tile with content")
}

fn assert_success(output: &Output, context: &str) {
    assert!(
        output.status.success(),
        "{context} failed: status={:?}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn object_type_filters_mixed_feature_glb_and_accepts_union() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let mixed_feature = serde_json::json!({
        "type": "CityJSONFeature",
        "id": "mixed-feature",
        "CityObjects": {
            "building-object": {
                "type": "Building",
                "attributes": {"name": "Building A"},
                "geometry": [{
                    "type": "MultiSurface",
                    "lod": "1",
                    "boundaries": [[[0, 1, 2], [0, 2, 3]]]
                }]
            },
            "water-object": {
                "type": "WaterBody",
                "attributes": {"name": "Canal Segment"},
                "geometry": [{
                    "type": "MultiSurface",
                    "lod": "1",
                    "boundaries": [[[4, 5, 6], [4, 6, 7]]]
                }]
            }
        },
        "vertices": [
            [0, 0, 0],
            [4, 0, 0],
            [4, 4, 0],
            [0, 4, 0],
            [10, 0, 0],
            [13, 0, 0],
            [13, 3, 0],
            [10, 3, 0]
        ]
    })
    .to_string();
    let dataset = write_ndjson_dataset("object-type", &metadata, &[mixed_feature]);

    let building_output = unique_test_dir("object-type-building");
    let output = run_tyler(&dataset, &building_output, &["--object-type", "Building"]);
    assert_success(&output, "building-only object type run");

    let building_glb = first_glb_json(&building_output);
    assert_eq!(material_count(&building_glb), 1);
    assert_eq!(metadata_feature_count(&building_glb), 1);

    let union_output = unique_test_dir("object-type-union");
    let output = run_tyler(
        &dataset,
        &union_output,
        &["--object-type", "Building", "--object-type", "WaterBody"],
    );
    assert_success(&output, "union object type run");

    let union_glb = first_glb_json(&union_output);
    assert_eq!(material_count(&union_glb), 2);
    assert_eq!(metadata_feature_count(&union_glb), 2);
}

#[test]
fn object_type_missing_type_fails() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("object-type-missing", &metadata, &[features]);
    let output_dir = unique_test_dir("object-type-missing-output");

    let output = run_tyler(&dataset, &output_dir, &["--object-type", "Road"]);

    assert!(
        !output.status.success(),
        "expected missing object type selection to fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Did not find any CityJSONFeatures of type"),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn lod_pruning_reduces_generated_glb_content() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_lod_building_part.city.jsonl");
    let dataset = write_ndjson_dataset("lod-pruning", &metadata, &[features]);

    let control_output = unique_test_dir("lod-pruning-control");
    let output = run_tyler(
        &dataset,
        &control_output,
        &["--object-type", "BuildingPart"],
    );
    assert_success(&output, "lod pruning control run");

    let control_glb = read_glb_json(
        &fs::read(
            find_first_glb(&control_output.join("t"))
                .expect("expected at least one GLB for the control run"),
        )
        .expect("read control glb"),
    );
    let control_primitive_count = control_glb["meshes"][0]["primitives"]
        .as_array()
        .expect("control primitives should exist")
        .len();
    let control_byte_length = control_glb["buffers"][0]["byteLength"]
        .as_u64()
        .expect("control buffer length should exist");

    let pruned_output = unique_test_dir("lod-pruning-pruned");
    let output = run_tyler(
        &dataset,
        &pruned_output,
        &[
            "--object-type",
            "BuildingPart",
            "--lod-building-part",
            "2.2",
        ],
    );
    assert_success(&output, "lod pruning run");

    let pruned_glb = read_glb_json(
        &fs::read(
            find_first_glb(&pruned_output.join("t"))
                .expect("expected at least one GLB for the pruned run"),
        )
        .expect("read pruned glb"),
    );
    let pruned_primitive_count = pruned_glb["meshes"][0]["primitives"]
        .as_array()
        .expect("pruned primitives should exist")
        .len();
    let pruned_byte_length = pruned_glb["buffers"][0]["byteLength"]
        .as_u64()
        .expect("pruned buffer length should exist");

    assert!(pruned_byte_length < control_byte_length);
    assert!(pruned_primitive_count <= control_primitive_count);
}

#[test]
fn include_parent_attributes_copies_parent_values_into_selected_child_objects() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = serde_json::json!({
        "type": "CityJSONFeature",
        "id": "building-parent",
        "CityObjects": {
            "building-parent": {
                "type": "Building",
                "attributes": {
                    "parent_only": "parent",
                    "shared": "parent",
                    "levels": 7
                },
                "children": ["building-part"]
            },
            "building-part": {
                "type": "BuildingPart",
                "parents": ["building-parent"],
                "attributes": {
                    "child_only": "child",
                    "shared": "child"
                },
                "geometry": [{
                    "type": "MultiSurface",
                    "lod": "1",
                    "boundaries": [[[0, 1, 2], [0, 2, 3]]]
                }]
            }
        },
        "vertices": [
            [0, 0, 0],
            [4, 0, 0],
            [4, 4, 0],
            [0, 4, 0]
        ]
    })
    .to_string();
    let dataset = write_ndjson_dataset("parent-attributes", &metadata, &[feature]);

    let control_output = unique_test_dir("parent-attributes-control");
    let output = run_tyler(
        &dataset,
        &control_output,
        &[
            "--object-type",
            "BuildingPart",
            "--3dtiles-metadata-class",
            "parcel",
        ],
    );
    assert_success(&output, "parent attribute control run");

    let control_glb = read_glb_json(
        &fs::read(
            find_first_glb(&control_output.join("t"))
                .expect("expected at least one GLB for the control run"),
        )
        .expect("read control glb"),
    );
    let control_properties = control_glb["extensions"]["EXT_structural_metadata"]["schema"]
        ["classes"]["parcel"]["properties"]
        .as_object()
        .expect("structural metadata properties should exist");
    assert!(!control_properties.contains_key("parent_only"));

    let inherited_output = unique_test_dir("parent-attributes-inherited");
    let output = run_tyler(
        &dataset,
        &inherited_output,
        &[
            "--object-type",
            "BuildingPart",
            "--include-parent-attributes",
            "--3dtiles-metadata-class",
            "parcel",
        ],
    );
    assert_success(&output, "parent attribute inheritance run");

    let inherited_glb = read_glb_json(
        &fs::read(
            find_first_glb(&inherited_output.join("t"))
                .expect("expected at least one GLB for the inheritance run"),
        )
        .expect("read inherited glb"),
    );
    let inherited_properties = inherited_glb["extensions"]["EXT_structural_metadata"]["schema"]
        ["classes"]["parcel"]["properties"]
        .as_object()
        .expect("structural metadata properties should exist");
    assert!(inherited_properties.contains_key("parent_only"));
    assert!(inherited_properties.contains_key("child_only"));
    assert!(inherited_properties.contains_key("shared"));
}

#[test]
fn tileset_flags_change_output_artifacts_and_numeric_fields() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("tileset-flags", &metadata, &[features]);

    let default_output = unique_test_dir("tileset-default");
    let output = run_tyler(&dataset, &default_output, &["--qtree-capacity", "1"]);
    assert_success(&output, "default tileset run");

    let default_tileset = read_json(&default_output.join("tileset.json"));
    let default_root = &default_tileset["root"];
    let default_geometric_error = default_root["geometricError"]
        .as_f64()
        .expect("default geometricError should be numeric");
    let default_region = region_from_bounding_volume(&default_root["boundingVolume"]);
    assert!(default_root["content"]["boundingVolume"].is_null());

    let with_content_bv = unique_test_dir("tileset-content-bv");
    let output = run_tyler(
        &dataset,
        &with_content_bv,
        &["--3dtiles-content-add-bv", "--qtree-capacity", "1"],
    );
    assert_success(&output, "content bounding volume run");

    let with_content_bv_tileset = read_json(&with_content_bv.join("tileset.json"));
    let with_content_bv_root = &with_content_bv_tileset["root"];
    let content_tile = leaf_content_tile(with_content_bv_root);
    assert!(content_tile["content"]["boundingVolume"].is_object());

    let content_region = region_from_bounding_volume(&content_tile["content"]["boundingVolume"]);
    let tile_region = region_from_bounding_volume(&content_tile["boundingVolume"]);
    assert!(tile_region[0] <= content_region[0]);
    assert!(tile_region[1] <= content_region[1]);
    assert!(tile_region[2] >= content_region[2]);
    assert!(tile_region[3] >= content_region[3]);
    assert!(tile_region[4] <= content_region[4]);
    assert!(tile_region[5] >= content_region[5]);

    let modified_output = unique_test_dir("tileset-modified");
    let output = run_tyler(
        &dataset,
        &modified_output,
        &[
            "--3dtiles-content-add-bv",
            "--3dtiles-content-bv-from-tile",
            "--qtree-capacity",
            "1",
            "--grid-minz=-10",
            "--grid-maxz=10",
        ],
    );
    assert_success(&output, "modified tileset run");

    let modified_tileset = read_json(&modified_output.join("tileset.json"));
    let modified_root = &modified_tileset["root"];
    let modified_geometric_error = modified_root["geometricError"]
        .as_f64()
        .expect("modified geometricError should be numeric");
    let modified_region = region_from_bounding_volume(&modified_root["boundingVolume"]);
    assert!(
        (default_geometric_error - modified_geometric_error).abs() <= f64::EPSILON,
        "grid height and content bounding-volume settings should not change traversal error"
    );
    assert!(modified_region[4] <= default_region[4]);
    assert!(modified_region[5] >= default_region[5]);

    let modified_content_tile = leaf_content_tile(modified_root);
    let modified_content_region =
        region_from_bounding_volume(&modified_content_tile["content"]["boundingVolume"]);
    let modified_tile_region =
        region_from_bounding_volume(&modified_content_tile["boundingVolume"]);
    assert_eq!(
        modified_content_region, modified_tile_region,
        "with --3dtiles-content-bv-from-tile the content BV should match the tile BV"
    );
}

#[test]
fn metadata_class_and_color_flags_propagate_into_glb_json() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("glb-metadata", &metadata, &[features]);
    let output_dir = unique_test_dir("glb-metadata-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--object-type",
            "Building",
            "--3dtiles-metadata-class",
            "parcel",
            "--color-building",
            "#010203",
        ],
    );
    assert_success(&output, "metadata/color run");

    let glb_path = find_first_glb(&output_dir.join("t")).expect("expected at least one GLB");
    let glb_json = read_glb_json(&fs::read(glb_path).expect("read glb"));

    let extensions_used = glb_json["extensionsUsed"]
        .as_array()
        .expect("glTF should declare extensionsUsed");
    assert!(extensions_used
        .iter()
        .any(|value| value.as_str() == Some("EXT_structural_metadata")));

    let structural_metadata = glb_json["extensions"]["EXT_structural_metadata"]
        .as_object()
        .expect("structural metadata extension should exist");
    let schema = structural_metadata["schema"]
        .as_object()
        .expect("structural metadata schema should exist");
    let classes = schema["classes"]
        .as_object()
        .expect("structural metadata classes should exist");
    assert!(classes.contains_key("parcel"));
    assert_eq!(classes["parcel"]["name"].as_str(), Some("parcel"));
    assert_eq!(
        structural_metadata["propertyTables"][0]["class"].as_str(),
        Some("parcel")
    );

    let material = &glb_json["materials"][0];
    let base_color = material["pbrMetallicRoughness"]["baseColorFactor"]
        .as_array()
        .expect("baseColorFactor should be present");
    let expected = [1.0 / 255.0, 2.0 / 255.0, 3.0 / 255.0, 1.0];
    for (component, expected) in base_color.iter().zip(expected) {
        let actual = component
            .as_f64()
            .expect("baseColorFactor component should be numeric");
        assert!(
            (actual - expected).abs() < 1.0e-6,
            "expected {expected}, got {actual}"
        );
    }
}

#[test]
fn implicit_tiling_writes_subtrees_and_marks_tileset_as_implicit() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("implicit-tiling", &metadata, &[features]);
    let output_dir = unique_test_dir("implicit-tiling-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--object-type",
            "Building",
            "--3dtiles-implicit",
            "--3dtiles-metadata-class",
            "cityobject",
        ],
    );
    assert_success(&output, "implicit tiling run");

    let tileset = read_json(&output_dir.join("tileset.json"));
    assert!(tileset["root"]["implicitTiling"].is_object());

    let mut subtree_files = Vec::new();
    collect_paths_with_extension(&output_dir.join("subtrees"), "subtree", &mut subtree_files);
    assert!(
        !subtree_files.is_empty(),
        "expected implicit tiling subtree outputs under {}",
        output_dir.join("subtrees").display()
    );
}

#[test]
fn tileset_only_skips_glb_tile_directory() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("tileset-only", &metadata, &[features]);
    let output_dir = unique_test_dir("tileset-only-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--3dtiles-tileset-only",
            "--3dtiles-metadata-class",
            "cityobject",
        ],
    );
    assert_success(&output, "tileset-only run");

    assert!(output_dir.join("tileset.json").is_file());
    assert!(
        !output_dir.join("t").exists(),
        "tileset-only output should not create GLB tile directories"
    );
}
