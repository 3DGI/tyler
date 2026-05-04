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

fn read_json_line_type(line: &str) -> String {
    let value: Value = serde_json::from_str(line).expect("parse json line");
    value["type"]
        .as_str()
        .expect("json line should have a string type")
        .to_string()
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

fn assert_success(output: &Output, context: &str) {
    assert!(
        output.status.success(),
        "{context} failed: status={:?}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn collect_paths_with_suffix(dir: &Path, suffix: &str, out: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                collect_paths_with_suffix(&path, suffix, out);
            } else if path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with(suffix))
            {
                out.push(path);
            }
        }
    }
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

fn zero_grid_vertex_counts(path: &Path) {
    let rewritten = fs::read_to_string(path)
        .expect("read grid tsv")
        .lines()
        .enumerate()
        .map(|(index, line)| {
            if index <= 1 || line.trim().is_empty() {
                return line.to_string();
            }
            let mut parts = line.splitn(3, '\t');
            let cell_id = parts.next().expect("cell id");
            let _nr_items = parts.next().expect("nr_items");
            let wkt = parts.next().expect("wkt");
            format!("{cell_id}\t0\t{wkt}")
        })
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(path, format!("{rewritten}\n")).expect("rewrite grid tsv");
}

#[test]
fn debug_dump_data_writes_bincode_and_intermediary_cityjson() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("debug-dump-data", &metadata, &[feature]);
    let output_dir = unique_test_dir("debug-dump-data-output");

    let output = run_tyler(&dataset, &output_dir, &["--debug-dump-data"]);
    assert_success(&output, "debug dump data run");

    let debug_dir = output_dir.join("debug");
    assert!(debug_dir.join("world.bincode").is_file());
    assert!(debug_dir.join("quadtree.bincode").is_file());
    assert!(debug_dir.join("tiles_results.bincode").is_file());

    let mut cityjson_inputs = Vec::new();
    collect_paths_with_suffix(
        &debug_dir.join("inputs"),
        ".city.jsonl",
        &mut cityjson_inputs,
    );
    assert!(
        !cityjson_inputs.is_empty(),
        "expected intermediary CityJSON dumps under {}",
        debug_dir.join("inputs").display()
    );
}

#[test]
fn format_cityjson_writes_cityjson_tiles() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("format-cityjson", &metadata, &[feature]);
    let output_dir = unique_test_dir("format-cityjson-output");

    let output = run_tyler(&dataset, &output_dir, &["--format", "cityjson"]);
    assert_success(&output, "CityJSON format run");

    let mut cityjson_tiles = Vec::new();
    collect_paths_with_suffix(&output_dir.join("t"), ".city.json", &mut cityjson_tiles);
    assert!(
        !cityjson_tiles.is_empty(),
        "expected CityJSON tiles under {}",
        output_dir.join("t").display()
    );
    let root = read_json(&cityjson_tiles[0]);
    assert_eq!(root["type"], "CityJSON");
    assert!(!output_dir.join("tileset.json").exists());
}

#[test]
fn format_cityjsonseq_writes_feature_stream_tiles_matching_debug_shape() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("format-cityjsonseq", &metadata, &[feature]);
    let output_dir = unique_test_dir("format-cityjsonseq-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &["--format", "cityjsonseq", "--debug-dump-data"],
    );
    assert_success(&output, "CityJSONSeq format run");

    let mut output_streams = Vec::new();
    collect_paths_with_suffix(&output_dir.join("t"), ".city.jsonl", &mut output_streams);
    assert!(
        !output_streams.is_empty(),
        "expected CityJSONSeq tiles under {}",
        output_dir.join("t").display()
    );
    let mut debug_streams = Vec::new();
    collect_paths_with_suffix(
        &output_dir.join("debug").join("inputs"),
        ".city.jsonl",
        &mut debug_streams,
    );
    assert!(
        !debug_streams.is_empty(),
        "expected debug CityJSONSeq streams under {}",
        output_dir.join("debug").join("inputs").display()
    );

    let output_items = fs::read_to_string(&output_streams[0]).expect("read output stream");
    let debug_items = fs::read_to_string(&debug_streams[0]).expect("read debug stream");
    let output_types = output_items
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(read_json_line_type)
        .collect::<Vec<_>>();
    let debug_types = debug_items
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(read_json_line_type)
        .collect::<Vec<_>>();
    assert_eq!(output_types, debug_types);
    assert_eq!(output_types.first().map(String::as_str), Some("CityJSON"));
    assert!(output_types.iter().any(|kind| kind == "CityJSONFeature"));
    assert!(!output_dir.join("tileset.json").exists());
}

#[test]
fn debug_dump_grid_and_grid_features_write_tsv_exports() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("debug-dump-grid", &metadata, &[feature]);
    let output_dir = unique_test_dir("debug-dump-grid-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--debug-dump-grid",
            "--debug-dump-grid-features",
            "--qtree-capacity",
            "1",
        ],
    );
    assert_success(&output, "debug dump grid run");

    let debug_dir = output_dir.join("debug");
    assert!(debug_dir.join("grid.tsv").is_file());
    assert!(debug_dir.join("features.tsv").is_file());
    assert!(debug_dir.join("quadtree_level-0.tsv").is_file());
    assert!(debug_dir.join("tileset_level-0.tsv").is_file());
}

#[test]
fn debug_load_grid_uses_loaded_grid_for_quadtree_computation() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("debug-load-grid", &metadata, &[feature]);

    let seeded_output = unique_test_dir("debug-load-grid-seeded");
    let output = run_tyler(
        &dataset,
        &seeded_output,
        &[
            "--debug-dump-grid",
            "--debug-dump-grid-features",
            "--qtree-capacity",
            "1",
        ],
    );
    assert_success(&output, "seed grid dump run");

    let seeded_tileset = read_json(&seeded_output.join("tileset.json"));
    assert!(
        seeded_tileset["root"]["children"].is_array(),
        "baseline run should split the quadtree"
    );

    let replay_grid_dir = unique_test_dir("debug-load-grid-replay");
    let replay_grid = replay_grid_dir.join("grid.tsv");
    fs::copy(seeded_output.join("debug").join("grid.tsv"), &replay_grid).expect("copy grid");
    fs::copy(
        seeded_output.join("debug").join("features.tsv"),
        replay_grid_dir.join("features.tsv"),
    )
    .expect("copy features");
    zero_grid_vertex_counts(&replay_grid);

    let replay_output = unique_test_dir("debug-load-grid-output");
    let output = run_tyler(
        &dataset,
        &replay_output,
        &[
            "--debug-load-grid",
            replay_grid.to_str().expect("utf8 replay grid path"),
            "--debug-dump-grid",
            "--qtree-capacity",
            "1",
        ],
    );
    assert_success(&output, "debug load grid replay run");

    let replay_tileset = read_json(&replay_output.join("tileset.json"));
    assert!(
        replay_tileset["root"].get("children").is_none(),
        "loaded zero-count grid should collapse the quadtree to the root tile"
    );

    let mut quadtree_levels = Vec::new();
    collect_paths_with_suffix(&replay_output.join("debug"), ".tsv", &mut quadtree_levels);
    let quadtree_level_files = quadtree_levels
        .into_iter()
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with("quadtree_level-"))
        })
        .collect::<Vec<_>>();
    assert_eq!(
        quadtree_level_files.len(),
        1,
        "expected the loaded grid to produce a single-level quadtree"
    );
}

#[test]
fn debug_3dtiles_tileset_only_skips_glb_conversion() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let features = read_fixture("cityjson-convert/tests/data/multi_feature_types.city.jsonl");
    let dataset = write_ndjson_dataset("debug-tileset-only", &metadata, &[features]);
    let output_dir = unique_test_dir("debug-tileset-only-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &["--debug-3dtiles-tileset-only", "--qtree-capacity", "1"],
    );
    assert_success(&output, "tileset-only run");

    assert!(output_dir.join("tileset.json").is_file());
    assert!(
        !output_dir.join("t").exists(),
        "tileset-only mode should skip GLB tile output"
    );
}

#[test]
fn object_attributes_filter_and_type_glb_metadata_schema() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = serde_json::json!({
        "type": "CityJSONFeature",
        "id": "attribute-types",
        "transform": {
            "scale": [1.0, 1.0, 1.0],
            "translate": [0.0, 0.0, 0.0]
        },
        "CityObjects": {
            "building": {
                "type": "Building",
                "attributes": {
                    "as_text": 7,
                    "as_bool": "true",
                    "as_int": 9.0,
                    "as_float": 3,
                    "ignored": "drop-me"
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
    let dataset = write_ndjson_dataset("object-attributes", &metadata, &[feature]);
    let output_dir = unique_test_dir("object-attributes-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--object-type",
            "Building",
            "--object-attributes",
            "as_text:string,as_bool:bool,as_int:int,as_float:float",
        ],
    );
    assert_success(&output, "object attributes run");

    let glb_path = find_first_glb(&output_dir.join("t")).expect("expected at least one GLB");
    let glb_json = read_glb_json(&fs::read(glb_path).expect("read glb"));
    let properties = glb_json["extensions"]["EXT_structural_metadata"]["schema"]["classes"]
        ["citymodel"]["properties"]
        .as_object()
        .expect("structural metadata schema should exist");

    assert_eq!(properties.len(), 4);
    assert_eq!(properties["as_text"]["type"].as_str(), Some("STRING"));
    assert_eq!(properties["as_bool"]["type"].as_str(), Some("SCALAR"));
    assert_eq!(
        properties["as_bool"]["componentType"].as_str(),
        Some("INT8")
    );
    assert_eq!(
        properties["as_int"]["componentType"].as_str(),
        Some("INT32")
    );
    assert_eq!(
        properties["as_float"]["componentType"].as_str(),
        Some("FLOAT32")
    );
    assert!(!properties.contains_key("ignored"));
}

#[test]
fn object_type_building_single_tile_tileset_keeps_positive_root_geometric_error() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = serde_json::json!({
        "type": "CityJSONFeature",
        "id": "single-building",
        "transform": {
            "scale": [1.0, 1.0, 1.0],
            "translate": [0.0, 0.0, 0.0]
        },
        "CityObjects": {
            "building": {
                "type": "Building",
                "attributes": {
                    "name": "Single building"
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
    let dataset = write_ndjson_dataset("object-type-building-single", &metadata, &[feature]);
    let output_dir = unique_test_dir("object-type-building-single-output");

    let output = run_tyler(&dataset, &output_dir, &["--object-type", "Building"]);
    assert_success(&output, "single building tileset run");

    let tileset = read_json(&output_dir.join("tileset.json"));
    let root = &tileset["root"];
    assert!(root["content"].is_object(), "root tile should keep content");
    assert!(
        root.get("children").is_none(),
        "single-building dataset should produce a single root tile"
    );
    let root_geometric_error = root["geometricError"]
        .as_f64()
        .expect("root geometricError should be numeric");
    assert!(
        root_geometric_error.abs() <= f64::EPSILON,
        "root tile geometricError can remain zero"
    );
    let tileset_geometric_error = tileset["geometricError"]
        .as_f64()
        .expect("tileset geometricError should be numeric");
    assert!(
        tileset_geometric_error > 0.0,
        "tileset geometricError should stay positive for a single-tile tileset"
    );
    assert!(tileset_geometric_error > f64::EPSILON);

    let unpruned_tileset = read_json(&output_dir.join("tileset_unpruned.json"));
    let unpruned_root = &unpruned_tileset["root"];
    assert!(unpruned_root["content"].is_object());
    let unpruned_root_geometric_error = unpruned_root["geometricError"]
        .as_f64()
        .expect("unpruned root geometricError should be numeric");
    assert!(
        unpruned_root_geometric_error.abs() <= f64::EPSILON,
        "unpruned root geometricError can remain zero"
    );
    let unpruned_tileset_geometric_error = unpruned_tileset["geometricError"]
        .as_f64()
        .expect("unpruned tileset geometricError should be numeric");
    assert!(
        unpruned_tileset_geometric_error > 0.0,
        "unpruned tileset geometricError should stay positive"
    );
    assert!(unpruned_tileset_geometric_error > f64::EPSILON);
}
