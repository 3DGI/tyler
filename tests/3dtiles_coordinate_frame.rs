use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
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

fn assert_tile_bounding_volumes_are_regions(tile: &Value) {
    assert!(
        tile["boundingVolume"]["region"].is_array(),
        "tile boundingVolume should be a region: {tile}"
    );

    if let Some(content) = tile.get("content") {
        if let Some(bounding_volume) = content.get("boundingVolume") {
            assert!(
                bounding_volume["region"].is_array(),
                "content boundingVolume should be a region: {content}"
            );
        }
    }

    if let Some(children) = tile.get("children").and_then(Value::as_array) {
        for child in children {
            assert_tile_bounding_volumes_are_regions(child);
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

fn assert_local_enu_tileset(root: &Value) {
    let root_transform = root["transform"]
        .as_array()
        .expect("root tile should have an ENU-to-ECEF transform");
    assert_eq!(root_transform.len(), 16);
    assert!((root_transform[0].as_f64().unwrap() - 1.0).abs() > 1.0e-6);
    assert!(root_transform[12].as_f64().unwrap().abs() > 1.0);

    assert_tile_bounding_volumes_are_regions(root);
    assert_eq!(
        root["refine"].as_str(),
        Some("ADD"),
        "implicit root should use additive refinement so parent content remains visible"
    );
}

fn assert_local_enu_glb(glb_path: &Path) {
    let glb_json = read_glb_json(&fs::read(glb_path).expect("GLB should be readable"));
    let node_matrix = glb_json["nodes"][0]["matrix"]
        .as_array()
        .expect("GLB root node should carry a matrix");
    let translation = [
        node_matrix[12].as_f64().unwrap(),
        node_matrix[13].as_f64().unwrap(),
        node_matrix[14].as_f64().unwrap(),
    ];
    assert!(
        translation
            .iter()
            .all(|component| component.abs() < 100_000.0),
        "GLB node translation should be local-scale, got {translation:?}"
    );
}

#[test]
fn debug_replay_writes_epsg4979_regions_and_local_enu_glbs() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let dataset = unique_test_dir("coordinate-frame-dataset");
    let metadata =
        fs::read_to_string(repo.join("resources/data/3dbag_x00.city.json")).expect("read metadata");
    let feature = fs::read_to_string(repo.join("resources/data/3dbag_feature_x71.city.jsonl"))
        .expect("read feature");
    fs::write(
        dataset.join("source.city.jsonl"),
        format!("{metadata}\n{feature}\n"),
    )
    .expect("write ndjson source");

    let seeded_output = unique_test_dir("coordinate-frame-seeded");
    let output = unique_test_dir("coordinate-frame");

    let seed_status = Command::new(env!("CARGO_BIN_EXE_tyler"))
        .env("RUST_LOG", "debug")
        .arg(&dataset)
        .arg("--output")
        .arg(&seeded_output)
        .arg("--3dtiles-implicit")
        .arg("--3dtiles-metadata-class")
        .arg("building")
        .status()
        .expect("tyler should seed debug data");
    assert!(
        seed_status.success(),
        "tyler debug seeding failed: {seed_status}"
    );

    let debug_data = seeded_output.join("debug");
    assert!(
        debug_data.is_dir(),
        "seed run should create debug bincode data"
    );

    let replay_data = unique_test_dir("coordinate-frame-replay");
    for filename in ["world.bincode", "quadtree.bincode"] {
        fs::copy(debug_data.join(filename), replay_data.join(filename)).expect("copy replay data");
    }
    let mut previous: Vec<(PathBuf, Vec<u8>)> = Vec::new();
    for suffix in ["lod2.2", "LoD_3-v1"] {
        let status = Command::new(env!("CARGO_BIN_EXE_tyler"))
            .arg("--output-suffix")
            .arg(suffix)
            .arg(&dataset)
            .arg("--output")
            .arg(&output)
            .arg("--debug-load-data")
            .arg(&replay_data)
            .arg("--3dtiles-implicit")
            .arg("--3dtiles-metadata-class")
            .arg("building")
            .status()
            .expect("tyler should run");
        assert!(status.success(), "tyler debug replay failed: {status}");

        let tileset_path = output.join(format!("tileset-{suffix}.json"));
        let tileset: Value = serde_json::from_slice(
            &fs::read(&tileset_path).expect("tileset.json should be generated"),
        )
        .expect("tileset.json should parse");

        let root = &tileset["root"];
        assert_local_enu_tileset(root);
        assert_eq!(
            root["content"]["uri"].as_str(),
            Some(format!("t/{{level}}/{{x}}/{{level}}-{{x}}-{{y}}-{suffix}.glb").as_str())
        );

        let glb_path = output.join(format!("t/0/0/0-0-0-{suffix}.glb"));
        assert!(find_first_glb(&output.join("t")).is_some());
        assert_eq!(
            glb_path.strip_prefix(&output).expect("GLB under output"),
            Path::new(&format!("t/0/0/0-0-0-{suffix}.glb"))
        );
        assert_local_enu_glb(&glb_path);
        let subtree_template = root["implicitTiling"]["subtrees"]["uri"]
            .as_str()
            .expect("subtree template");
        assert_eq!(
            subtree_template,
            format!("subtrees/{{level}}/{{x}}/{{y}}-{suffix}.subtree")
        );
        let subtree_path = output.join(
            subtree_template
                .replace("{level}", "0")
                .replace("{x}", "0")
                .replace("{y}", "0"),
        );
        assert!(fs::read(&subtree_path)
            .expect("read subtree")
            .starts_with(b"subt"));
        for (path, bytes) in &previous {
            assert_eq!(&fs::read(path).expect("previous export exists"), bytes);
        }
        for path in [tileset_path, glb_path, subtree_path] {
            let bytes = fs::read(&path).expect("read published file");
            previous.push((path, bytes));
        }
    }
}
