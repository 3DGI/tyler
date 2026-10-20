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
    let result = command.output().expect("run tyler");

    // Re-emit subprocess stdout/stderr to test process
    // This makes it controlled by cargo test's --show-output / --nocapture
    if !result.stdout.is_empty() {
        println!("{}", String::from_utf8_lossy(&result.stdout));
    }
    if !result.stderr.is_empty() {
        eprintln!("{}", String::from_utf8_lossy(&result.stderr));
    }

    result
}

fn run_tyler_with_rust_log(dataset: &Path, output: &Path, rust_log: &str, args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_tyler"));
    command.env("RUST_LOG", rust_log);
    command.arg(dataset).arg("--output").arg(output);
    for arg in args {
        command.arg(arg);
    }
    let result = command.output().expect("run tyler");

    // Re-emit subprocess stdout/stderr to test process
    // This makes it controlled by cargo test's --show-output / --nocapture
    if !result.stdout.is_empty() {
        println!("{}", String::from_utf8_lossy(&result.stdout));
    }
    if !result.stderr.is_empty() {
        eprintln!("{}", String::from_utf8_lossy(&result.stderr));
    }

    result
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read json file")).expect("parse json file")
}

fn first_content_uri(tile: &Value) -> Option<&str> {
    tile["content"]["uri"].as_str().or_else(|| {
        tile["children"]
            .as_array()
            .and_then(|children| children.iter().find_map(first_content_uri))
    })
}

fn parse_tsv(text: &str) -> Vec<Vec<String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(text.as_bytes());
    reader
        .records()
        .map(|record| {
            record
                .expect("parse TSV record")
                .iter()
                .map(ToString::to_string)
                .collect()
        })
        .collect()
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

fn assert_full_tile_file_names(root: &Path, paths: &[PathBuf]) {
    for path in paths {
        let relative = path
            .strip_prefix(root)
            .expect("tile is under tile directory");
        let components = relative
            .components()
            .map(|component| component.as_os_str().to_str().expect("UTF-8 path"))
            .collect::<Vec<_>>();
        assert_eq!(
            components.len(),
            3,
            "unexpected tile path: {}",
            path.display()
        );
        let [level, x, file_name] = components.as_slice() else {
            unreachable!("checked tile path component count");
        };
        let prefix = format!("{level}-{x}-");
        let y = file_name
            .strip_prefix(&prefix)
            .and_then(|name| name.split('.').next())
            .expect("tile filename should include level-x-y");
        assert!(
            !y.is_empty() && y.chars().all(|digit| digit.is_ascii_digit()),
            "tile filename should contain full numeric tile ID: {}",
            path.display()
        );
    }
}

fn assert_tsv_content_paths(output_dir: &Path, rows: &[Vec<String>]) {
    let tile_id_ix = rows[0]
        .iter()
        .position(|column| column == "tile_id")
        .expect("metadata.tsv should contain tile_id");
    let content_path_ix = rows[0]
        .iter()
        .position(|column| column == "content_path")
        .expect("metadata.tsv should contain content_path");
    for row in rows.iter().skip(1) {
        let mut tile_id = row[tile_id_ix].split('/');
        let level = tile_id.next().expect("tile ID level");
        let x = tile_id.next().expect("tile ID x");
        let y = tile_id.next().expect("tile ID y");
        assert!(tile_id.next().is_none(), "tile_id should be level/x/y");
        assert_eq!(
            row[content_path_ix],
            format!("t/{level}/{x}/{level}-{x}-{y}/cityobjects.tsv")
        );
        assert!(output_dir.join(&row[content_path_ix]).is_file());
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
fn format_tsv_writes_tile_tables_and_aggregate_metadata() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("format-tsv", &metadata, &[feature]);
    let output_dir = unique_test_dir("format-tsv-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--format",
            "tsv",
            "--tsv-include-null-rows",
            "--tsv-include-hierarchy",
            "--tsv-include-cityjson-ordinal",
            "--tsv-include-semantics",
        ],
    );
    assert_success(&output, "TSV format run");

    let mut cityobject_tables = Vec::new();
    collect_paths_with_suffix(
        &output_dir.join("t"),
        "cityobjects.tsv",
        &mut cityobject_tables,
    );
    assert!(
        !cityobject_tables.is_empty(),
        "expected CityObject TSV tables under {}",
        output_dir.join("t").display()
    );
    let cityobjects = fs::read_to_string(&cityobject_tables[0]).expect("read cityobjects.tsv");
    let cityobjects_header = cityobjects.lines().next().expect("cityobjects header");
    assert!(cityobjects_header.contains("cityobject_id"));
    assert!(cityobjects_header.contains("cityobject_ix"));
    assert!(!cityobjects_header.contains("parents"));
    assert!(!cityobjects_header.contains("children"));

    let mut cityobject_hierarchy_tables = Vec::new();
    collect_paths_with_suffix(
        &output_dir.join("t"),
        "cityobject_hierarchy.tsv",
        &mut cityobject_hierarchy_tables,
    );
    assert!(
        !cityobject_hierarchy_tables.is_empty(),
        "expected CityObject hierarchy TSV tables under {}",
        output_dir.join("t").display()
    );
    let cityobject_hierarchy =
        fs::read_to_string(&cityobject_hierarchy_tables[0]).expect("read cityobject_hierarchy.tsv");
    let hierarchy_rows = parse_tsv(&cityobject_hierarchy);
    assert_eq!(hierarchy_rows[0].as_slice(), ["parent_id", "child_id"]);

    let mut semantic_tables = Vec::new();
    collect_paths_with_suffix(&output_dir.join("t"), "semantics.tsv", &mut semantic_tables);
    assert!(
        !semantic_tables.is_empty(),
        "expected split semantic TSV tables under {}",
        output_dir.join("t").display()
    );

    let aggregate_metadata =
        fs::read_to_string(output_dir.join("metadata.tsv")).expect("read metadata.tsv");
    let metadata_header = aggregate_metadata.lines().next().expect("metadata header");
    assert!(metadata_header.starts_with("tile_id	content_path	geographical_extent	identifier"));
    assert!(metadata_header.contains("geographical_extent"));
    assert!(aggregate_metadata
        .lines()
        .skip(1)
        .any(|line| line.contains("t/")));

    let rows = parse_tsv(&aggregate_metadata);
    assert!(
        rows.len() >= 2,
        "metadata.tsv should contain at least one data row"
    );
    let extent_wkt_ix = rows[0]
        .iter()
        .position(|column| column == "geographical_extent")
        .expect("metadata.tsv should contain geographical_extent");
    assert_tsv_content_paths(&output_dir, &rows);
    assert!(
        rows[1][extent_wkt_ix].starts_with("POLYGON(("),
        "geographical_extent should contain polygon WKT"
    );
    assert!(!output_dir.join(".tyler-tsv-metadata").exists());
    assert!(!output_dir
        .join("metadata/cjindex-metadata.city.json")
        .exists());
    assert!(!output_dir.join("tileset.json").exists());
}

#[test]
fn format_gpkg_writes_tile_databases_and_aggregate_metadata() {
    let metadata = read_fixture("resources/data/3dbag_x00.city.json");
    let feature = read_fixture("resources/data/3dbag_feature_x71.city.jsonl");
    let dataset = write_ndjson_dataset("format-gpkg", &metadata, &[feature]);
    let output_dir = unique_test_dir("format-gpkg-output");
    let legacy_metadata_dir = output_dir.join(".tyler-gpkg-metadata");
    fs::create_dir_all(&legacy_metadata_dir).expect("create legacy metadata directory");
    fs::write(legacy_metadata_dir.join("fragment.gpkg"), b"stale")
        .expect("write stale metadata fragment");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--format",
            "gpkg",
            "--tiles-dir",
            "tiles",
            "--gpkg-split-lod",
            "--gpkg-include-semantics",
            "--gpkg-include-hierarchy",
            "--gpkg-include-address",
        ],
    );
    assert_success(&output, "GeoPackage format run");

    let mut tile_databases = Vec::new();
    collect_paths_with_suffix(&output_dir.join("tiles"), ".gpkg", &mut tile_databases);
    assert_full_tile_file_names(&output_dir.join("tiles"), &tile_databases);
    assert!(
        !tile_databases.is_empty(),
        "expected GeoPackage tiles under {}",
        output_dir.join("tiles").display()
    );
    let tile = rusqlite::Connection::open(&tile_databases[0]).expect("open tile GeoPackage");
    let tables = tile
        .prepare("SELECT table_name FROM gpkg_contents ORDER BY table_name")
        .expect("prepare tile table query")
        .query_map([], |row| row.get::<_, String>(0))
        .expect("query tile tables")
        .collect::<rusqlite::Result<Vec<_>>>()
        .expect("read tile tables");
    assert!(tables.iter().any(|table| table == "semantics"));
    assert!(tables.iter().any(|table| table == "addresses"));
    assert!(tables.iter().any(|table| table == "cityobject_hierarchy"));
    assert!(tables.iter().any(|table| table == "semantic_hierarchy"));

    let aggregate_path = output_dir.join("metadata.gpkg");
    let aggregate =
        rusqlite::Connection::open(&aggregate_path).expect("open aggregate metadata GeoPackage");
    let rows = aggregate
        .prepare("SELECT tile_id, content_path, geographical_extent FROM metadata ORDER BY id")
        .expect("prepare aggregate query")
        .query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, Vec<u8>>(2)?,
            ))
        })
        .expect("query aggregate metadata")
        .collect::<rusqlite::Result<Vec<_>>>()
        .expect("read aggregate metadata");
    assert!(!rows.is_empty());
    assert_eq!(
        rows.len(),
        tile_databases.len(),
        "metadata.gpkg should contain exactly one row per tile GeoPackage"
    );
    for (tile_id, content_path, extent) in rows {
        let mut tile_id_components = tile_id.split('/');
        let level = tile_id_components.next().expect("tile ID level");
        let x = tile_id_components.next().expect("tile ID x");
        let y = tile_id_components.next().expect("tile ID y");
        assert!(
            tile_id_components.next().is_none(),
            "tile_id should be level/x/y"
        );
        assert_eq!(
            content_path,
            format!("tiles/{level}/{x}/{level}-{x}-{y}.gpkg")
        );
        assert!(output_dir.join(&content_path).is_file());
        assert!(extent.starts_with(b"GP"));
    }
    let geometry_registration: (String, i32) = aggregate
        .query_row(
            "SELECT column_name, srs_id FROM gpkg_geometry_columns WHERE table_name = 'metadata'",
            [],
            |row| Ok((row.get(0)?, row.get(1)?)),
        )
        .expect("read aggregate geometry registration");
    assert_eq!(geometry_registration.0, "geographical_extent");
    assert!(!output_dir.join(".tyler-gpkg-metadata").exists());
    assert!(!output_dir.join("tileset.json").exists());
}

/// Issue <https://github.com/3DGI/tyler/issues/137>
#[test]
fn cityjsonfeature_buildingpart_filter_does_not_duplicate_cityobjects() {
    let dataset = unique_test_dir("cityjsonfeature-buildingpart-filter");
    fs::copy(
        repo_path("tests/data/cjindex_cityjsonfeature_alias/source.city.jsonl"),
        dataset.join("source.city.jsonl"),
    )
    .expect("copy CityJSONFeature alias fixture");
    let output_dir = unique_test_dir("cityjsonfeature-buildingpart-filter-output");

    let output = run_tyler_with_rust_log(
        &dataset,
        &output_dir,
        "debug",
        &[
            "--format",
            "cityjson",
            "--color-building-part",
            "#ff0000",
            "--object-type",
            "BuildingPart",
        ],
    );
    assert_success(&output, "CityJSON BuildingPart fixture run");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("duplicate CityObject id") && !stdout.contains("duplicate CityObject id"),
        "Tyler should not duplicate CityObjects while filtering the fixture\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );

    let mut cityjson_tiles = Vec::new();
    collect_paths_with_suffix(&output_dir.join("t"), ".city.json", &mut cityjson_tiles);
    assert!(
        !cityjson_tiles.is_empty(),
        "expected CityJSON tiles under {}\nstdout:\n{}\nstderr:\n{}",
        output_dir.join("t").display(),
        stdout,
        stderr
    );
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
        &[
            "--format",
            "cityjsonseq",
            "--tiles-dir",
            "tiles",
            "--debug-dump-data",
        ],
    );
    assert_success(&output, "CityJSONSeq format run");

    let mut output_streams = Vec::new();
    collect_paths_with_suffix(
        &output_dir.join("tiles"),
        ".city.jsonl",
        &mut output_streams,
    );
    assert_full_tile_file_names(&output_dir.join("tiles"), &output_streams);
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
    assert_full_tile_file_names(&output_dir.join("debug").join("inputs"), &debug_streams);

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
    let dataset = unique_test_dir("debug-3dtiles-tileset");
    fs::copy(
        "cityjson-convert/tests/data/multi_feature_types.city.jsonl",
        dataset.join("multi_feature_types.city.jsonl"),
    )
    .unwrap();
    let output_dir = unique_test_dir("debug-tileset-only-output");

    let output = run_tyler(
        &dataset,
        &output_dir,
        &[
            "--debug-3dtiles-tileset-only",
            "--3dtiles-implicit",
            "--tiles-dir",
            "tiles",
            "--qtree-capacity",
            "1",
        ],
    );
    assert_success(&output, "tileset-only run");

    assert!(output_dir.join("tileset.json").is_file());
    let tileset = read_json(&output_dir.join("tileset.json"));
    assert_eq!(
        first_content_uri(&tileset["root"]),
        Some("tiles/{level}/{x}/{level}-{x}-{y}.glb")
    );
    assert!(
        !output_dir.join("tiles").exists(),
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
    assert_full_tile_file_names(&output_dir.join("t"), std::slice::from_ref(&glb_path));
    let tileset = read_json(&output_dir.join("tileset.json"));
    let uri = first_content_uri(&tileset["root"]).expect("explicit tileset content URI");
    let relative_path = glb_path
        .strip_prefix(&output_dir)
        .expect("GLB path under output directory")
        .to_string_lossy()
        .replace('\\', "/");
    assert_eq!(uri, relative_path);
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
