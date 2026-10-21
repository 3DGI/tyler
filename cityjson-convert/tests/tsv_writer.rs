use std::fs;
use std::path::PathBuf;

use cityjson_convert::{
    convert_to_tsv, tabulate_addresses, tabulate_cityobject_hierarchy, tabulate_cityobjects,
    tabulate_model_metadata, tabulate_semantic_hierarchy, tabulate_semantic_primitives,
    write_addresses_tsv, write_cityobject_hierarchy_tsv, write_cityobjects_tsv, write_metadata_tsv,
    write_semantic_hierarchy_tsv, write_semantics_tsv, TsvExportOptions, TsvWriteOptions,
};
use cityjson_lib::cityjson_types::v2_0::OwnedAttributeValue;
use cityjson_lib::json;

fn parse_tsv(bytes: &[u8]) -> Vec<Vec<String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes);
    let headers = reader
        .headers()
        .expect("read TSV headers")
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let mut rows = vec![headers];
    rows.extend(reader.records().map(|record| {
        record
            .expect("read TSV row")
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
    }));
    rows
}

fn temp_output_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "cityjson_convert_{name}_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn writes_cityobjects_tsv_with_case_insensitive_collisions() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building":{
                    "type":"Building",
                    "attributes":{
                        "eindRegistratie":"first",
                        "eindregistratie":"second"
                    }
                }
            },
            "vertices":[]
        }"#,
    )
    .expect("parse collision CityJSON");
    let table = tabulate_cityobjects(&model).unwrap();

    let mut bytes = Vec::new();
    write_cityobjects_tsv(&table, &TsvWriteOptions::default(), &mut bytes).unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(
        rows[0],
        [
            "cityobject_id",
            "cityobject_type",
            "attributes__eindRegistratie",
            "attributes__eindregistratie__2"
        ]
    );
    assert_eq!(rows[1][2], "first");
    assert_eq!(rows[1][3], "second");
}

#[test]
fn omits_geometry_ref_attributes_from_cityobjects_tsv() {
    let mut model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building":{
                    "type":"Building",
                    "geometry":[{
                        "type":"MultiSurface",
                        "lod":"1",
                        "boundaries":[[[0,1,2,0]]]
                    }]
                }
            },
            "vertices":[[0,0,0],[1,0,0],[0,1,0]]
        }"#,
    )
    .expect("parse inline CityJSON");
    let geometry_handle = model
        .cityobjects()
        .iter()
        .next()
        .and_then(|(_, object)| object.geometry())
        .and_then(|geometries| geometries.first().copied())
        .expect("geometry handle");
    let (_, cityobject) = model
        .cityobjects_mut()
        .iter_mut()
        .next()
        .expect("cityobject");
    cityobject.attributes_mut().insert(
        "location".to_string(),
        OwnedAttributeValue::Geometry(geometry_handle),
    );

    let table = tabulate_cityobjects(&model).expect("tabulate CityObjects");
    let mut bytes = Vec::new();
    write_cityobjects_tsv(&table, &TsvWriteOptions::default(), &mut bytes).unwrap();
    let rows = parse_tsv(&bytes);

    assert!(!rows[0].contains(&"attributes__location".to_string()));
}

#[test]
fn writes_include_address_tsv_with_dynamic_columns() {
    let mut model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building":{
                    "type":"Building",
                    "geometry":[{
                        "type":"MultiPoint",
                        "lod":"1",
                        "boundaries":[0]
                    }]
                }
            },
            "vertices":[[4,5,6]]
        }"#,
    )
    .expect("parse inline CityJSON");
    let geometry_handle = model
        .cityobjects()
        .iter()
        .next()
        .and_then(|(_, object)| object.geometry())
        .and_then(|geometries| geometries.first().copied())
        .expect("geometry handle");
    let (_, cityobject) = model
        .cityobjects_mut()
        .iter_mut()
        .next()
        .expect("cityobject");
    cityobject.extra_mut().insert(
        "address".to_string(),
        OwnedAttributeValue::Vec(vec![OwnedAttributeValue::Map(
            std::collections::HashMap::from([
                (
                    "location".to_string(),
                    OwnedAttributeValue::Geometry(geometry_handle),
                ),
                (
                    "street".to_string(),
                    OwnedAttributeValue::String("Main Street".to_string()),
                ),
                (
                    "houseNumber".to_string(),
                    OwnedAttributeValue::String("7".to_string()),
                ),
            ]),
        )]),
    );

    let cityobjects = tabulate_cityobjects(&model).expect("tabulate CityObjects");
    let mut cityobjects_bytes = Vec::new();
    write_cityobjects_tsv(
        &cityobjects,
        &TsvWriteOptions {
            include_null_rows: true,
            include_hierarchy: false,
            include_cityjson_ordinal: false,
        },
        &mut cityobjects_bytes,
    )
    .unwrap();
    let cityobjects_rows = parse_tsv(&cityobjects_bytes);
    assert_eq!(cityobjects_rows[0], ["cityobject_id", "cityobject_type"]);

    let table = tabulate_addresses(&model).expect("tabulate addresses");
    let mut bytes = Vec::new();
    write_addresses_tsv(
        &table,
        &TsvWriteOptions {
            include_null_rows: false,
            include_hierarchy: false,
            include_cityjson_ordinal: true,
        },
        &mut bytes,
    )
    .unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[0],
        ["cityobject_id", "cityobject_type", "houseNumber", "street"]
    );
    assert_eq!(rows[1], ["building", "Building", "7", "Main Street"]);
    assert!(!rows[0].contains(&"location".to_string()));
    assert!(!rows[0].contains(&"geom".to_string()));
}

#[test]
fn writes_cityobjects_tsv_with_filtering_ordinal_hierarchy_and_json_cells() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{
                "building":{
                    "type":"Building",
                    "children":["room"],
                    "attributes":{"name":"Library","scores":[1,2]}
                },
                "room":{
                    "type":"BuildingRoom",
                    "parents":["building"],
                    "attributes":{"name":null}
                }
            },
            "vertices":[]
        }"#,
    )
    .unwrap();
    let table = tabulate_cityobjects(&model).unwrap();

    let mut bytes = Vec::new();
    write_cityobjects_tsv(
        &table,
        &TsvWriteOptions {
            include_null_rows: false,
            include_hierarchy: true,
            include_cityjson_ordinal: true,
        },
        &mut bytes,
    )
    .unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[0],
        [
            "cityobject_id",
            "cityobject_type",
            "cityobject_ix",
            "attributes__name",
            "attributes__scores",
        ]
    );
    assert_eq!(rows[1][0], "building");
    assert_eq!(rows[1][2], "0");
    assert_eq!(rows[1][3], "Library");
    assert_eq!(rows[1][4], "[1,2]");

    let hierarchy = tabulate_cityobject_hierarchy(&model).unwrap();
    let mut bytes = Vec::new();
    write_cityobject_hierarchy_tsv(&hierarchy, &mut bytes).unwrap();
    let rows = parse_tsv(&bytes);
    assert_eq!(rows, [["parent_id", "child_id"], ["building", "room"]]);

    let mut bytes = Vec::new();
    write_cityobjects_tsv(
        &table,
        &TsvWriteOptions {
            include_null_rows: true,
            include_hierarchy: false,
            include_cityjson_ordinal: false,
        },
        &mut bytes,
    )
    .unwrap();
    assert_eq!(parse_tsv(&bytes).len(), 3);
}

#[test]
fn writes_metadata_tsv_with_fixed_fields_extent_and_extra() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{},
            "vertices":[],
            "metadata":{
                "identifier":"dataset-1",
                "referenceDate":"2026-01-02",
                "referenceSystem":"EPSG:7415",
                "title":"Demo",
                "geographicalExtent":[1.0,2.0,3.0,4.0,5.0,6.0],
                "+quality":{"score":7}
            }
        }"#,
    )
    .unwrap();
    let table = tabulate_model_metadata(&model).unwrap();

    let mut bytes = Vec::new();
    write_metadata_tsv(&table, &mut bytes).unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(rows.len(), 2);
    assert!(rows[0].contains(&"geographical_extent_wkt".to_string()));
    assert!(!rows[0].contains(&"geographical_extent".to_string()));
    assert!(!rows[0].contains(&"geographical_extent_wkb".to_string()));
    assert!(rows[0].contains(&"+quality__score".to_string()));
    assert!(!rows[0].contains(&"metadata_extra__+quality__score".to_string()));
    assert_eq!(rows[1][0], "dataset-1");
    assert_eq!(rows[1][4], "POLYGON((1 2, 4 2, 4 5, 1 5, 1 2))");
    assert_eq!(rows[1].last().unwrap(), "7");
}

#[test]
fn writes_empty_metadata_extent_cell_when_extent_is_absent() {
    let model = json::from_slice(
        br#"{
            "type":"CityJSON",
            "version":"2.0",
            "CityObjects":{},
            "vertices":[],
            "metadata":{"identifier":"dataset-1"}
        }"#,
    )
    .unwrap();
    let table = tabulate_model_metadata(&model).unwrap();

    let mut bytes = Vec::new();
    write_metadata_tsv(&table, &mut bytes).unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(rows[1][4], "");
}

#[test]
fn writes_semantics_tsv_with_primitive_rows() {
    let model = json::from_slice(semantic_fixture()).unwrap();
    let semantics = tabulate_semantic_primitives(&model).unwrap();

    let mut bytes = Vec::new();
    write_semantics_tsv(
        &semantics,
        &TsvWriteOptions {
            include_null_rows: true,
            include_hierarchy: true,
            include_cityjson_ordinal: true,
        },
        &mut bytes,
    )
    .unwrap();
    let rows = parse_tsv(&bytes);
    assert_eq!(rows.len(), 3);
    assert_eq!(
        rows[0],
        [
            "cityobject_id",
            "geometry_id",
            "semantic_id",
            "primitive_ix",
            "geometry_type",
            "geometry_lod",
            "semantic_type",
            "attribute__slope",
        ]
    );
    assert_eq!(
        rows[1],
        [
            "building",
            "0",
            "0",
            "0",
            "MultiSurface",
            "2.0",
            "RoofSurface",
            "30"
        ]
    );
    assert_eq!(rows[2][2], "1");
    assert_eq!(rows[2][6], "WallSurface");

    let hierarchy = tabulate_semantic_hierarchy(&model);
    let mut hierarchy_bytes = Vec::new();
    write_semantic_hierarchy_tsv(&hierarchy, &mut hierarchy_bytes).unwrap();
    let rows = parse_tsv(&hierarchy_bytes);
    assert_eq!(rows, [["parent_id", "child_id"], ["0", "1"]]);
}

#[test]
fn writes_semantics_tsv_filtered_to_assigned_primitives() {
    let model = json::from_slice(semantic_fixture()).unwrap();
    let semantics = tabulate_semantic_primitives(&model).unwrap();

    let mut bytes = Vec::new();
    write_semantics_tsv(
        &semantics,
        &TsvWriteOptions {
            include_null_rows: false,
            include_hierarchy: true,
            include_cityjson_ordinal: true,
        },
        &mut bytes,
    )
    .unwrap();
    let rows = parse_tsv(&bytes);

    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0][0], "cityobject_id");
    assert_eq!(rows[0][2], "semantic_id");
    assert_eq!(rows[0][3], "primitive_ix");
    assert!(rows[0].contains(&"semantic_type".to_string()));
    assert!(rows[0].contains(&"attribute__slope".to_string()));
    assert_eq!(rows[1][0], "building");
    assert_eq!(rows[1][2], "0");
    assert_eq!(rows[1][6], "RoofSurface");
    assert_eq!(rows[1].last().unwrap(), "30");
    assert_eq!(rows[2][2], "1");
}

#[test]
fn converts_model_to_tsv_file_outputs() {
    let model = json::from_slice(semantic_fixture()).unwrap();
    let dir = temp_output_dir("tsv_convert");
    let output = dir.join("nested").join("model.tsv");

    let mut previous = Vec::new();
    for suffix in [None, Some("lod2.2"), Some("LoD_3-v1")] {
        convert_to_tsv(
            &model,
            &output,
            &TsvExportOptions {
                output_suffix: suffix.map(str::to_owned),
                include_null_rows: true,
                include_hierarchy: true,
                include_cityjson_ordinal: true,
                include_metadata: true,
                include_semantics: true,
                include_address: true,
            },
        )
        .unwrap();
        for (path, bytes) in &previous {
            assert_eq!(&fs::read(path).unwrap(), bytes);
        }
        for table in [
            "",
            "_metadata",
            "_semantics",
            "_addresses",
            "_cityobject_hierarchy",
            "_semantic_hierarchy",
        ] {
            let label = suffix
                .map(|suffix| format!("-{suffix}"))
                .unwrap_or_default();
            let path = dir.join("nested").join(format!("model{table}{label}.tsv"));
            let bytes = fs::read(&path).unwrap();
            assert!(!bytes.is_empty());
            previous.push((path, bytes));
        }
    }

    fs::remove_dir_all(dir).unwrap();
}

fn semantic_fixture() -> &'static [u8] {
    br#"{
        "type":"CityJSON",
        "version":"2.0",
        "CityObjects":{
            "building":{
                "type":"Building",
                "attributes":{"name":"Library"},
                "geometry":[{
                    "type":"MultiSurface",
                    "lod":"2.0",
                    "boundaries":[[[0,1,2]],[[0,2,3]],[[1,2,3]]],
                    "semantics":{
                        "surfaces":[
                            {"type":"RoofSurface","children":[1],"slope":30},
                            {"type":"WallSurface","parent":0}
                        ],
                        "values":[0,1,null]
                    }
                }]
            }
        },
        "vertices":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],
        "metadata":{"identifier":"semantic-fixture"}
    }"#
}
