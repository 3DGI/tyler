// Copyright 2023 Balázs Dukai, Ravi Peters
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
use std::path::{Path, PathBuf};

use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Cli {
    /// Input dataset root. Must be a cjindex-compatible dataset in NDJSON,
    /// CityJSON, or feature-files layout.
    #[arg(value_parser = existing_canonical_path)]
    pub input: PathBuf,
    /// Directory for the output.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Output format.
    #[arg(long, value_enum, default_value_t)]
    pub format: crate::OutputFormatKind,
    /// The CityObject type to include
    /// (https://www.cityjson.org/specs/2.0.1/#city-objects).
    /// You can specify it multiple times.
    #[arg(long, value_enum)]
    pub object_type: Option<Vec<crate::parser::CityObjectType>>,
    /// The CityObject attribute name and value type to include as feature attribute when the
    /// output is 3D Tiles. Format: <attribute_name>:<attribute_type>,<attribute_name>:<attribute_type>
    /// eg: 'name1:string,name2:bool'.
    /// Possible value types are, 'bool', 'int', 'float', 'string'.
    #[arg(long, value_delimiter = ',')]
    pub object_attributes: Option<Vec<String>>,
    /// The metadata class to assign to the property table when the output is
    /// 3D Tiles (https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata#class).
    #[arg(long = "3dtiles-metadata-class", default_value = "citymodel")]
    pub cesium3dtiles_metadata_class: String,
    /// Create implicit tiling when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    /// By default, explicit tiling is created for the 3D Tiles output.
    #[arg(long = "3dtiles-implicit")]
    pub cesium3dtiles_implicit: bool,
    /// Generate and write the Tileset only, without exporting the glTF tiles, when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    #[arg(long = "debug-3dtiles-tileset-only")]
    pub debug_cesium3dtiles_tileset_only: bool,
    /// Use the tile boundingVolume as the content boundingVolume, instead of calculating the content boundingVolume from the data.
    #[arg(long = "3dtiles-content-bv-from-tile")]
    pub cesium3dtiles_content_bv_from_tile: bool,
    /// Clip tile content geometry to the tile bounding box before writing the GLB.
    #[arg(long = "3dtiles-content-clip-to-tile-bounds")]
    pub cesium3dtiles_content_clip_to_tile_bounds: bool,
    /// Add the boundingVolume of the content for the the tiles that have content.
    #[arg(long = "3dtiles-content-add-bv")]
    pub cesium3dtiles_content_add_bv: bool,
    /// Set the geometric error factor (see 3D Tiles specification) for internal tiles.
    /// Internal tile geometricError is computed as tile width multiplied by this factor.
    /// Higher values make detailed content visible earlier when zooming in.
    #[arg(long = "3dtiles-geometric-error-factor", default_value = "0.024")]
    pub cesium3dtiles_geometric_error_factor: f64,
    /// Set the 2D cell size for the grid that is used for constructing the quadtree.
    /// In input units (eg. meters). Note that the cell size will be adjusted so that it is
    /// possible to construct a tightly fit square, containing 4^n cells. The final cell size will
    /// larger than this value.
    #[arg(long, default_value = "250")]
    pub grid_cellsize: Option<u32>,
    /// Load grid.tsv debug output and use it for quadtree construction. If a sibling
    /// features.tsv exists, it is also loaded.
    #[arg(long = "debug-load-grid", value_parser = existing_canonical_path)]
    pub debug_load_grid: Option<PathBuf>,
    /// Limit the minimum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely small z
    /// coordinates. In input units (eg. meters).
    #[arg(long)]
    pub grid_minz: Option<i32>,
    /// Limit the maximum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely large z
    /// coordinates. In input units (eg. meters).
    #[arg(long)]
    pub grid_maxz: Option<i32>,
    /// Export the grid into .tsv files in the debug output directory.
    #[arg(long = "debug-dump-grid")]
    pub debug_dump_grid: bool,
    /// Export the grid, and also the feature centroids into .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long = "debug-dump-grid-features")]
    pub debug_dump_grid_features: bool,
    /// Load instances from this directory.
    /// In debug mode, tyler writes the generated world, quadtree etc. instances to .bincode files, which later can be used for debugging.
    /// When this argument is specified, tyler will load the instances from the .bincode files that are available in the directory.
    #[arg(long, value_parser = existing_canonical_path)]
    pub debug_load_data: Option<PathBuf>,
    /// The maximum number of vertices in a leaf of the quadtree.
    #[arg(long, default_value = "42000")]
    pub qtree_capacity: Option<usize>,
    /// Dump debug data for the current run, including `.bincode` snapshots and
    /// per-tile intermediary CityJSONFeature streams.
    #[arg(long = "debug-dump-data")]
    pub debug_dump_data: bool,
    /// Copy inherited parent attributes onto geometry-bearing CityObjects before GLB conversion.
    /// This is useful for 3DBAG/roofer-style parent-child models where attributes live on the parent.
    #[arg(long)]
    pub include_parent_attributes: bool,
    /// Compute smooth vertex normals.
    #[arg(long)]
    pub smooth_normals: bool,
    /// LoD to use in output for Building features
    #[arg(long)]
    pub lod_building: Option<String>,
    /// LoD to use in output for building_part features
    #[arg(long)]
    pub lod_building_part: Option<String>,
    /// LoD to use in output for building_installation features
    #[arg(long)]
    pub lod_building_installation: Option<String>,
    /// LoD to use in output for tin_relief features
    #[arg(long)]
    pub lod_tin_relief: Option<String>,
    /// LoD to use in output for road features
    #[arg(long)]
    pub lod_road: Option<String>,
    /// LoD to use in output for railway features
    #[arg(long)]
    pub lod_railway: Option<String>,
    /// LoD to use in output for transport_square features
    #[arg(long)]
    pub lod_transport_square: Option<String>,
    /// LoD to use in output for water_body features
    #[arg(long)]
    pub lod_water_body: Option<String>,
    /// LoD to use in output for plant_cover features
    #[arg(long)]
    pub lod_plant_cover: Option<String>,
    /// LoD to use in output for solitary_vegetation_object features
    #[arg(long)]
    pub lod_solitary_vegetation_object: Option<String>,
    /// LoD to use in output for land_use features
    #[arg(long)]
    pub lod_land_use: Option<String>,
    /// LoD to use in output for city_furniture features
    #[arg(long)]
    pub lod_city_furniture: Option<String>,
    /// LoD to use in output for bridge features
    #[arg(long)]
    pub lod_bridge: Option<String>,
    /// LoD to use in output for bridge_part features
    #[arg(long)]
    pub lod_bridge_part: Option<String>,
    /// LoD to use in output for bridge_installation features
    #[arg(long)]
    pub lod_bridge_installation: Option<String>,
    /// LoD to use in output for bridge_construction_element features
    #[arg(long)]
    pub lod_bridge_construction_element: Option<String>,
    /// LoD to use in output for tunnel features
    #[arg(long)]
    pub lod_tunnel: Option<String>,
    /// LoD to use in output for tunnel_part features
    #[arg(long)]
    pub lod_tunnel_part: Option<String>,
    /// LoD to use in output for tunnel_installation features
    #[arg(long)]
    pub lod_tunnel_installation: Option<String>,
    /// LoD to use in output for lod_generic_city_object features
    #[arg(long)]
    pub lod_generic_city_object: Option<String>,
    /// Color for Building features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_building: Option<String>,
    /// Color for BuildingPart features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_building_part: Option<String>,
    /// Color for BuildingInstallation features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_building_installation: Option<String>,
    /// Color for TINRelief features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_tin_relief: Option<String>,
    /// Color for Road features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_road: Option<String>,
    /// Color for Railway features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_railway: Option<String>,
    /// Color for TransportSquare features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_transport_square: Option<String>,
    /// Color for WaterBody features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_water_body: Option<String>,
    /// Color for PlantCover features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_plant_cover: Option<String>,
    /// Color for SolitaryVegetationObject features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_solitary_vegetation_object: Option<String>,
    /// Color for LandUse features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_land_use: Option<String>,
    /// Color for CityFurniture features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_city_furniture: Option<String>,
    /// Color for Bridge features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_bridge: Option<String>,
    /// Color for BridgePart features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_bridge_part: Option<String>,
    /// Color for BridgeInstallation features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_bridge_installation: Option<String>,
    /// Color for BridgeConstructionElement features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_bridge_construction_element: Option<String>,
    /// Color for Tunnel features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_tunnel: Option<String>,
    /// Color for TunnelPart features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_tunnel_part: Option<String>,
    /// Color for TunnelInstallation features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_tunnel_installation: Option<String>,
    /// Color for GenericCityObject features specified as a hex rgb-color value, eg. #FF0000 is red.
    #[arg(long, value_parser = hex_color)]
    pub color_generic_city_object: Option<String>,
    // The number of levels to export as content from the quadtree.
    // Counted from the leaves.
    // #[arg(long, default_value = "0")]
    // pub qtree_export_levels: Option<u16>,
    // /// The criteria to check for the quadtree leaf capacity.
    // #[arg(long, value_enum, default_value = "vertices")]
    // pub qtree_criteria: Option<crate::spatial_structs::QuadTreeCriteria>,
    // /// Path to the python interpreter (>=3.8) to use for generating CityJSON tiles.
    // /// The interpreter must have a recent cjio (https://github.com/cityjson/cjio)
    // /// installed.
    // #[arg(long, value_parser = existing_path)]
    // pub exe_python: Option<PathBuf>,
}

fn existing_canonical_path(s: &str) -> Result<PathBuf, String> {
    if let Ok(c) = Path::new(s).canonicalize() {
        if c.exists() {
            Ok(c)
        } else {
            Err(format!("path {:?} does not exist", &c))
        }
    } else {
        Err(format!("could not resolve the path {:?}", s))
    }
}

/// Checks is `s` constains a 6 digit hexadecimal value preceded by a '#', eg. #FF0000
fn hex_color(s: &str) -> Result<String, String> {
    if s.len() != 7 || !s.starts_with('#') {
        return Err(String::from(
            "Input must be a 6-digit hexadecimal value preceded by a '#'",
        ));
    }
    let hex_digits = &s[1..];
    if !hex_digits.chars().all(|c| c.is_ascii_hexdigit()) {
        return Err(String::from(
            "Input must be a 6-digit hexadecimal value preceded by a '#'",
        ));
    }
    Ok(String::from(s))
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn resource_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("data")
            .join(name)
    }

    fn unique_test_dir(prefix: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("tyler-{prefix}-{unique}"));
        fs::create_dir_all(&path).expect("create test dir");
        path
    }

    fn legacy_dataset_dir() -> PathBuf {
        let dataset_dir = unique_test_dir("cli-legacy");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        fs::copy(
            resource_path("3dbag_x00.city.json"),
            dataset_dir.join("metadata.city.json"),
        )
        .expect("copy metadata");
        fs::copy(
            resource_path("3dbag_feature_x71.city.jsonl"),
            features_dir.join("sample.city.jsonl"),
        )
        .expect("copy feature");
        dataset_dir
    }

    fn required_args(input_dir: &Path) -> Vec<String> {
        vec![
            "tyler".to_string(),
            input_dir.display().to_string(),
            "-o".to_string(),
            env!("CARGO_MANIFEST_DIR").to_string(),
        ]
    }

    fn dataset_args() -> Vec<String> {
        vec![
            "tyler".to_string(),
            concat!(env!("CARGO_MANIFEST_DIR"), "/resources/data").to_string(),
            "-o".to_string(),
            env!("CARGO_MANIFEST_DIR").to_string(),
        ]
    }

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert()
    }

    /// CLI arguments review item 7: the CityJSON v2.0 `--object-type` parser
    /// still accepts repeated occurrences and stores all requested types.
    #[test]
    fn verify_object_types() {
        let dataset_dir = legacy_dataset_dir();
        let types: Vec<String> = vec![
            "--object-type".to_string(),
            "Building".to_string(),
            "--object-type".to_string(),
            "PlantCover".to_string(),
        ];
        let mut args = required_args(&dataset_dir);
        args.extend(types);
        let cli = Cli::try_parse_from(args).unwrap();
        let otypes = &cli.object_type.unwrap();
        assert!(otypes.contains(&crate::parser::CityObjectType::Building));
        assert!(otypes.contains(&crate::parser::CityObjectType::PlantCover));
    }

    /// CLI arguments review item 10: `--include-parent-attributes` remains the
    /// supported flag after removing `--bag3d-attributes-per-part`.
    #[test]
    fn verify_include_parent_attributes_flag() {
        let dataset_dir = legacy_dataset_dir();
        let mut args = required_args(&dataset_dir);
        args.push("--include-parent-attributes".to_string());
        let cli = Cli::try_parse_from(args).unwrap();
        assert!(cli.include_parent_attributes);
    }

    /// CLI arguments review item 12: `--3dtiles-geometric-error-factor`
    /// replaces the old flag name and still parses the configured float value.
    #[test]
    fn verify_geometric_error_factor() {
        let dataset_dir = legacy_dataset_dir();
        let mut args = required_args(&dataset_dir);
        args.extend([
            "--3dtiles-geometric-error-factor".to_string(),
            "0.05".to_string(),
        ]);
        let cli = Cli::try_parse_from(args).unwrap();
        assert!((cli.cesium3dtiles_geometric_error_factor - 0.05).abs() < f64::EPSILON);
    }

    /// CLI arguments review item 8: `--object-attributes` now accepts a
    /// comma-separated `name:type` mapping list and splits it into entries.
    #[test]
    fn verify_object_attributes_csv_mapping() {
        let dataset_dir = legacy_dataset_dir();
        let mut args = required_args(&dataset_dir);
        args.extend([
            "--object-attributes".to_string(),
            "name:string,levels:int,occupied:bool".to_string(),
        ]);
        let cli = Cli::try_parse_from(args).unwrap();
        assert_eq!(
            cli.object_attributes,
            Some(vec![
                "name:string".to_string(),
                "levels:int".to_string(),
                "occupied:bool".to_string()
            ])
        );
    }

    #[test]
    fn verify_single_input_arg_for_dataset_root() {
        let cli = Cli::try_parse_from(dataset_args()).unwrap();
        assert_eq!(
            cli.input,
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("resources")
                .join("data")
                .canonicalize()
                .unwrap()
        );
    }

    #[test]
    fn reject_legacy_input_flags() {
        let dataset_dir = legacy_dataset_dir();
        let args = vec![
            "tyler".to_string(),
            "--metadata".to_string(),
            dataset_dir.join("metadata.city.json").display().to_string(),
            "--features".to_string(),
            dataset_dir.join("features").display().to_string(),
            "-o".to_string(),
            env!("CARGO_MANIFEST_DIR").to_string(),
        ];
        assert!(Cli::try_parse_from(args).is_err());
    }
}
