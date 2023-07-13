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

#[derive(Parser)]
#[command(author, version, about)]
pub struct Cli {
    /// Main CityJSON file (.city.json), containing the coordinate reference system and
    /// transformation properties.
    #[arg(short, long, value_parser = existing_canonical_path)]
    pub metadata: PathBuf,
    /// Directory of CityJSONFeatures (.city.jsonl). The directory and all its
    /// subdirectories are searched recursively for feature files.
    #[arg(short, long, value_parser = existing_canonical_path)]
    pub features: PathBuf,
    /// Directory for the output.
    #[arg(short, long)]
    pub output: PathBuf,
    // /// Output format.
    // #[arg(long, value_enum)]
    // pub format: crate::Formats,
    /// The CityObject type to use for the 3D Tiles
    /// (https://www.cityjson.org/specs/1.1.3/#the-different-city-objects).
    /// You can specify it multiple times.
    #[arg(long, value_enum)]
    pub object_type: Option<Vec<crate::parser::CityObjectType>>,
    /// The CityObject attribute name and value type to include as feature attribute when the
    /// output is 3D Tiles. Format: <attribute_name>:<attribute_type> eg: 'name1:string'.
    /// Possible value types are, 'bool', 'int', 'float', 'string'.
    /// You can specify it multiple times.
    #[arg(long)]
    pub object_attribute: Option<Vec<String>>,
    /// The CityObject attribute
    /// The metadata class to assign to the property table when the output is
    /// 3D Tiles (https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata#class).
    #[arg(long = "3dtiles-metadata-class")]
    pub cesium3dtiles_metadata_class: Option<String>,
    /// Create implicit tiling when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    /// By default, explicit tiling is created for the 3D Tiles output.
    #[arg(long = "3dtiles-implicit")]
    pub cesium3dtiles_implicit: bool,
    /// Generate and write the Tileset only, without exporting the glTF tiles, when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    #[arg(long = "3dtiles-tileset-only")]
    pub cesium3dtiles_tileset_only: bool,
    /// Set the geometric error (see 3D Tiles specification) on the parent nodes of leafs. This controls at what
    /// camera distance leaf nodes become visible. Higher values make content visible earlier when zooming in.
    #[arg(long, short = 'e', default_value = "12")]
    pub geometric_error_above_leaf: Option<f64>,
    /// Set the 2D cell size for the grid that is used for constructing the quadtree. In input units (eg. meters).
    #[arg(long, default_value = "250")]
    pub grid_cellsize: Option<u16>,
    /// Generate the quadtree directly from a grid.tsv file, skipping the extent computation and feature indexing. A grid.tsv file is created with the --grid-export option. Used for debugging.
    #[arg(long)]
    pub grid_file: Option<String>,
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
    /// Export the grid into .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long)]
    pub grid_export: bool,
    /// Export the grid, and also the feature centroids into .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long)]
    pub grid_export_features: bool,
    /// The maximum number of vertices in a leaf of the quadtree.
    #[arg(long, default_value = "42000")]
    pub qtree_capacity: Option<usize>,
    /// Path to the geoflow executable for clipping and exporting the gltf files.
    #[arg(long, value_parser = existing_path)]
    pub exe_geof: Option<PathBuf>,
    /// Use mesh simplification to reduce the number of vertices per object by this fraction. Value should be a float between 0.0 (100% reduction) and 1.0 (do not use simplification). Ignored for building object types.
    #[arg(long, default_value = "0.05")]
    pub reduce_vertices: Option<f64>,
    /// Wait for the tile conversion process to finish, or terminate it if it is not finished after the provided number of seconds.
    #[arg(long)]
    pub timeout: Option<u64>,
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

/// We don't want to canonicalize paths to executables, especially a python exe from a
/// virtualenv, because the symlink would get resolved and we would end up with a path
/// to the python interpreter that was used for creating the virtualenv, and not the
/// interpreter that links to the virtualenv.
fn existing_path(s: &str) -> Result<PathBuf, String> {
    let p = Path::new(s).to_path_buf();
    if p.exists() {
        Ok(p)
    } else {
        Err(format!("path {:?} does not exist", &p))
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

    fn required_args() -> Vec<&'static str> {
        vec![
            "tyler",
            "-m",
            "metadata.city.json",
            "-f",
            env!("CARGO_MANIFEST_DIR"),
            "-o",
            env!("CARGO_MANIFEST_DIR"),
            "--format",
            "3dtiles",
        ]
    }

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert()
    }

    /// Can we pass multiple CityObject types?
    #[test]
    fn verify_object_types() {
        let mut types: Vec<&'static str> =
            vec!["--object-type", "Building", "--object-type", "PlantCover"];
        let mut args = required_args();
        args.append(&mut types);
        let cli = Cli::try_parse_from(args).unwrap();
        let otypes = &cli.object_type.unwrap();
        assert!(otypes.contains(&crate::parser::CityObjectType::Building));
        assert!(otypes.contains(&crate::parser::CityObjectType::PlantCover));
    }
}
