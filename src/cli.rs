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

use anyhow::Result;
use clap::Parser;

use crate::material::MaterialConfig;

/// 3D Tiles specification version for tile output.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum TilesVersion {
    /// 3D Tiles 1.1 — GLB content with EXT_mesh_features + EXT_structural_metadata (default)
    #[default]
    #[value(name = "1.1")]
    V1_1,
    /// 3D Tiles 1.0 — B3DM content with batch table
    #[value(name = "1.0")]
    V1_0,
}

impl TilesVersion {
    /// File extension for tile content files.
    pub fn extension(&self) -> &'static str {
        match self {
            TilesVersion::V1_1 => "glb",
            TilesVersion::V1_0 => "b3dm",
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
#[command(group(
    clap::ArgGroup::new("input")
        .required(true)
        .multiple(true)
        .args(["buildings", "trees"]),
))]
pub struct Cli {
    /// CityJSONL file with building features (e.g. from roofer).
    #[arg(long, value_parser = existing_canonical_path, display_order = 1)]
    pub buildings: Option<PathBuf>,
    /// GeoParquet file with tree (SolitaryVegetationObject) 3D polygons.
    /// When combined with --buildings, the trees are reprojected and quantized
    /// to the buildings' coordinate reference system.
    #[arg(long, value_parser = existing_canonical_path, display_order = 2)]
    pub trees: Option<PathBuf>,
    /// Column name in the GeoParquet file to use as the feature ID.
    /// If not specified, sequential IDs (tree-00000, tree-00001, ...) are generated.
    #[arg(long, display_order = 3)]
    pub tree_id_column: Option<String>,
    /// Directory for the output.
    #[arg(short, long, display_order = 4)]
    pub output: PathBuf,
    /// The CityObject type to use for the 3D Tiles
    /// (https://www.cityjson.org/specs/1.1.3/#the-different-city-objects).
    /// You can specify it multiple times. If not set, all types are accepted.
    #[arg(long, value_enum, display_order = 5)]
    pub object_type: Option<Vec<crate::parser::CityObjectType>>,
    /// The metadata class to assign to the property table when the output is
    /// 3D Tiles (https://github.com/CesiumGS/glTF/tree/3d-tiles-next/extensions/2.0/Vendor/EXT_structural_metadata#class).
    #[arg(long = "3dtiles-metadata-class", display_order = 10)]
    pub cesium3dtiles_metadata_class: Option<String>,
    /// Create implicit tiling when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    /// By default, explicit tiling is created for the 3D Tiles output.
    #[arg(long = "3dtiles-implicit", display_order = 11)]
    pub cesium3dtiles_implicit: bool,
    /// Generate and write the Tileset only, without exporting the glTF tiles, when the output format is 3D Tiles (https://docs.ogc.org/cs/22-025r4/22-025r4.html#toc31).
    #[arg(long = "3dtiles-tileset-only", display_order = 12)]
    pub cesium3dtiles_tileset_only: bool,
    /// Use the tile boundingVolume as the content boundingVolume, instead of calculating the content boundingVolume from the data.
    #[arg(long = "3dtiles-content-bv-from-tile", display_order = 13)]
    pub cesium3dtiles_content_bv_from_tile: bool,
    /// Add the boundingVolume of the content for the the tiles that have content.
    #[arg(long = "3dtiles-content-add-bv", display_order = 14)]
    pub cesium3dtiles_content_add_bv: bool,
    /// 3D Tiles version: "1.1" (default, GLB with EXT_mesh_features) or "1.0" (B3DM with batch table).
    #[arg(long = "3dtiles-version", value_enum, default_value_t = TilesVersion::V1_1, display_order = 9)]
    pub tiles_version: TilesVersion,
    /// Set the geometric error (see 3D Tiles specification) on the parent nodes of leafs. This controls at what
    /// camera distance leaf nodes become visible. Higher values make content visible earlier when zooming in.
    #[arg(long, short = 'e', default_value = "12", display_order = 15)]
    pub geometric_error_above_leaf: Option<f64>,
    /// Set the 2D cell size for the grid that is used for constructing the quadtree.
    /// In input units (eg. meters). Note that the cell size will be adjusted so that it is
    /// possible to construct a tightly fit square, containing 4^n cells. The final cell size will
    /// larger than this value.
    #[arg(long, default_value = "250", display_order = 16)]
    pub grid_cellsize: Option<u32>,
    /// Limit the minimum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely small z
    /// coordinates. In input units (eg. meters).
    #[arg(long, display_order = 17)]
    pub grid_minz: Option<i32>,
    /// Limit the maximum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely large z
    /// coordinates. In input units (eg. meters).
    #[arg(long, display_order = 18)]
    pub grid_maxz: Option<i32>,
    /// Export the grid into .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long, display_order = 50)]
    pub grid_export: bool,
    /// Export the grid, and also the feature centroids into .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long, display_order = 51)]
    pub grid_export_features: bool,
    /// Load instances from this directory.
    /// In debug mode, tyler writes the generated world, quadtree etc. instances to .bincode files, which later can be used for debugging.
    /// When this argument is specified, tyler will load the instances from the .bincode files that are available in the directory.
    #[arg(long, value_parser = existing_canonical_path, display_order = 52)]
    pub debug_load_data: Option<PathBuf>,
    /// The maximum number of vertices in a leaf of the quadtree.
    #[arg(long, default_value = "42000", display_order = 19)]
    pub qtree_capacity: Option<usize>,

    /// Path to a TOML file with per-CityObjectType PBR material configuration.
    /// Defines base_color, metallic_factor, and roughness_factor per type.
    /// CLI --glb-color-* flags override base_color from this file.
    #[arg(long = "material-config", value_parser = existing_canonical_path, display_order = 29)]
    pub material_config: Option<PathBuf>,

    // --- Per-CityObjectType GLB colors ---
    // Each is a hex color (#RRGGBB). The type name in the flag is case-insensitive.
    // If not specified for a type, the default pink (#FFC0CB) is used.

    /// GLB color for Building features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-building", value_parser = hex_color, display_order = 30)]
    pub glb_color_building: Option<String>,
    /// GLB color for BuildingPart features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-buildingpart", value_parser = hex_color, display_order = 30)]
    pub glb_color_building_part: Option<String>,
    /// GLB color for BuildingInstallation features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-buildinginstallation", value_parser = hex_color, display_order = 30)]
    pub glb_color_building_installation: Option<String>,
    /// GLB color for TINRelief features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-tinrelief", value_parser = hex_color, display_order = 30)]
    pub glb_color_tin_relief: Option<String>,
    /// GLB color for Road features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-road", value_parser = hex_color, display_order = 30)]
    pub glb_color_road: Option<String>,
    /// GLB color for Railway features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-railway", value_parser = hex_color, display_order = 30)]
    pub glb_color_railway: Option<String>,
    /// GLB color for TransportSquare features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-transportsquare", value_parser = hex_color, display_order = 30)]
    pub glb_color_transport_square: Option<String>,
    /// GLB color for WaterBody features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-waterbody", value_parser = hex_color, display_order = 30)]
    pub glb_color_water_body: Option<String>,
    /// GLB color for PlantCover features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-plantcover", value_parser = hex_color, display_order = 30)]
    pub glb_color_plant_cover: Option<String>,
    /// GLB color for SolitaryVegetationObject features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-solitaryvegetationobject", value_parser = hex_color, display_order = 30)]
    pub glb_color_solitary_vegetation_object: Option<String>,
    /// GLB color for LandUse features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-landuse", value_parser = hex_color, display_order = 30)]
    pub glb_color_land_use: Option<String>,
    /// GLB color for CityFurniture features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-cityfurniture", value_parser = hex_color, display_order = 30)]
    pub glb_color_city_furniture: Option<String>,
    /// GLB color for Bridge features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-bridge", value_parser = hex_color, display_order = 30)]
    pub glb_color_bridge: Option<String>,
    /// GLB color for BridgePart features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-bridgepart", value_parser = hex_color, display_order = 30)]
    pub glb_color_bridge_part: Option<String>,
    /// GLB color for BridgeInstallation features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-bridgeinstallation", value_parser = hex_color, display_order = 30)]
    pub glb_color_bridge_installation: Option<String>,
    /// GLB color for BridgeConstructionElement features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-bridgeconstructionelement", value_parser = hex_color, display_order = 30)]
    pub glb_color_bridge_construction_element: Option<String>,
    /// GLB color for Tunnel features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-tunnel", value_parser = hex_color, display_order = 30)]
    pub glb_color_tunnel: Option<String>,
    /// GLB color for TunnelPart features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-tunnelpart", value_parser = hex_color, display_order = 30)]
    pub glb_color_tunnel_part: Option<String>,
    /// GLB color for TunnelInstallation features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-tunnelinstallation", value_parser = hex_color, display_order = 30)]
    pub glb_color_tunnel_installation: Option<String>,
    /// GLB color for GenericCityObject features, hex #RRGGBB (default #FFC0CB).
    #[arg(long = "glb-color-genericcityobject", value_parser = hex_color, display_order = 30)]
    pub glb_color_generic_city_object: Option<String>,
}

use crate::parser::CityObjectType;

impl Cli {
    /// Build a MaterialConfig from the TOML file (if provided) and CLI color overrides.
    pub fn build_material_config(&self) -> Result<MaterialConfig> {
        let mut config = match &self.material_config {
            Some(path) => MaterialConfig::from_toml(path)?,
            None => MaterialConfig::default_config(),
        };

        let cli_colors: Vec<(CityObjectType, &Option<String>)> = vec![
            (CityObjectType::Building, &self.glb_color_building),
            (CityObjectType::BuildingPart, &self.glb_color_building_part),
            (CityObjectType::BuildingInstallation, &self.glb_color_building_installation),
            (CityObjectType::TINRelief, &self.glb_color_tin_relief),
            (CityObjectType::Road, &self.glb_color_road),
            (CityObjectType::Railway, &self.glb_color_railway),
            (CityObjectType::TransportSquare, &self.glb_color_transport_square),
            (CityObjectType::WaterBody, &self.glb_color_water_body),
            (CityObjectType::PlantCover, &self.glb_color_plant_cover),
            (CityObjectType::SolitaryVegetationObject, &self.glb_color_solitary_vegetation_object),
            (CityObjectType::LandUse, &self.glb_color_land_use),
            (CityObjectType::CityFurniture, &self.glb_color_city_furniture),
            (CityObjectType::Bridge, &self.glb_color_bridge),
            (CityObjectType::BridgePart, &self.glb_color_bridge_part),
            (CityObjectType::BridgeInstallation, &self.glb_color_bridge_installation),
            (CityObjectType::BridgeConstructiveElement, &self.glb_color_bridge_construction_element),
            (CityObjectType::Tunnel, &self.glb_color_tunnel),
            (CityObjectType::TunnelPart, &self.glb_color_tunnel_part),
            (CityObjectType::TunnelInstallation, &self.glb_color_tunnel_installation),
            (CityObjectType::GenericCityObject, &self.glb_color_generic_city_object),
        ];
        config.apply_cli_color_overrides(&cli_colors);

        Ok(config)
    }
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

    fn required_args() -> Vec<&'static str> {
        vec![
            "tyler-glb",
            "--buildings", env!("CARGO_MANIFEST_DIR"),
            "-o",
            env!("CARGO_MANIFEST_DIR"),
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

    #[test]
    fn verify_per_type_color() {
        let mut args = required_args();
        args.extend(&["--glb-color-building", "#FF0000", "--glb-color-solitaryvegetationobject", "#00FF00"]);
        let cli = Cli::try_parse_from(args).unwrap();
        let config = cli.build_material_config().unwrap();
        let building_color = config.color_map[&crate::parser::CityObjectType::Building];
        assert!((building_color[0] - 1.0).abs() < 0.01); // red = 1.0
        assert!(building_color[1] < 0.01); // green = 0.0
        let tree_color = config.color_map[&crate::parser::CityObjectType::SolitaryVegetationObject];
        assert!(tree_color[0] < 0.01); // red = 0.0
        assert!((tree_color[1] - 1.0).abs() < 0.01); // green = 1.0
    }

    #[test]
    fn verify_tiles_version_default() {
        let args = required_args();
        let cli = Cli::try_parse_from(args).unwrap();
        assert_eq!(cli.tiles_version, crate::cli::TilesVersion::V1_1);
    }

    #[test]
    fn verify_tiles_version_1_0() {
        let mut args = required_args();
        args.extend(&["--3dtiles-version", "1.0"]);
        let cli = Cli::try_parse_from(args).unwrap();
        assert_eq!(cli.tiles_version, crate::cli::TilesVersion::V1_0);
    }
}
