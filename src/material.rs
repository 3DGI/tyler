use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Context, Result};
use log::warn;
use serde::Deserialize;

use crate::parser::CityObjectType;

/// Default GLB color (pink) used when no per-type color is specified.
pub const DEFAULT_BASE_COLOR: &str = "#FFC0CB";
/// Default metallic factor matching current behavior (0.0 = dielectric).
pub const DEFAULT_METALLIC: f32 = 0.0;
/// Default roughness factor matching current behavior (128/255 ≈ 0.502).
pub const DEFAULT_ROUGHNESS: f32 = 128.0 / 255.0;

/// Convert a validated hex color string (#RRGGBB) to [f32; 4] RGBA with alpha = 1.0.
pub fn hex_to_rgba(hex: &str) -> Result<[f32; 4]> {
    if hex.len() != 7 || !hex.starts_with('#') {
        bail!("Invalid hex color format '{}': expected #RRGGBB", hex);
    }
    let hex_digits = &hex[1..];
    if !hex_digits.chars().all(|c| c.is_ascii_hexdigit()) {
        bail!("Invalid hex color '{}': contains non-hex characters", hex);
    }
    let r = u8::from_str_radix(&hex_digits[0..2], 16)? as f32 / 255.0;
    let g = u8::from_str_radix(&hex_digits[2..4], 16)? as f32 / 255.0;
    let b = u8::from_str_radix(&hex_digits[4..6], 16)? as f32 / 255.0;
    Ok([r, g, b, 1.0])
}

/// PBR material parameters for a single CityObjectType (TOML deserialization target).
#[derive(Debug, Clone, Deserialize)]
pub struct PbrMaterialParams {
    /// Base color as hex string "#RRGGBB".
    pub base_color: Option<String>,
    /// Metallic factor, 0.0 (dielectric) to 1.0 (metal).
    pub metallic_factor: Option<f32>,
    /// Roughness factor, 0.0 (smooth/glossy) to 1.0 (rough/matte).
    pub roughness_factor: Option<f32>,
}

/// Top-level TOML config file structure.
/// Section names must match CityObjectType variant names (PascalCase).
#[derive(Debug, Clone, Deserialize)]
pub struct MaterialConfigFile {
    /// Fallback PBR parameters for types not explicitly listed.
    pub defaults: Option<PbrMaterialParams>,
    /// Per-CityObjectType material definitions (captured via serde flatten).
    #[serde(flatten)]
    pub types: HashMap<String, PbrMaterialParams>,
}

/// Resolved material configuration used at runtime.
#[derive(Debug, Clone)]
pub struct MaterialConfig {
    /// Per-CityObjectType base color as RGBA [f32; 4].
    pub color_map: HashMap<CityObjectType, [f32; 4]>,
    /// Global metallic factor for the GLB material (from [defaults]).
    pub metallic_factor: f32,
    /// Global roughness factor for the GLB material (from [defaults]).
    pub roughness_factor: f32,
}

fn clamp_factor(value: f32, name: &str, section: &str) -> f32 {
    if value < 0.0 || value > 1.0 {
        let clamped = value.clamp(0.0, 1.0);
        warn!(
            "[{}] {}: {} is out of range [0.0, 1.0], clamped to {}",
            section, name, value, clamped
        );
        clamped
    } else {
        value
    }
}

fn parse_city_object_type(name: &str) -> Option<CityObjectType> {
    match name {
        "Bridge" => Some(CityObjectType::Bridge),
        "BridgePart" => Some(CityObjectType::BridgePart),
        "BridgeInstallation" => Some(CityObjectType::BridgeInstallation),
        "BridgeConstructiveElement" => Some(CityObjectType::BridgeConstructiveElement),
        "BridgeRoom" => Some(CityObjectType::BridgeRoom),
        "BridgeFurniture" => Some(CityObjectType::BridgeFurniture),
        "Building" => Some(CityObjectType::Building),
        "BuildingPart" => Some(CityObjectType::BuildingPart),
        "BuildingInstallation" => Some(CityObjectType::BuildingInstallation),
        "BuildingConstructiveElement" => Some(CityObjectType::BuildingConstructiveElement),
        "BuildingFurniture" => Some(CityObjectType::BuildingFurniture),
        "BuildingStorey" => Some(CityObjectType::BuildingStorey),
        "BuildingRoom" => Some(CityObjectType::BuildingRoom),
        "BuildingUnit" => Some(CityObjectType::BuildingUnit),
        "CityFurniture" => Some(CityObjectType::CityFurniture),
        "LandUse" => Some(CityObjectType::LandUse),
        "OtherConstruction" => Some(CityObjectType::OtherConstruction),
        "PlantCover" => Some(CityObjectType::PlantCover),
        "SolitaryVegetationObject" => Some(CityObjectType::SolitaryVegetationObject),
        "TINRelief" => Some(CityObjectType::TINRelief),
        "Tunnel" => Some(CityObjectType::Tunnel),
        "TunnelPart" => Some(CityObjectType::TunnelPart),
        "TunnelInstallation" => Some(CityObjectType::TunnelInstallation),
        "WaterBody" => Some(CityObjectType::WaterBody),
        "Road" => Some(CityObjectType::Road),
        "Railway" => Some(CityObjectType::Railway),
        "Waterway" => Some(CityObjectType::Waterway),
        "TransportSquare" => Some(CityObjectType::TransportSquare),
        "GenericCityObject" => Some(CityObjectType::GenericCityObject),
        _ => None,
    }
}

impl MaterialConfig {
    /// Build from a TOML file path.
    pub fn from_toml(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read material config file: {:?}", path))?;
        let config_file: MaterialConfigFile = toml::from_str(&content)
            .with_context(|| format!("Failed to parse material config file: {:?}", path))?;

        // Resolve effective defaults: built-in <- [defaults] section
        let default_base_color = config_file
            .defaults
            .as_ref()
            .and_then(|d| d.base_color.as_deref())
            .unwrap_or(DEFAULT_BASE_COLOR);
        let default_rgba = hex_to_rgba(default_base_color)
            .context("Invalid base_color in [defaults] section")?;
        let global_metallic = config_file
            .defaults
            .as_ref()
            .and_then(|d| d.metallic_factor)
            .map(|v| clamp_factor(v, "metallic_factor", "defaults"))
            .unwrap_or(DEFAULT_METALLIC);
        let global_roughness = config_file
            .defaults
            .as_ref()
            .and_then(|d| d.roughness_factor)
            .map(|v| clamp_factor(v, "roughness_factor", "defaults"))
            .unwrap_or(DEFAULT_ROUGHNESS);

        // Build color map from per-type sections
        let mut color_map = HashMap::new();
        for (name, params) in &config_file.types {
            match parse_city_object_type(name) {
                Some(cotype) => {
                    let rgba = match &params.base_color {
                        Some(hex) => hex_to_rgba(hex)
                            .with_context(|| format!("Invalid base_color in [{}]", name))?,
                        None => default_rgba,
                    };
                    color_map.insert(cotype, rgba);
                }
                None => {
                    warn!(
                        "Unknown CityObjectType '{}' in material config, skipping",
                        name
                    );
                }
            }
        }

        // Fill in any types not explicitly listed with the effective default color
        for cotype in all_city_object_types() {
            color_map.entry(cotype).or_insert(default_rgba);
        }

        Ok(MaterialConfig {
            color_map,
            metallic_factor: global_metallic,
            roughness_factor: global_roughness,
        })
    }

    /// Build a default config matching current behavior (all pink, metallic=0.0, roughness≈0.502).
    pub fn default_config() -> Self {
        let default_rgba = hex_to_rgba(DEFAULT_BASE_COLOR).unwrap();
        let mut color_map = HashMap::new();
        for cotype in all_city_object_types() {
            color_map.insert(cotype, default_rgba);
        }
        MaterialConfig {
            color_map,
            metallic_factor: DEFAULT_METALLIC,
            roughness_factor: DEFAULT_ROUGHNESS,
        }
    }

    /// Apply CLI --glb-color-* overrides on top of current config.
    /// CLI flags take precedence over TOML per-type base_color.
    pub fn apply_cli_color_overrides(&mut self, cli_colors: &[(CityObjectType, &Option<String>)]) {
        for (cotype, opt_hex) in cli_colors {
            if let Some(hex) = opt_hex {
                // CLI hex values are already validated by clap's value_parser
                let rgba = hex_to_rgba(hex).unwrap();
                self.color_map.insert(*cotype, rgba);
            }
        }
    }
}

fn all_city_object_types() -> Vec<CityObjectType> {
    vec![
        CityObjectType::Bridge,
        CityObjectType::BridgePart,
        CityObjectType::BridgeInstallation,
        CityObjectType::BridgeConstructiveElement,
        CityObjectType::BridgeRoom,
        CityObjectType::BridgeFurniture,
        CityObjectType::Building,
        CityObjectType::BuildingPart,
        CityObjectType::BuildingInstallation,
        CityObjectType::BuildingConstructiveElement,
        CityObjectType::BuildingFurniture,
        CityObjectType::BuildingStorey,
        CityObjectType::BuildingRoom,
        CityObjectType::BuildingUnit,
        CityObjectType::CityFurniture,
        CityObjectType::LandUse,
        CityObjectType::OtherConstruction,
        CityObjectType::PlantCover,
        CityObjectType::SolitaryVegetationObject,
        CityObjectType::TINRelief,
        CityObjectType::Tunnel,
        CityObjectType::TunnelPart,
        CityObjectType::TunnelInstallation,
        CityObjectType::WaterBody,
        CityObjectType::Road,
        CityObjectType::Railway,
        CityObjectType::Waterway,
        CityObjectType::TransportSquare,
        CityObjectType::GenericCityObject,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hex_to_rgba() {
        let rgba = hex_to_rgba("#FF0000").unwrap();
        assert!((rgba[0] - 1.0).abs() < 0.01);
        assert!(rgba[1] < 0.01);
        assert!(rgba[2] < 0.01);
        assert!((rgba[3] - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_hex_to_rgba_invalid() {
        assert!(hex_to_rgba("FF0000").is_err());
        assert!(hex_to_rgba("#GGGGGG").is_err());
        assert!(hex_to_rgba("#FFF").is_err());
    }

    #[test]
    fn test_default_config_matches_current_behavior() {
        let config = MaterialConfig::default_config();
        let expected_pink = hex_to_rgba(DEFAULT_BASE_COLOR).unwrap();
        assert_eq!(config.color_map[&CityObjectType::Building], expected_pink);
        assert!((config.metallic_factor - 0.0).abs() < f32::EPSILON);
        assert!((config.roughness_factor - 128.0 / 255.0).abs() < 0.001);
    }

    #[test]
    fn test_parse_city_object_type() {
        assert_eq!(
            parse_city_object_type("Building"),
            Some(CityObjectType::Building)
        );
        assert_eq!(
            parse_city_object_type("GenericCityObject"),
            Some(CityObjectType::GenericCityObject)
        );
        assert_eq!(parse_city_object_type("Unknown"), None);
    }

    #[test]
    fn test_from_toml_basic() {
        let toml_str = r##"
[defaults]
base_color = "#808080"
metallic_factor = 0.1
roughness_factor = 0.9

[Building]
base_color = "#FF0000"

[WaterBody]
base_color = "#0000FF"
metallic_factor = 0.3
roughness_factor = 0.2
"##;
        let dir = std::env::temp_dir();
        let path = dir.join("test_material_config.toml");
        std::fs::write(&path, toml_str).unwrap();

        let config = MaterialConfig::from_toml(&path).unwrap();
        std::fs::remove_file(&path).ok();

        // Building gets red from TOML
        let building_color = config.color_map[&CityObjectType::Building];
        assert!((building_color[0] - 1.0).abs() < 0.01);
        assert!(building_color[1] < 0.01);

        // WaterBody gets blue from TOML
        let water_color = config.color_map[&CityObjectType::WaterBody];
        assert!(water_color[0] < 0.01);
        assert!((water_color[2] - 1.0).abs() < 0.01);

        // Unlisted types get the [defaults] gray
        let road_color = config.color_map[&CityObjectType::Road];
        let expected_gray = hex_to_rgba("#808080").unwrap();
        assert_eq!(road_color, expected_gray);

        // Global metallic/roughness from [defaults]
        assert!((config.metallic_factor - 0.1).abs() < 0.01);
        assert!((config.roughness_factor - 0.9).abs() < 0.01);
    }

    #[test]
    fn test_cli_overrides_toml() {
        let mut config = MaterialConfig::default_config();
        let red = Some("#FF0000".to_string());
        let none: Option<String> = None;
        let overrides: Vec<(CityObjectType, &Option<String>)> = vec![
            (CityObjectType::Building, &red),
            (CityObjectType::Road, &none),
        ];
        config.apply_cli_color_overrides(&overrides);

        // Building overridden to red
        let building_color = config.color_map[&CityObjectType::Building];
        assert!((building_color[0] - 1.0).abs() < 0.01);
        assert!(building_color[1] < 0.01);

        // Road unchanged (None override)
        let road_color = config.color_map[&CityObjectType::Road];
        let default_pink = hex_to_rgba(DEFAULT_BASE_COLOR).unwrap();
        assert_eq!(road_color, default_pink);
    }
}
