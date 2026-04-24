use std::collections::BTreeMap;
use std::path::PathBuf;

use cityjson_lib::json;
use clap::Parser;

use cityjson_convert::{convert_to_glb, ExportOptions, GeometryPlacement};

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input `CityJSON` file (.city.json).
    input: PathBuf,
    /// Path to the output GLB file.
    #[arg(short, long)]
    output: PathBuf,
    /// Default PBR base color for the generated GLB.
    #[arg(long = "native-glb-color", default_value = "#FFC0CB")]
    native_glb_color: String,
    /// Disable geometry quantization in the generated GLB.
    #[arg(long = "no-quantization", default_value_t = false)]
    no_quantization: bool,
    /// Share vertices and average incident face normals.
    #[arg(long = "smooth-normals", default_value_t = false)]
    smooth_normals: bool,
    /// Disable `EXT_meshopt_compression` in the generated GLB.
    #[arg(long = "no-meshopt-compression", default_value_t = false)]
    no_meshopt_compression: bool,
    /// Metadata class name for `EXT_structural_metadata`.
    #[arg(long = "3dtiles-metadata-class", default_value = "cityobject")]
    metadata_class_name: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    let model = json::from_file(&cli.input)?;
    let options = ExportOptions {
        native_glb_color: cli.native_glb_color,
        metadata_class_name: cli.metadata_class_name,
        feature_type_colors: BTreeMap::default(),
        geometry_placement: GeometryPlacement::SourceCoordinates,
        clip_bbox: None,
        clip_geographic_region: None,
        smooth_normals: cli.smooth_normals,
        quantize_geometry: !cli.no_quantization,
        meshopt_compression: !cli.no_meshopt_compression,
    };

    convert_to_glb(&model, cli.output, &options)?;
    Ok(())
}
