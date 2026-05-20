use std::collections::BTreeMap;
use std::io::Cursor;
use std::path::PathBuf;

use anyhow::{bail, Context};
use cityjson_convert::{
    convert_to_cityjson, convert_to_cityjsonseq, convert_to_glb, convert_to_obj,
    CityJsonSeqExportOptions, ExportOptions, GeometryPlacement, JsonExportOptions,
    ObjExportOptions,
};
use cityjson_lib::json;
use clap::{Parser, ValueEnum};
use log::info;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
#[clap(rename_all = "lower")]
enum OutputFormat {
    #[default]
    Glb,
    Obj,
    Cityjson,
    Cityjsonseq,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input `CityJSON` file (.city.json) or CityJSONSeq/CityJSONFeature stream.
    input: PathBuf,
    /// Path to the output file.
    #[arg(short, long)]
    output: PathBuf,
    /// Output format.
    #[arg(long, value_enum, default_value_t)]
    format: OutputFormat,
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

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.format {
        OutputFormat::Glb => {
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

            info!("Converting to GLB");
            convert_to_glb(&model, &cli.output, &options)?;
            info!("GLB written to {}", cli.output.display());
        }
        OutputFormat::Obj => {
            let model = json::from_file(&cli.input)?;
            info!("Converting to OBJ");
            convert_to_obj(&model, &cli.output, &ObjExportOptions::default())?;
            info!("OBJ written to {}", cli.output.display());
        }
        OutputFormat::Cityjson => {
            let model = json::from_file(&cli.input)?;
            info!("Converting to CityJSON");
            convert_to_cityjson(&model, &cli.output, &JsonExportOptions::default())?;
            info!("CityJSON written to {}", cli.output.display());
        }
        OutputFormat::Cityjsonseq => {
            let (base_root, feature_models) = read_cityjsonseq_or_feature_stream(&cli.input)?;
            if feature_models.is_empty() {
                bail!(
                    "--format cityjsonseq requires CityJSONSeq or CityJSONFeature-stream input; single CityJSON decomposition is not implemented"
                );
            }
            info!("Converting to CityJSONSeq");
            convert_to_cityjsonseq(
                &base_root,
                &feature_models,
                &cli.output,
                &CityJsonSeqExportOptions::default(),
            )?;
            info!("CityJSONSeq written to {}", cli.output.display());
        }
    }
    Ok(())
}

fn read_cityjsonseq_or_feature_stream(
    input: &PathBuf,
) -> anyhow::Result<(cityjson_lib::CityModel, Vec<cityjson_lib::CityModel>)> {
    let bytes = std::fs::read(input).with_context(|| format!("read {}", input.display()))?;
    let mut stream = serde_json::Deserializer::from_slice(&bytes).into_iter::<serde_json::Value>();
    let Some(first) = stream.next().transpose()? else {
        bail!("input stream is empty");
    };
    let first_bytes = serde_json::to_vec(&first)?;
    match json::probe(&first_bytes)?.kind() {
        cityjson_lib::json::RootKind::CityJSON => {
            if stream.next().is_none() {
                bail!(
                    "--format cityjsonseq requires CityJSONSeq or CityJSONFeature-stream input; single CityJSON decomposition is not implemented"
                );
            }
            let base_root = json::from_slice(&first_bytes)?;
            let feature_models =
                json::read_cityjsonseq(Cursor::new(bytes))?.collect::<Result<Vec<_>, _>>()?;
            Ok((base_root, feature_models))
        }
        cityjson_lib::json::RootKind::CityJSONFeature => {
            let mut feature_models = vec![json::from_feature_slice(&first_bytes)?];
            for item in stream {
                let item = item?;
                let item_bytes = serde_json::to_vec(&item)?;
                feature_models.push(json::from_feature_slice(&item_bytes)?);
            }
            let base_root = feature_models
                .first()
                .cloned()
                .ok_or_else(|| anyhow::anyhow!("input stream is empty"))?;
            Ok((base_root, feature_models))
        }
    }
}
