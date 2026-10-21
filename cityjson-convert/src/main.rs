use std::collections::BTreeMap;
use std::io::Cursor;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context};
use cityjson_convert::{
    convert_to_cityjson, convert_to_cityjsonseq, convert_to_glb, convert_to_obj, convert_to_tsv,
    CityJsonSeqExportOptions, ExportOptions, GeometryPlacement, GpkgExportOptions,
    JsonExportOptions, ObjExportOptions, TsvExportOptions,
};
use cityjson_lib::json;
use clap::{Args, Parser, ValueEnum};
use log::info;

const CJINDEX_PAGE_SIZE: usize = 65_536;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
#[clap(rename_all = "lower")]
enum OutputFormat {
    #[default]
    Glb,
    Obj,
    Cityjson,
    Cityjsonseq,
    Tsv,
    Gpkg,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input `CityJSON` file, `CityJSONSeq`/`CityJSONFeature` stream, or dataset directory.
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
    // GeoPackage options.
    #[command(flatten)]
    gpkg: GpkgCliOptions,
    /// TSV row-shape options.
    #[command(flatten)]
    tsv_rows: TsvRowCliOptions,
    /// TSV metadata and semantics options.
    #[command(flatten)]
    tsv_files: TsvFileCliOptions,
}

#[derive(Args, Debug, Default)]
struct TsvRowCliOptions {
    /// Include rows whose dynamic TSV attributes are all null.
    #[arg(long = "tsv-include-null-rows", default_value_t = false)]
    tsv_include_null_rows: bool,
    /// Include parent/child hierarchy columns in TSV outputs.
    #[arg(long = "tsv-include-hierarchy", default_value_t = false)]
    tsv_include_hierarchy: bool,
    /// Include the source `CityJSON` ordinal column in TSV outputs.
    #[arg(long = "tsv-include-cityjson-ordinal", default_value_t = false)]
    tsv_include_cityjson_ordinal: bool,
}

#[derive(Args, Debug, Default)]
struct TsvFileCliOptions {
    /// Write a separate TSV metadata file.
    #[arg(long = "tsv-include-metadata", default_value_t = false)]
    tsv_include_metadata: bool,
    /// Write a separate TSV semantics file joined from primitive assignments.
    #[arg(long = "tsv-include-semantics", default_value_t = false)]
    tsv_include_semantics: bool,
    /// Write `CityObject` extra.address values to addresses.tsv.
    #[arg(long = "tsv-include-address", default_value_t = false)]
    tsv_include_address: bool,
}

#[derive(Args, Debug, Default)]
struct GpkgCliOptions {
    #[command(flatten)]
    semantics: GpkgSemanticCliOptions,
    #[command(flatten)]
    address: GpkgAddressCliOptions,
    #[command(flatten)]
    metadata: GpkgMetadataCliOptions,
    /// Split feature layers by `CityObject` type, geometry family, and `CityJSON` `LoD`.
    #[arg(long = "gpkg-split-lod", default_value_t = false)]
    split_lod: bool,
}

#[derive(Args, Debug, Default)]
struct GpkgSemanticCliOptions {
    /// Add a `semantics` feature layer containing semantic primitive rows.
    #[arg(long = "gpkg-include-semantics", default_value_t = false)]
    include_semantics: bool,
    /// Add the non-spatial `cityobject_hierarchy` and `semantic_hierarchy` tables.
    #[arg(long = "gpkg-include-hierarchy", default_value_t = false)]
    include_hierarchy: bool,
}

#[derive(Args, Debug, Default)]
struct GpkgAddressCliOptions {
    /// Add an `addresses` feature layer from `CityObject.extra.address` values.
    #[arg(long = "gpkg-include-address", default_value_t = false)]
    include_address: bool,
}

#[derive(Args, Debug, Default)]
struct GpkgMetadataCliOptions {
    /// Export source `CityJSON` metadata through the `GeoPackage` metadata extension.
    #[arg(long = "gpkg-include-metadata", default_value_t = false)]
    include_source: bool,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.format {
        OutputFormat::Glb => {
            let model = read_model(&cli.input)?;
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
            let model = read_model(&cli.input)?;
            info!("Converting to OBJ");
            convert_to_obj(&model, &cli.output, &ObjExportOptions::default())?;
            info!("OBJ written to {}", cli.output.display());
        }
        OutputFormat::Cityjson => {
            let model = read_model(&cli.input)?;
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
        OutputFormat::Gpkg => {
            let model = read_model(&cli.input)?;
            let options = GpkgExportOptions {
                include_semantics: cli.gpkg.semantics.include_semantics,
                include_address: cli.gpkg.address.include_address,
                include_hierarchy: cli.gpkg.semantics.include_hierarchy,
                include_metadata: cli.gpkg.metadata.include_source,
                split_lod: cli.gpkg.split_lod,
            };
            info!("Converting to GeoPackage");
            cityjson_convert::convert_to_gpkg(&model, &cli.output, &options)?;
            info!("GeoPackage written to {}", cli.output.display());
        }
        OutputFormat::Tsv => {
            let model = json::from_file(&cli.input)?;
            let options = TsvExportOptions {
                output_suffix: None,
                include_null_rows: cli.tsv_rows.tsv_include_null_rows,
                include_hierarchy: cli.tsv_rows.tsv_include_hierarchy,
                include_cityjson_ordinal: cli.tsv_rows.tsv_include_cityjson_ordinal,
                include_metadata: cli.tsv_files.tsv_include_metadata,
                include_semantics: cli.tsv_files.tsv_include_semantics,
                include_address: cli.tsv_files.tsv_include_address,
            };
            info!("Converting to TSV");
            convert_to_tsv(&model, &cli.output, &options)?;
            info!("TSV files written to {}", cli.output.display());
        }
    }
    Ok(())
}

fn read_model(input: &Path) -> anyhow::Result<cityjson_lib::CityModel> {
    if !input.is_dir() {
        return json::from_file(input).map_err(Into::into);
    }

    let (_, feature_models) = read_indexed_dataset(input)?;
    cityjson_lib::ops::merge(feature_models).map_err(Into::into)
}

fn read_cityjsonseq_or_feature_stream(
    input: &Path,
) -> anyhow::Result<(cityjson_lib::CityModel, Vec<cityjson_lib::CityModel>)> {
    if input.is_dir() {
        return read_indexed_dataset(input);
    }

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

fn read_indexed_dataset(
    input: &Path,
) -> anyhow::Result<(cityjson_lib::CityModel, Vec<cityjson_lib::CityModel>)> {
    let resolved = cityjson_index::resolve_dataset(input, None)
        .with_context(|| format!("resolve dataset {}", input.display()))?;
    let inspection = resolved.inspect()?;
    let mut city_index =
        cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)?;
    if !inspection.index.exists || inspection.index.fresh != Some(true) {
        info!(
            "Rebuilding cjindex sidecar at {}",
            resolved.index_path.display()
        );
        city_index.reindex()?;
    }

    let metadata = city_index.metadata()?;
    let base_document = metadata
        .first()
        .ok_or_else(|| anyhow::anyhow!("dataset does not contain any source metadata"))?;
    let base_root = json::from_slice(&serde_json::to_vec(base_document.as_ref())?)?;
    let mut feature_models = Vec::new();
    let mut after_record_id = None;
    loop {
        let page =
            city_index.package_ref_page_after_record_id(after_record_id, CJINDEX_PAGE_SIZE)?;
        if page.is_empty() {
            break;
        }
        after_record_id = page.last().map(|package| package.record_id);
        feature_models.extend(
            city_index
                .read_packages(&page)?
                .into_iter()
                .map(|package| package.model),
        );
    }
    if feature_models.is_empty() {
        bail!("dataset does not contain any features");
    }

    Ok((base_root, feature_models))
}
