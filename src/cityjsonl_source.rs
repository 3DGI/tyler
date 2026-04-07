//! CityJSONL source — reads a CityJSONL file and splits it into a metadata file and
//! individual feature JSONL files compatible with Tyler's tiling pipeline.
//!
//! This is a Rust port of `roofer2tyler.py`. It works for any CityObject type, not just
//! buildings.

use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use log::info;

use crate::parser::{CityJSONFeatureVertices, Transform};

/// Metadata extracted from the first line of a CityJSONL file.
pub struct CityJsonlMetadata {
    /// The raw JSON string of the metadata line.
    #[allow(dead_code)]
    pub raw: String,
    /// The parsed transform (scale + translate).
    pub transform: Transform,
    /// The CRS string (e.g. "https://www.opengis.net/def/crs/EPSG/0/7415").
    pub reference_system: String,
}

/// Required keys on the metadata (first) line.
const REQUIRED_METADATA_KEYS: &[&str] = &[
    "type",
    "version",
    "transform",
    "metadata",
    "CityObjects",
    "vertices",
];

/// Required keys on each feature line.
#[allow(dead_code)]
const REQUIRED_FEATURE_KEYS: &[&str] = &["CityObjects", "id", "type", "vertices"];

/// Validate that `obj` contains all `required` keys.
fn validate_keys(
    obj: &serde_json::Map<String, serde_json::Value>,
    required: &[&str],
    label: &str,
) -> Result<()> {
    for &key in required {
        if !obj.contains_key(key) {
            bail!("{label} is missing required key '{key}'");
        }
    }
    Ok(())
}

/// Process a CityJSONL file: extract metadata and split features into individual files.
///
/// Returns the paths to the generated metadata file and features directory.
#[allow(dead_code)]
pub fn process_cityjsonl(
    jsonl_path: &Path,
    output_dir: &Path,
) -> Result<(PathBuf, PathBuf, CityJsonlMetadata)> {
    if !jsonl_path.exists() {
        bail!("Input file {} does not exist", jsonl_path.display());
    }

    let file =
        fs::File::open(jsonl_path).with_context(|| format!("opening {}", jsonl_path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // --- First line: metadata ---
    let first_line = lines
        .next()
        .context("CityJSONL file is empty; expected metadata line")?
        .context("reading metadata line")?;

    let metadata_value: serde_json::Value =
        serde_json::from_str(&first_line).context("parsing metadata JSON")?;
    let metadata_obj = metadata_value
        .as_object()
        .context("metadata line must be a JSON object")?;
    validate_keys(metadata_obj, REQUIRED_METADATA_KEYS, "Metadata")?;

    // Extract transform and CRS for callers that need them.
    let transform: Transform =
        serde_json::from_value(metadata_obj["transform"].clone()).context("parsing transform")?;

    let reference_system = metadata_obj
        .get("metadata")
        .and_then(|m| m.get("referenceSystem"))
        .and_then(|r| r.as_str())
        .unwrap_or("")
        .to_string();

    let metadata_path = output_dir.join("metadata.city.json");
    let features_dir = output_dir.join("features");

    // Clean up previous features directory if it exists.
    if features_dir.exists() {
        fs::remove_dir_all(&features_dir)
            .with_context(|| format!("removing existing {}", features_dir.display()))?;
    }
    fs::create_dir_all(&features_dir)
        .with_context(|| format!("creating {}", features_dir.display()))?;

    // --- Remaining lines: features ---
    let mut feature_count: usize = 0;
    for line_result in lines {
        let line = line_result.context("reading feature line")?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let value: serde_json::Value =
            serde_json::from_str(trimmed).context("parsing feature JSON")?;
        let obj = value
            .as_object()
            .context("feature line must be a JSON object")?;
        validate_keys(obj, REQUIRED_FEATURE_KEYS, "Feature")?;

        let feature_id = obj
            .get("id")
            .and_then(|v| v.as_str())
            .context("Feature 'id' must be a non-empty string")?;
        if feature_id.is_empty() {
            bail!("Feature 'id' must be a non-empty string");
        }

        let feature_path = features_dir.join(format!("{feature_id}.jsonl"));
        let mut f = fs::File::create(&feature_path)
            .with_context(|| format!("creating {}", feature_path.display()))?;
        writeln!(f, "{trimmed}").context("writing feature")?;

        feature_count += 1;
    }

    if feature_count == 0 {
        bail!(
            "{} does not contain any CityObject feature records",
            jsonl_path.display()
        );
    }

    // Write the metadata file (pretty-printed, matching Python behaviour).
    let pretty =
        serde_json::to_string_pretty(&metadata_value).context("serializing metadata JSON")?;
    fs::write(&metadata_path, pretty)
        .with_context(|| format!("writing {}", metadata_path.display()))?;

    info!(
        "Wrote {} features to {}",
        feature_count,
        features_dir.display()
    );

    let meta = CityJsonlMetadata {
        raw: first_line,
        transform,
        reference_system,
    };

    Ok((metadata_path, features_dir, meta))
}

/// Load a CityJSONL file entirely into memory: metadata + all features parsed.
///
/// Unlike [`process_cityjsonl`], this never writes individual files to disk.
/// Each feature line is parsed directly into [`CityJSONFeatureVertices`].
pub fn load_cityjsonl_to_memory(
    jsonl_path: &Path,
) -> Result<(CityJsonlMetadata, Vec<CityJSONFeatureVertices>)> {
    if !jsonl_path.exists() {
        bail!("Input file {} does not exist", jsonl_path.display());
    }

    let file_size = fs::metadata(jsonl_path)
        .with_context(|| format!("stat {}", jsonl_path.display()))?
        .len() as usize;

    let file =
        fs::File::open(jsonl_path).with_context(|| format!("opening {}", jsonl_path.display()))?;
    let reader = BufReader::new(file);

    // --- First line: metadata ---
    let mut lines = reader.lines();
    let first_line = lines
        .next()
        .context("CityJSONL file is empty; expected metadata line")?
        .context("reading metadata line")?;

    let metadata_value: serde_json::Value =
        serde_json::from_str(&first_line).context("parsing metadata JSON")?;
    let metadata_obj = metadata_value
        .as_object()
        .context("metadata line must be a JSON object")?;
    validate_keys(metadata_obj, REQUIRED_METADATA_KEYS, "Metadata")?;

    let transform: Transform =
        serde_json::from_value(metadata_obj["transform"].clone()).context("parsing transform")?;

    let reference_system = metadata_obj
        .get("metadata")
        .and_then(|m| m.get("referenceSystem"))
        .and_then(|r| r.as_str())
        .unwrap_or("")
        .to_string();

    // --- Remaining lines: features (parsed into memory) ---
    // Estimate capacity: average feature ~5-8KB
    let mut features: Vec<CityJSONFeatureVertices> =
        Vec::with_capacity(file_size / 5000);

    for line_result in lines {
        let line = line_result.context("reading feature line")?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let cf: CityJSONFeatureVertices =
            serde_json::from_str(trimmed).context("parsing feature into CityJSONFeatureVertices")?;
        features.push(cf);
    }

    if features.is_empty() {
        bail!(
            "{} does not contain any CityObject feature records",
            jsonl_path.display()
        );
    }

    info!(
        "Loaded {} features into memory from {}",
        features.len(),
        jsonl_path.display()
    );

    let meta = CityJsonlMetadata {
        raw: first_line,
        transform,
        reference_system,
    };

    Ok((meta, features))
}
