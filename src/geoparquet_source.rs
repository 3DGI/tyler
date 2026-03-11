//! GeoParquet source — reads 3D polygon geometries from a GeoParquet file and converts
//! them into CityJSONFeature JSONL files.
//!
//! Uses Apache `parquet` crate for Parquet I/O and `geozero` for WKB geometry parsing.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use arrow::array::{Array, BinaryArray, LargeBinaryArray, StringArray};
use arrow::datatypes::DataType;
use log::{debug, info, warn};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use serde_json::{json, Map, Value};

use crate::parser::Transform;
use crate::transform_align::TransformAligner;

/// A single 3D polygon ring (closed ring of [x, y, z] vertices).
type Ring3D = Vec<[f64; 3]>;

/// A polygon: outer ring + optional inner rings.
type Polygon3D = Vec<Ring3D>;

/// A tree (or generic feature) parsed from GeoParquet.
struct GeoFeature {
    id: String,
    polygons: Vec<Polygon3D>,
    attributes: Map<String, Value>,
}

/// Extract the EPSG code from GeoParquet file-level metadata.
///
/// GeoParquet stores a JSON object under the `geo` key in Parquet's key-value metadata.
/// The CRS is in `columns.<geom_col>.crs` as PROJJSON containing `id.code`.
fn extract_epsg_from_geo_metadata(geo_json: &str) -> Result<u16> {
    let geo: Value = serde_json::from_str(geo_json).context("parsing 'geo' metadata JSON")?;

    // Find the primary geometry column name.
    let primary_column = geo
        .get("primary_column")
        .and_then(|v| v.as_str())
        .unwrap_or("geometry");

    // Navigate to the CRS definition.
    let crs = geo
        .get("columns")
        .and_then(|c| c.get(primary_column))
        .and_then(|c| c.get("crs"));

    let crs = match crs {
        Some(c) => c,
        None => bail!("No CRS found in GeoParquet metadata for column '{primary_column}'"),
    };

    // PROJJSON: look for id.authority + id.code
    if let Some(id) = crs.get("id") {
        if let (Some(auth), Some(code)) = (
            id.get("authority").and_then(|a| a.as_str()),
            id.get("code").and_then(|c| c.as_u64()),
        ) {
            if auth == "EPSG" {
                return Ok(code as u16);
            }
        }
    }

    // Fallback: look in properties → identifiers
    if let Some(ids) = crs.get("$schema").and_then(|_| crs.get("id")) {
        if let Some(code) = ids.get("code").and_then(|c| c.as_u64()) {
            return Ok(code as u16);
        }
    }

    bail!("Could not extract EPSG code from GeoParquet CRS metadata")
}

/// Extract the geometry column name from GeoParquet metadata.
fn geometry_column_name(geo_json: &str) -> Result<String> {
    let geo: Value = serde_json::from_str(geo_json).context("parsing 'geo' metadata")?;
    Ok(geo
        .get("primary_column")
        .and_then(|v| v.as_str())
        .unwrap_or("geometry")
        .to_string())
}

/// Parse WKB bytes into 3D polygons.
///
/// Supports PolygonZ (WKB type 1003 / 0x80000003) and MultiPolygonZ (1006 / 0x80000006).
/// Uses manual WKB parsing for zero-copy performance — geozero's `GeomProcessor` trait
/// adds abstraction we don't need for this simple case.
fn parse_wkb_3d(wkb: &[u8]) -> Result<Vec<Polygon3D>> {
    if wkb.len() < 5 {
        bail!("WKB too short: {} bytes", wkb.len());
    }

    let is_little_endian = wkb[0] == 1;
    let geom_type = read_u32(wkb, 1, is_little_endian);

    // Mask out EWKB SRID flag and Z/M flags to get the base + Z type
    let base_type = geom_type & 0x0000FFFF;
    let has_z = (geom_type & 0x80000000) != 0 || base_type >= 1000;

    // Normalize to ISO type numbers
    let iso_type = if base_type >= 1000 {
        base_type - 1000
    } else {
        base_type
    };

    match iso_type {
        3 => {
            // PolygonZ
            if !has_z && base_type == 3 {
                // Plain 2D polygon, treat Z as 0
            }
            let (polygon, _) = parse_polygon(wkb, 5, is_little_endian, has_z)?;
            Ok(vec![polygon])
        }
        6 => {
            // MultiPolygonZ
            let num_polygons = read_u32(wkb, 5, is_little_endian) as usize;
            let mut offset = 9;
            let mut polygons = Vec::with_capacity(num_polygons);
            for _ in 0..num_polygons {
                if offset + 5 > wkb.len() {
                    bail!("WKB truncated in MultiPolygon");
                }
                let sub_le = wkb[offset] == 1;
                // Skip byte-order + type (5 bytes)
                offset += 5;
                let (polygon, new_offset) = parse_polygon(wkb, offset, sub_le, has_z)?;
                polygons.push(polygon);
                offset = new_offset;
            }
            Ok(polygons)
        }
        _ => bail!("Unsupported WKB geometry type: {geom_type} (iso: {iso_type})"),
    }
}

fn parse_polygon(
    wkb: &[u8],
    start: usize,
    le: bool,
    has_z: bool,
) -> Result<(Polygon3D, usize)> {
    let num_rings = read_u32(wkb, start, le) as usize;
    let mut offset = start + 4;
    let mut rings = Vec::with_capacity(num_rings);
    let coord_size = if has_z { 3 } else { 2 };
    let bytes_per_coord = coord_size * 8;

    for _ in 0..num_rings {
        if offset + 4 > wkb.len() {
            bail!("WKB truncated reading ring point count");
        }
        let num_points = read_u32(wkb, offset, le) as usize;
        offset += 4;
        let needed = num_points * bytes_per_coord;
        if offset + needed > wkb.len() {
            bail!("WKB truncated reading ring coordinates");
        }
        let mut ring = Vec::with_capacity(num_points);
        for _ in 0..num_points {
            let x = read_f64(wkb, offset, le);
            let y = read_f64(wkb, offset + 8, le);
            let z = if has_z {
                read_f64(wkb, offset + 16, le)
            } else {
                0.0
            };
            ring.push([x, y, z]);
            offset += bytes_per_coord;
        }
        rings.push(ring);
    }

    Ok((rings, offset))
}

#[inline]
fn read_u32(buf: &[u8], offset: usize, le: bool) -> u32 {
    let bytes: [u8; 4] = buf[offset..offset + 4].try_into().unwrap();
    if le {
        u32::from_le_bytes(bytes)
    } else {
        u32::from_be_bytes(bytes)
    }
}

#[inline]
fn read_f64(buf: &[u8], offset: usize, le: bool) -> f64 {
    let bytes: [u8; 8] = buf[offset..offset + 8].try_into().unwrap();
    if le {
        f64::from_le_bytes(bytes)
    } else {
        f64::from_be_bytes(bytes)
    }
}

/// Compute a simple centroid from the first polygon's outer ring.
#[allow(dead_code)]
fn centroid_2d(polygons: &[Polygon3D]) -> Option<[f64; 2]> {
    let ring = polygons.first()?.first()?;
    if ring.is_empty() {
        return None;
    }
    let (sx, sy) = ring
        .iter()
        .fold((0.0, 0.0), |(ax, ay), v| (ax + v[0], ay + v[1]));
    let n = ring.len() as f64;
    Some([sx / n, sy / n])
}

/// Convert 3D polygons into CityJSONFeature JSON and write to `features_dir`.
///
/// The `aligner` handles CRS reprojection (if needed) and quantization.
fn write_tree_feature(
    feature: &GeoFeature,
    aligner: &TransformAligner,
    city_object_type: &str,
    features_dir: &Path,
) -> Result<()> {
    // Collect all unique vertices and build boundary indices.
    let mut vertices_qc: Vec<[i64; 3]> = Vec::new();
    let mut boundaries: Vec<Value> = Vec::new();

    for polygon in &feature.polygons {
        for ring in polygon {
            let mut ring_indices: Vec<Value> = Vec::with_capacity(ring.len());
            for vertex in ring {
                let qc = aligner.align_and_quantize(vertex)?;
                let idx = vertices_qc.len();
                vertices_qc.push(qc);
                ring_indices.push(json!(idx));
            }
            // CityJSON boundary: surface = [ring], so we wrap ring_indices
            boundaries.push(json!([ring_indices]));
        }
    }

    // Compute height from Z range of quantized vertices.
    let transform = aligner.ref_transform();
    let z_values: Vec<f64> = vertices_qc
        .iter()
        .map(|v| v[2] as f64 * transform.scale[2] + transform.translate[2])
        .collect();
    let height = z_values
        .iter()
        .cloned()
        .reduce(f64::max)
        .unwrap_or(0.0)
        - z_values
            .iter()
            .cloned()
            .reduce(f64::min)
            .unwrap_or(0.0);

    // Build attributes: include height + any user-supplied attributes from the Parquet row.
    let mut attributes = feature.attributes.clone();
    attributes.insert("height".to_string(), json!(height));

    let city_object = json!({
        "type": city_object_type,
        "attributes": attributes,
        "geometry": [{
            "type": "MultiSurface",
            "lod": "2",
            "boundaries": boundaries,
        }],
    });

    let cityjson_feature = json!({
        "type": "CityJSONFeature",
        "id": feature.id,
        "CityObjects": {
            &feature.id: city_object,
        },
        "vertices": vertices_qc,
    });

    let feature_path = features_dir.join(format!("{}.jsonl", feature.id));
    let mut f =
        fs::File::create(&feature_path).with_context(|| format!("creating {}", feature_path.display()))?;
    let line = serde_json::to_string(&cityjson_feature).context("serializing CityJSONFeature")?;
    writeln!(f, "{line}").context("writing feature file")?;

    Ok(())
}

/// Read a GeoParquet file and write CityJSONFeature JSONL files into `features_dir`.
///
/// * `parquet_path` — Path to the GeoParquet file.
/// * `features_dir` — Existing directory where feature JSONL files will be written.
/// * `ref_transform` — The reference transform to quantize vertices into.
/// * `ref_epsg` — The EPSG code of the reference CRS.
/// * `city_object_type` — The CityJSON object type (e.g. "SolitaryVegetationObject").
/// * `id_column` — Optional column name to use as feature ID. Defaults to row index.
pub fn process_geoparquet(
    parquet_path: &Path,
    features_dir: &Path,
    ref_transform: Transform,
    ref_epsg: u16,
    city_object_type: &str,
    id_column: Option<&str>,
) -> Result<usize> {
    if !parquet_path.exists() {
        bail!(
            "GeoParquet file does not exist: {}",
            parquet_path.display()
        );
    }

    let file = fs::File::open(parquet_path)
        .with_context(|| format!("opening {}", parquet_path.display()))?;

    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .context("building Parquet reader")?;

    // Extract GeoParquet metadata.
    let kv_metadata = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .context("Parquet file has no key-value metadata")?;

    let geo_meta = kv_metadata
        .iter()
        .find(|kv| kv.key == "geo")
        .and_then(|kv| kv.value.as_ref())
        .context("No 'geo' key in Parquet metadata — is this a GeoParquet file?")?;

    let src_epsg = extract_epsg_from_geo_metadata(geo_meta)?;
    let geom_col = geometry_column_name(geo_meta)?;
    debug!(
        "GeoParquet: source EPSG:{}, geometry column: '{}'",
        src_epsg, geom_col
    );

    let aligner = TransformAligner::new(src_epsg, ref_epsg, ref_transform)?;

    let schema = builder.schema().clone();
    let reader = builder.build().context("building Parquet batch reader")?;

    // Find the geometry column index.
    let geom_idx = schema
        .index_of(&geom_col)
        .with_context(|| format!("geometry column '{geom_col}' not found in schema"))?;

    // Find optional ID column index.
    let id_col_idx = id_column.and_then(|name| schema.index_of(name).ok());

    // Find attribute columns (everything except geometry).
    let attr_indices: Vec<(usize, String)> = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != geom_idx)
        .map(|(i, f)| (i, f.name().clone()))
        .collect();

    let mut feature_count: usize = 0;
    let mut row_idx: usize = 0;

    for batch_result in reader {
        let batch = batch_result.context("reading Parquet record batch")?;
        let geom_array = batch.column(geom_idx);
        let num_rows = batch.num_rows();

        for i in 0..num_rows {
            // Read geometry WKB.
            let wkb_bytes = extract_binary(geom_array, i)?;
            if wkb_bytes.is_empty() {
                warn!("Row {row_idx}: empty geometry, skipping");
                row_idx += 1;
                continue;
            }

            // Parse WKB into 3D polygons.
            let polygons = match parse_wkb_3d(wkb_bytes) {
                Ok(p) => p,
                Err(e) => {
                    warn!("Row {row_idx}: failed to parse WKB: {e}, skipping");
                    row_idx += 1;
                    continue;
                }
            };

            // Determine feature ID.
            let id = if let Some(col_idx) = id_col_idx {
                extract_string(batch.column(col_idx), i)
                    .unwrap_or_else(|| format!("tree-{row_idx:05}"))
            } else {
                format!("tree-{row_idx:05}")
            };

            // Collect attribute values.
            let mut attributes = Map::new();
            for (col_idx, col_name) in &attr_indices {
                if id_col_idx == Some(*col_idx) {
                    continue; // Don't duplicate the ID column as an attribute.
                }
                if let Some(val) = extract_value(batch.column(*col_idx), i) {
                    attributes.insert(col_name.clone(), val);
                }
            }

            let feature = GeoFeature {
                id,
                polygons,
                attributes,
            };

            write_tree_feature(&feature, &aligner, city_object_type, features_dir)?;
            feature_count += 1;
            row_idx += 1;
        }
    }

    info!(
        "Wrote {} {} features from {}",
        feature_count,
        city_object_type,
        parquet_path.display()
    );

    Ok(feature_count)
}

/// Process a GeoParquet file in standalone mode (no buildings reference).
///
/// Reads the GeoParquet, computes the data's bounding box to derive a CityJSON
/// transform, writes feature files and a `metadata.city.json`.
///
/// Returns `(metadata_path, features_dir)`.
pub fn process_geoparquet_standalone(
    parquet_path: &Path,
    output_dir: &Path,
    city_object_type: &str,
    id_column: Option<&str>,
) -> Result<(PathBuf, PathBuf)> {
    use crate::parser::Transform;
    use crate::transform_align::quantize_vertex;

    if !parquet_path.exists() {
        bail!(
            "GeoParquet file does not exist: {}",
            parquet_path.display()
        );
    }

    let file = fs::File::open(parquet_path)
        .with_context(|| format!("opening {}", parquet_path.display()))?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).context("building Parquet reader")?;

    let kv_metadata = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .context("Parquet file has no key-value metadata")?;
    let geo_meta = kv_metadata
        .iter()
        .find(|kv| kv.key == "geo")
        .and_then(|kv| kv.value.as_ref())
        .context("No 'geo' key in Parquet metadata — is this a GeoParquet file?")?;

    let src_epsg = extract_epsg_from_geo_metadata(geo_meta)?;
    let geom_col = geometry_column_name(geo_meta)?;
    info!(
        "Standalone GeoParquet: EPSG:{}, geometry column: '{}'",
        src_epsg, geom_col
    );

    let schema = builder.schema().clone();
    let reader = builder.build().context("building Parquet batch reader")?;

    let geom_idx = schema
        .index_of(&geom_col)
        .with_context(|| format!("geometry column '{geom_col}' not found in schema"))?;
    let id_col_idx = id_column.and_then(|name| schema.index_of(name).ok());
    let attr_indices: Vec<(usize, String)> = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != geom_idx)
        .map(|(i, f)| (i, f.name().clone()))
        .collect();

    // First pass: parse all features and compute bounding box.
    let mut features: Vec<GeoFeature> = Vec::new();
    let mut bbox_min = [f64::INFINITY; 3];
    let mut bbox_max = [f64::NEG_INFINITY; 3];
    let mut row_idx: usize = 0;

    for batch_result in reader {
        let batch = batch_result.context("reading Parquet record batch")?;
        let geom_array = batch.column(geom_idx);
        let num_rows = batch.num_rows();

        for i in 0..num_rows {
            let wkb_bytes = extract_binary(geom_array, i)?;
            if wkb_bytes.is_empty() {
                warn!("Row {row_idx}: empty geometry, skipping");
                row_idx += 1;
                continue;
            }
            let polygons = match parse_wkb_3d(wkb_bytes) {
                Ok(p) => p,
                Err(e) => {
                    warn!("Row {row_idx}: failed to parse WKB: {e}, skipping");
                    row_idx += 1;
                    continue;
                }
            };
            // Update bounding box.
            for polygon in &polygons {
                for ring in polygon {
                    for v in ring {
                        bbox_min[0] = bbox_min[0].min(v[0]);
                        bbox_min[1] = bbox_min[1].min(v[1]);
                        bbox_min[2] = bbox_min[2].min(v[2]);
                        bbox_max[0] = bbox_max[0].max(v[0]);
                        bbox_max[1] = bbox_max[1].max(v[1]);
                        bbox_max[2] = bbox_max[2].max(v[2]);
                    }
                }
            }
            let id = if let Some(col_idx) = id_col_idx {
                extract_string(batch.column(col_idx), i)
                    .unwrap_or_else(|| format!("tree-{row_idx:05}"))
            } else {
                format!("tree-{row_idx:05}")
            };
            let mut attributes = Map::new();
            for (col_idx, col_name) in &attr_indices {
                if id_col_idx == Some(*col_idx) {
                    continue;
                }
                if let Some(val) = extract_value(batch.column(*col_idx), i) {
                    attributes.insert(col_name.clone(), val);
                }
            }
            features.push(GeoFeature {
                id,
                polygons,
                attributes,
            });
            row_idx += 1;
        }
    }

    if features.is_empty() {
        bail!(
            "No valid features found in {}",
            parquet_path.display()
        );
    }

    // Derive transform: millimetre precision, translate to bbox min.
    let transform = Transform {
        scale: [0.001, 0.001, 0.001],
        translate: [bbox_min[0].floor(), bbox_min[1].floor(), bbox_min[2].floor()],
    };

    let features_dir = output_dir.join("features");
    if features_dir.exists() {
        fs::remove_dir_all(&features_dir)?;
    }
    fs::create_dir_all(&features_dir)?;

    // Second pass: write feature files using the computed transform.
    for feature in &features {
        let mut vertices_qc: Vec<[i64; 3]> = Vec::new();
        let mut boundaries: Vec<Value> = Vec::new();

        for polygon in &feature.polygons {
            for ring in polygon {
                let mut ring_indices: Vec<Value> = Vec::with_capacity(ring.len());
                for vertex in ring {
                    let qc = quantize_vertex(vertex, &transform);
                    let idx = vertices_qc.len();
                    vertices_qc.push(qc);
                    ring_indices.push(json!(idx));
                }
                boundaries.push(json!([ring_indices]));
            }
        }

        let z_values: Vec<f64> = vertices_qc
            .iter()
            .map(|v| v[2] as f64 * transform.scale[2] + transform.translate[2])
            .collect();
        let height = z_values.iter().cloned().reduce(f64::max).unwrap_or(0.0)
            - z_values.iter().cloned().reduce(f64::min).unwrap_or(0.0);

        let mut attributes = feature.attributes.clone();
        attributes.insert("height".to_string(), json!(height));

        let city_object = json!({
            "type": city_object_type,
            "attributes": attributes,
            "geometry": [{
                "type": "MultiSurface",
                "lod": "2",
                "boundaries": boundaries,
            }],
        });
        let cityjson_feature = json!({
            "type": "CityJSONFeature",
            "id": feature.id,
            "CityObjects": { &feature.id: city_object },
            "vertices": vertices_qc,
        });

        let feature_path = features_dir.join(format!("{}.jsonl", feature.id));
        let mut f = fs::File::create(&feature_path)
            .with_context(|| format!("creating {}", feature_path.display()))?;
        let line = serde_json::to_string(&cityjson_feature).context("serializing feature")?;
        writeln!(f, "{line}").context("writing feature file")?;
    }

    // Write metadata.city.json
    let ref_system = format!("https://www.opengis.net/def/crs/EPSG/0/{src_epsg}");
    let metadata_value = json!({
        "type": "CityJSON",
        "version": "2.0",
        "transform": {
            "scale": transform.scale,
            "translate": transform.translate,
        },
        "metadata": {
            "referenceSystem": ref_system,
        },
        "CityObjects": {},
        "vertices": [],
    });
    let metadata_path = output_dir.join("metadata.city.json");
    let pretty = serde_json::to_string_pretty(&metadata_value)?;
    fs::write(&metadata_path, &pretty)?;

    info!(
        "Wrote {} {} features from {} (standalone)",
        features.len(),
        city_object_type,
        parquet_path.display()
    );

    Ok((metadata_path, features_dir))
}

/// Check if a tree centroid overlaps any building bounding box in the features directory.
/// Returns the number of overlaps detected (as warnings).
pub fn check_overlap(features_dir: &Path, parquet_path: &Path) -> Result<usize> {
    // Collect building bboxes from existing features.
    let mut building_bboxes: Vec<[f64; 4]> = Vec::new(); // [xmin, ymin, xmax, ymax]

    for entry in fs::read_dir(features_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(true, |e| e != "jsonl") {
            continue;
        }
        let content = fs::read_to_string(&path)?;
        let value: Value = match serde_json::from_str(content.trim()) {
            Ok(v) => v,
            Err(_) => continue,
        };
        // Check if this is a Building feature.
        if let Some(objs) = value.get("CityObjects").and_then(|v| v.as_object()) {
            for (_key, obj) in objs {
                if obj.get("type").and_then(|t| t.as_str()) == Some("Building") {
                    // Compute bbox from vertices.
                    if let Some(vertices) = value.get("vertices").and_then(|v| v.as_array()) {
                        let xs: Vec<f64> = vertices
                            .iter()
                            .filter_map(|v| v.as_array())
                            .filter_map(|a| a.first())
                            .filter_map(|v| v.as_f64())
                            .collect();
                        let ys: Vec<f64> = vertices
                            .iter()
                            .filter_map(|v| v.as_array())
                            .filter_map(|a| a.get(1))
                            .filter_map(|v| v.as_f64())
                            .collect();
                        if let (Some(&xmin), Some(&xmax), Some(&ymin), Some(&ymax)) = (
                            xs.iter().reduce(|a, b| if a < b { a } else { b }),
                            xs.iter().reduce(|a, b| if a > b { a } else { b }),
                            ys.iter().reduce(|a, b| if a < b { a } else { b }),
                            ys.iter().reduce(|a, b| if a > b { a } else { b }),
                        ) {
                            building_bboxes.push([xmin, ymin, xmax, ymax]);
                        }
                    }
                }
            }
        }
    }

    if building_bboxes.is_empty() {
        debug!("No building features found for overlap check");
        return Ok(0);
    }

    // Now check tree features from the GeoParquet file centroids stored in features dir.
    let mut overlaps = 0;
    for entry in fs::read_dir(features_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(true, |e| e != "jsonl") {
            continue;
        }
        let content = fs::read_to_string(&path)?;
        let value: Value = match serde_json::from_str(content.trim()) {
            Ok(v) => v,
            Err(_) => continue,
        };
        if let Some(objs) = value.get("CityObjects").and_then(|v| v.as_object()) {
            for (key, obj) in objs {
                if obj.get("type").and_then(|t| t.as_str()) == Some("SolitaryVegetationObject") {
                    if let Some(vertices) = value.get("vertices").and_then(|v| v.as_array()) {
                        // Compute centroid in quantized space.
                        let xs: Vec<f64> = vertices
                            .iter()
                            .filter_map(|v| v.as_array())
                            .filter_map(|a| a.first())
                            .filter_map(|v| v.as_f64())
                            .collect();
                        let ys: Vec<f64> = vertices
                            .iter()
                            .filter_map(|v| v.as_array())
                            .filter_map(|a| a.get(1))
                            .filter_map(|v| v.as_f64())
                            .collect();
                        if !xs.is_empty() && !ys.is_empty() {
                            let cx = xs.iter().sum::<f64>() / xs.len() as f64;
                            let cy = ys.iter().sum::<f64>() / ys.len() as f64;
                            for bbox in &building_bboxes {
                                if cx >= bbox[0] && cx <= bbox[2] && cy >= bbox[1] && cy <= bbox[3]
                                {
                                    warn!(
                                        "Tree '{}' centroid ({:.1}, {:.1}) overlaps a building bbox",
                                        key, cx, cy
                                    );
                                    overlaps += 1;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if overlaps > 0 {
        warn!(
            "{} tree(s) from {} have centroids overlapping building bounding boxes",
            overlaps,
            parquet_path.display()
        );
    }

    Ok(overlaps)
}

// --- Arrow column extraction helpers ---

/// Extract binary data from a Binary or LargeBinary column at row `i`.
fn extract_binary<'a>(array: &'a dyn Array, i: usize) -> Result<&'a [u8]> {
    if array.is_null(i) {
        return Ok(&[]);
    }
    if let Some(arr) = array.as_any().downcast_ref::<BinaryArray>() {
        return Ok(arr.value(i));
    }
    if let Some(arr) = array.as_any().downcast_ref::<LargeBinaryArray>() {
        return Ok(arr.value(i));
    }
    bail!(
        "Geometry column has unsupported type {:?}; expected Binary or LargeBinary",
        array.data_type()
    );
}

/// Try to extract a string from a column at row `i`.
fn extract_string(array: &dyn Array, i: usize) -> Option<String> {
    if array.is_null(i) {
        return None;
    }
    array
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| a.value(i).to_string())
}

/// Extract a generic JSON value from an Arrow column at row `i`.
fn extract_value(array: &dyn Array, i: usize) -> Option<Value> {
    if array.is_null(i) {
        return None;
    }
    match array.data_type() {
        DataType::Utf8 => array
            .as_any()
            .downcast_ref::<StringArray>()
            .map(|a| Value::String(a.value(i).to_string())),
        DataType::Int32 => {
            use arrow::array::Int32Array;
            array
                .as_any()
                .downcast_ref::<Int32Array>()
                .map(|a| json!(a.value(i)))
        }
        DataType::Int64 => {
            use arrow::array::Int64Array;
            array
                .as_any()
                .downcast_ref::<Int64Array>()
                .map(|a| json!(a.value(i)))
        }
        DataType::Float32 => {
            use arrow::array::Float32Array;
            array
                .as_any()
                .downcast_ref::<Float32Array>()
                .map(|a| json!(a.value(i)))
        }
        DataType::Float64 => {
            use arrow::array::Float64Array;
            array
                .as_any()
                .downcast_ref::<Float64Array>()
                .map(|a| json!(a.value(i)))
        }
        DataType::Boolean => {
            use arrow::array::BooleanArray;
            array
                .as_any()
                .downcast_ref::<BooleanArray>()
                .map(|a| json!(a.value(i)))
        }
        _ => None, // Unsupported types are silently skipped.
    }
}
