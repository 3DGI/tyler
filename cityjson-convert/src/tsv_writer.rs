use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use cityjson_lib::CityModel;
use csv::Terminator;

use crate::{
    tabular::{value_to_text_cell, TextCell},
    tabulate_addresses, tabulate_cityobject_hierarchy, tabulate_cityobjects,
    tabulate_model_metadata, tabulate_semantic_hierarchy, tabulate_semantic_primitives,
    AddressTable, CityObjectHierarchyTable, CityObjectTable, MetadataRow, MetadataTable,
    SemanticHierarchyTable, SemanticPrimitiveRow, SemanticPrimitiveTable, TiledMetadataTable,
    Value,
};

#[derive(Clone, Debug, Default)]
#[allow(clippy::struct_excessive_bools)]
pub struct TsvExportOptions {
    /// Optional label inserted last in every output basename.
    pub output_suffix: Option<String>,
    pub include_null_rows: bool,
    pub include_hierarchy: bool,
    pub include_cityjson_ordinal: bool,
    pub include_metadata: bool,
    pub include_semantics: bool,
    pub include_address: bool,
}

#[derive(Clone, Debug, Default)]
pub struct TsvWriteOptions {
    pub include_null_rows: bool,
    pub include_hierarchy: bool,
    pub include_cityjson_ordinal: bool,
}

/// Converts a `CityJSON` model to TSV files.
///
/// The requested output file receives the `CityObject` table. Optional tables
/// are written beside it using the output stem, for example `model_metadata.tsv`.
/// Pass an unsuffixed output path: `output_suffix` is appended last to the
/// primary and companion basenames, for example `model_metadata-lod2.2.tsv`.
///
/// # Errors
///
/// Returns an error when output files cannot be created or tabular values
/// cannot be resolved or serialized.
pub fn convert_to_tsv<P: AsRef<Path>>(
    model: &CityModel,
    output: P,
    options: &TsvExportOptions,
) -> Result<()> {
    let output = output.as_ref();
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }

    let write_options = TsvWriteOptions {
        include_null_rows: options.include_null_rows,
        include_hierarchy: options.include_hierarchy,
        include_cityjson_ordinal: options.include_cityjson_ordinal,
    };

    let cityobjects = tabulate_cityobjects(model)?;
    let primary_path = match options.output_suffix.as_deref() {
        Some(suffix) => output.with_file_name(crate::output_filename(
            output
                .file_stem()
                .and_then(|stem| stem.to_str())
                .unwrap_or("cityobjects"),
            "tsv",
            Some(suffix),
        )),
        None => output.to_path_buf(),
    };
    let mut file = File::create(primary_path)?;
    write_cityobjects_tsv(&cityobjects, &write_options, &mut file)?;

    if options.include_hierarchy {
        let hierarchy = tabulate_cityobject_hierarchy(model)?;
        let mut file = File::create(sidecar_path(
            output,
            "cityobject_hierarchy",
            options.output_suffix.as_deref(),
        ))?;
        write_cityobject_hierarchy_tsv(&hierarchy, &mut file)?;

        let hierarchy = tabulate_semantic_hierarchy(model);
        let mut file = File::create(sidecar_path(
            output,
            "semantic_hierarchy",
            options.output_suffix.as_deref(),
        ))?;
        write_semantic_hierarchy_tsv(&hierarchy, &mut file)?;
    }

    if options.include_metadata {
        let metadata = tabulate_model_metadata(model)?;
        let mut file = File::create(sidecar_path(
            output,
            "metadata",
            options.output_suffix.as_deref(),
        ))?;
        write_metadata_tsv(&metadata, &mut file)?;
    }

    if options.include_address {
        let addresses = tabulate_addresses(model)?;
        let mut file = File::create(sidecar_path(
            output,
            "addresses",
            options.output_suffix.as_deref(),
        ))?;
        write_addresses_tsv(&addresses, &write_options, &mut file)?;
    }

    if options.include_semantics {
        let semantics = tabulate_semantic_primitives(model)?;
        let mut file = File::create(sidecar_path(
            output,
            "semantics",
            options.output_suffix.as_deref(),
        ))?;
        write_semantics_tsv(&semantics, &write_options, &mut file)?;
    }

    Ok(())
}

fn sidecar_path(output: &Path, table: &str, output_suffix: Option<&str>) -> std::path::PathBuf {
    let stem = output
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("cityobjects");
    output.with_file_name(crate::output_filename(
        &format!("{stem}_{table}"),
        "tsv",
        output_suffix,
    ))
}

/// Writes `CityObject` rows as TSV.
///
/// # Errors
///
/// Returns an error when writing fails or a row value cannot be resolved.
pub fn write_cityobjects_tsv<W: Write>(
    table: &CityObjectTable<'_>,
    options: &TsvWriteOptions,
    writer: W,
) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    let mut header = vec!["cityobject_id".to_string(), "cityobject_type".to_string()];
    if options.include_cityjson_ordinal {
        header.push("cityobject_ix".to_string());
    }
    header.extend(
        table
            .schema()
            .columns
            .iter()
            .map(|column| column.name.clone()),
    );
    tsv.write_record(header)?;

    for row in table.rows() {
        let dynamic = dynamic_cells(table.model(), row.values()).with_context(|| {
            format!(
                "resolve dynamic values for CityObject {}",
                row.cityobject_id
            )
        })?;
        if !options.include_null_rows && all_null(&dynamic) {
            continue;
        }

        let mut record = vec![
            row.cityobject_id.to_string(),
            row.cityobject_type_name().to_string(),
        ];
        if options.include_cityjson_ordinal {
            record.push(row.cityobject_ix.to_string());
        }
        record.extend(dynamic.into_iter().map(|cell| cell.text));
        tsv.write_record(record)?;
    }

    tsv.flush()?;
    Ok(())
}

/// Writes `CityObject` hierarchy edges as TSV.
///
/// # Errors
///
/// Returns an error when writing fails.
pub fn write_cityobject_hierarchy_tsv<W: Write>(
    table: &CityObjectHierarchyTable<'_>,
    writer: W,
) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    tsv.write_record(["parent_id", "child_id"])?;
    for row in table.rows() {
        tsv.write_record([row.parent_id, row.child_id])?;
    }
    tsv.flush()?;
    Ok(())
}

/// Writes semantic hierarchy edges as TSV.
///
/// # Errors
///
/// Returns an error when writing fails.
pub fn write_semantic_hierarchy_tsv<W: Write>(
    table: &SemanticHierarchyTable,
    writer: W,
) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    tsv.write_record(["parent_id", "child_id"])?;
    for row in table.rows() {
        tsv.write_record([row.parent_id.to_string(), row.child_id.to_string()])?;
    }
    tsv.flush()?;
    Ok(())
}

/// Writes address rows as TSV.
///
/// # Errors
///
/// Returns an error when writing fails or address values cannot be resolved.
pub fn write_addresses_tsv<W: Write>(
    table: &AddressTable<'_>,
    _options: &TsvWriteOptions,
    writer: W,
) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    let mut header = vec!["cityobject_id".to_string(), "cityobject_type".to_string()];
    header.extend(
        table
            .schema()
            .columns
            .iter()
            .map(|column| column.name.clone()),
    );
    tsv.write_record(header)?;

    for row in table.rows() {
        let fixed = row.fixed();
        let mut record = vec![
            fixed.cityobject_id.to_string(),
            fixed.cityobject_type_name().to_string(),
        ];
        record.extend(
            dynamic_cells(table.model(), row.values())?
                .into_iter()
                .map(|cell| cell.text),
        );
        tsv.write_record(record)?;
    }

    tsv.flush()?;
    Ok(())
}

/// Writes metadata rows as TSV.
///
/// # Errors
///
/// Returns an error when writing fails or a row value cannot be resolved.
pub fn write_metadata_tsv<W: Write>(table: &MetadataTable<'_>, writer: W) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    let mut header = metadata_fixed_header();
    header.extend(
        table
            .schema()
            .columns
            .iter()
            .map(|column| column.name.clone()),
    );
    tsv.write_record(header)?;

    for row in table.rows() {
        let mut record = metadata_fixed_cells(row.fixed());
        record.extend(
            dynamic_cells(table.model(), row.values())?
                .into_iter()
                .map(|cell| cell.text),
        );
        tsv.write_record(record)?;
    }

    tsv.flush()?;
    Ok(())
}

/// Writes format-independent tiled metadata as TSV.
///
/// # Errors
///
/// Returns an error when writing fails or a metadata value cannot be resolved.
pub fn write_tiled_metadata_tsv<W: Write>(table: &TiledMetadataTable<'_>, writer: W) -> Result<()> {
    let metadata = table.metadata();
    let mut tsv = tsv_writer(writer);
    let mut header = vec![
        "tile_id".to_string(),
        "content_path".to_string(),
        "geographical_extent".to_string(),
    ];
    let mut metadata_header = metadata_fixed_header();
    metadata_header.remove(4);
    header.extend(metadata_header);
    header.extend(
        metadata
            .schema()
            .columns
            .iter()
            .map(|column| column.name.clone()),
    );
    tsv.write_record(header)?;

    let row = metadata
        .rows()
        .next()
        .context("metadata projection contains no rows")?;
    let mut fixed_cells = metadata_fixed_cells(row.fixed());
    fixed_cells.remove(4);
    let dynamic = dynamic_cells(metadata.model(), row.values())?
        .into_iter()
        .map(|cell| cell.text)
        .collect::<Vec<_>>();
    for tile in table.tiles() {
        let mut record = vec![
            tile.tile_id.clone(),
            tile.content_path.clone(),
            bbox_wkt_2d(tile.geographical_extent),
        ];
        record.extend(fixed_cells.iter().cloned());
        record.extend(dynamic.iter().cloned());
        tsv.write_record(record)?;
    }
    tsv.flush()?;
    Ok(())
}

/// Writes semantic primitive rows as TSV.
///
/// # Errors
///
/// Returns an error when writing fails or semantic attribute values cannot be
/// resolved.
pub fn write_semantics_tsv<W: Write>(
    table: &SemanticPrimitiveTable<'_>,
    _options: &TsvWriteOptions,
    writer: W,
) -> Result<()> {
    let mut tsv = tsv_writer(writer);
    let mut header = semantic_primitive_header();
    header.extend(
        table
            .schema()
            .columns
            .iter()
            .map(|column| column.name.clone()),
    );
    tsv.write_record(header)?;

    for row in table.rows() {
        let fixed = row.fixed();
        let dynamic = dynamic_cells(table.model(), row.values()).with_context(|| {
            format!(
                "resolve dynamic values for semantic primitive {} on CityObject {}",
                fixed.primitive_ix, fixed.cityobject_id
            )
        })?;

        let mut record = semantic_primitive_cells(fixed);
        record.extend(dynamic.into_iter().map(|cell| cell.text));
        tsv.write_record(record)?;
    }

    tsv.flush()?;
    Ok(())
}

fn tsv_writer<W: Write>(writer: W) -> csv::Writer<W> {
    csv::WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::Any(b'\n'))
        .from_writer(writer)
}

fn dynamic_cells<'value, 'model: 'value>(
    model: &'model CityModel,
    values: impl IntoIterator<Item = Result<Value<'value, 'model>>>,
) -> Result<Vec<TextCell>> {
    values
        .into_iter()
        .map(|value| value_to_text_cell(model, value?))
        .collect::<Result<Vec<_>>>()
}

fn all_null(cells: &[TextCell]) -> bool {
    cells.iter().all(|cell| cell.is_null)
}

fn metadata_fixed_header() -> Vec<String> {
    [
        "identifier",
        "reference_date",
        "reference_system",
        "title",
        "geographical_extent_wkt",
        "contact_name",
        "contact_email_address",
        "contact_role",
        "contact_website",
        "contact_type",
        "contact_phone",
        "contact_organization",
    ]
    .into_iter()
    .map(str::to_string)
    .collect()
}

fn metadata_fixed_cells(row: &MetadataRow<'_>) -> Vec<String> {
    vec![
        option_string_cell(row.identifier.as_deref()),
        option_string_cell(row.reference_date.as_deref()),
        option_string_cell(row.reference_system.as_deref()),
        option_string_cell(row.title.as_deref()),
        row.geographical_extent.map(bbox_wkt_2d).unwrap_or_default(),
        option_string_cell(row.contact_name.as_deref()),
        option_string_cell(row.contact_email_address.as_deref()),
        option_string_cell(row.contact_role.as_deref()),
        option_string_cell(row.contact_website.as_deref()),
        option_string_cell(row.contact_type.as_deref()),
        option_string_cell(row.contact_phone.as_deref()),
        option_string_cell(row.contact_organization.as_deref()),
    ]
}

fn bbox_wkt_2d([min_x, min_y, _, max_x, max_y, _]: [f64; 6]) -> String {
    format!(
        "POLYGON(({min_x} {min_y}, {max_x} {min_y}, {max_x} {max_y}, {min_x} {max_y}, {min_x} {min_y}))"
    )
}

fn semantic_primitive_header() -> Vec<String> {
    [
        "cityobject_id",
        "geometry_id",
        "semantic_id",
        "primitive_ix",
        "geometry_type",
        "geometry_lod",
        "semantic_type",
    ]
    .into_iter()
    .map(str::to_string)
    .collect()
}

fn semantic_primitive_cells(row: &SemanticPrimitiveRow<'_>) -> Vec<String> {
    vec![
        row.cityobject_id.to_string(),
        row.geometry_id.to_string(),
        optional_u64_cell(row.semantic_id),
        row.primitive_ix.to_string(),
        row.geometry_type.to_string(),
        option_string_cell(row.geometry_lod.as_deref()),
        row.semantic_type_name()
            .map(|semantic_type| semantic_type.to_string())
            .unwrap_or_default(),
    ]
}

fn option_string_cell(value: Option<&str>) -> String {
    value.unwrap_or_default().to_string()
}

fn optional_u64_cell(value: Option<u64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}
