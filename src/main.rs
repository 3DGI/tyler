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
#![allow(
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cloned_instead_of_copied,
    clippy::default_trait_access,
    clippy::doc_markdown,
    clippy::explicit_iter_loop,
    clippy::if_not_else,
    clippy::assigning_clones,
    clippy::manual_string_new,
    clippy::manual_assert,
    clippy::manual_is_multiple_of,
    clippy::manual_midpoint,
    clippy::match_bool,
    clippy::match_same_arms,
    clippy::needless_as_bytes,
    clippy::needless_borrows_for_generic_args,
    clippy::needless_pass_by_value,
    clippy::redundant_else,
    clippy::redundant_closure_for_method_calls,
    clippy::semicolon_if_nothing_returned,
    clippy::similar_names,
    clippy::single_match_else,
    clippy::struct_excessive_bools,
    clippy::struct_field_names,
    clippy::to_string_trait_impl,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::trivially_copy_pass_by_ref,
    clippy::type_complexity,
    clippy::unnecessary_fallible_conversions,
    clippy::unnecessary_debug_formatting,
    clippy::unnecessary_semicolon,
    clippy::unnecessary_to_owned,
    clippy::unnecessary_unwrap,
    clippy::unnecessary_wraps,
    clippy::uninlined_format_args,
    clippy::unreadable_literal,
    clippy::used_underscore_binding,
    clippy::used_underscore_items,
    clippy::stable_sort_primitive,
    clippy::useless_vec
)]
mod cli;
mod coordinates;
mod formats;
mod parser;
mod proj;
mod spatial_structs;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::coordinates::RootEnuFrame;
use crate::formats::cesium3dtiles::{Tile, TileId};
use crate::proj::Proj;
use cityjson_lib::cityjson_types::prelude::CityObjectHandle;
use clap::{Parser, ValueEnum};
use log::{debug, info, log_enabled, warn, Level};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
#[clap(rename_all = "lower")]
pub enum OutputFormatKind {
    #[value(name = "3dtiles")]
    #[default]
    Cesium3dTiles,
    Cityjson,
    Cityjsonseq,
}

#[derive(Default, Debug)]
struct DebugData {
    world: Option<PathBuf>,
    quadtree: Option<PathBuf>,
    tiles_results: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct PreparedInput {
    source: parser::InputSource,
    metadata_path: PathBuf,
    feature_base_document: Vec<u8>,
}

#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
struct TileCoord {
    level: u16,
    x: usize,
    y: usize,
}

impl TileCoord {
    fn new(level: u16, x: usize, y: usize) -> Self {
        Self { level, x, y }
    }
}

impl std::fmt::Display for TileCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}/{}", self.level, self.x, self.y)
    }
}

impl From<&spatial_structs::QuadTreeNodeId> for TileCoord {
    fn from(value: &spatial_structs::QuadTreeNodeId) -> Self {
        Self::new(value.level, value.x, value.y)
    }
}

impl From<&TileCoord> for spatial_structs::QuadTreeNodeId {
    fn from(value: &TileCoord) -> Self {
        Self::new(value.x, value.y, value.level)
    }
}

impl From<&TileId> for TileCoord {
    fn from(value: &TileId) -> Self {
        Self::new(value.level, value.x, value.y)
    }
}

impl From<&TileCoord> for TileId {
    fn from(value: &TileCoord) -> Self {
        Self::new(value.x, value.y, value.level)
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
enum TileClip {
    None,
    SourceBbox,
    GeographicRegion(GeographicBounds),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct TileExportJob {
    source_node_id: Option<TileCoord>,
    content_tile_coord: TileCoord,
    feature_ids: Vec<usize>,
    clip: TileClip,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
struct GeographicBounds {
    west: f64,
    south: f64,
    east: f64,
    north: f64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ObjectAttributeType {
    String,
    Bool,
    Int,
    Float,
}

fn build_glb_export_options(
    cli: &crate::cli::Cli,
    geometry_placement: cityjson_convert::GeometryPlacement,
    clip_bbox: Option<[f64; 6]>,
) -> cityjson_convert::ExportOptions {
    let mut feature_type_colors = BTreeMap::new();

    for (feature_type, color) in [
        ("Building", cli.color_building.as_ref()),
        ("BuildingPart", cli.color_building_part.as_ref()),
        (
            "BuildingInstallation",
            cli.color_building_installation.as_ref(),
        ),
        ("TINRelief", cli.color_tin_relief.as_ref()),
        ("Road", cli.color_road.as_ref()),
        ("Railway", cli.color_railway.as_ref()),
        ("TransportSquare", cli.color_transport_square.as_ref()),
        ("WaterBody", cli.color_water_body.as_ref()),
        ("PlantCover", cli.color_plant_cover.as_ref()),
        (
            "SolitaryVegetationObject",
            cli.color_solitary_vegetation_object.as_ref(),
        ),
        ("LandUse", cli.color_land_use.as_ref()),
        ("CityFurniture", cli.color_city_furniture.as_ref()),
        ("Bridge", cli.color_bridge.as_ref()),
        ("BridgePart", cli.color_bridge_part.as_ref()),
        ("BridgeInstallation", cli.color_bridge_installation.as_ref()),
        (
            "BridgeConstructiveElement",
            cli.color_bridge_construction_element.as_ref(),
        ),
        ("Tunnel", cli.color_tunnel.as_ref()),
        ("TunnelPart", cli.color_tunnel_part.as_ref()),
        ("TunnelInstallation", cli.color_tunnel_installation.as_ref()),
        ("GenericCityObject", cli.color_generic_city_object.as_ref()),
    ] {
        if let Some(color) = color {
            feature_type_colors.insert(feature_type.to_string(), color.clone());
        }
    }

    cityjson_convert::ExportOptions {
        native_glb_color: "#FFC0CB".to_string(),
        metadata_class_name: cli.cesium3dtiles_metadata_class.clone(),
        feature_type_colors,
        geometry_placement,
        clip_bbox,
        clip_geographic_region: None,
        smooth_normals: cli.smooth_normals,
        quantize_geometry: true,
        meshopt_compression: true,
    }
}

fn build_feature_type_lods(cli: &crate::cli::Cli) -> BTreeMap<String, String> {
    let mut feature_type_lods = BTreeMap::new();

    for (feature_type, lod) in [
        ("Building", cli.lod_building.as_ref()),
        ("BuildingPart", cli.lod_building_part.as_ref()),
        (
            "BuildingInstallation",
            cli.lod_building_installation.as_ref(),
        ),
        ("TINRelief", cli.lod_tin_relief.as_ref()),
        ("Road", cli.lod_road.as_ref()),
        ("Railway", cli.lod_railway.as_ref()),
        ("TransportSquare", cli.lod_transport_square.as_ref()),
        ("WaterBody", cli.lod_water_body.as_ref()),
        ("PlantCover", cli.lod_plant_cover.as_ref()),
        (
            "SolitaryVegetationObject",
            cli.lod_solitary_vegetation_object.as_ref(),
        ),
        ("LandUse", cli.lod_land_use.as_ref()),
        ("CityFurniture", cli.lod_city_furniture.as_ref()),
        ("Bridge", cli.lod_bridge.as_ref()),
        ("BridgePart", cli.lod_bridge_part.as_ref()),
        ("BridgeInstallation", cli.lod_bridge_installation.as_ref()),
        (
            "BridgeConstructiveElement",
            cli.lod_bridge_construction_element.as_ref(),
        ),
        ("Tunnel", cli.lod_tunnel.as_ref()),
        ("TunnelPart", cli.lod_tunnel_part.as_ref()),
        ("TunnelInstallation", cli.lod_tunnel_installation.as_ref()),
        ("GenericCityObject", cli.lod_generic_city_object.as_ref()),
    ] {
        if let Some(lod) = lod {
            feature_type_lods.insert(feature_type.to_string(), lod.clone());
        }
    }

    feature_type_lods
}

fn build_feature_filter(
    cityobject_types: Option<&Vec<parser::CityObjectType>>,
    feature_type_lods: &BTreeMap<String, String>,
) -> cityjson_index::FeatureFilter {
    cityjson_index::FeatureFilter {
        cityobject_types: cityobject_types.map(|types| {
            types
                .iter()
                .map(std::string::ToString::to_string)
                .collect::<BTreeSet<_>>()
        }),
        default_lod: cityjson_index::LodSelection::Highest,
        lods_by_type: feature_type_lods
            .iter()
            .map(|(feature_type, lod)| {
                (
                    feature_type.clone(),
                    cityjson_index::LodSelection::Exact(lod.clone()),
                )
            })
            .collect(),
    }
}

fn build_object_attribute_types(
    cli: &crate::cli::Cli,
) -> Result<BTreeMap<String, ObjectAttributeType>, Box<dyn std::error::Error>> {
    let mut attribute_types = BTreeMap::new();
    let Some(mappings) = cli.object_attributes.as_ref() else {
        return Ok(attribute_types);
    };

    for mapping in mappings {
        let (name, value_type) = mapping
            .split_once(':')
            .ok_or_else(|| format!("invalid object attribute mapping {mapping:?}"))?;
        if name.is_empty() {
            return Err(format!("object attribute name cannot be empty in {mapping:?}").into());
        }
        let value_type = match value_type {
            "string" => ObjectAttributeType::String,
            "bool" => ObjectAttributeType::Bool,
            "int" => ObjectAttributeType::Int,
            "float" => ObjectAttributeType::Float,
            _ => {
                return Err(
                    format!("invalid object attribute type {value_type:?} in {mapping:?}").into(),
                )
            }
        };
        attribute_types.insert(name.to_string(), value_type);
    }

    Ok(attribute_types)
}

fn should_dump_debug_data(cli: &crate::cli::Cli) -> bool {
    cli.debug_dump_data || log_enabled!(Level::Debug)
}

fn compute_root_enu_frame(
    world: &parser::World,
    quadtree: &spatial_structs::QuadTree,
) -> Result<RootEnuFrame, Box<dyn std::error::Error>> {
    let crs_from = format!("EPSG:{}", world.crs.to_epsg()?);
    let root_bbox = quadtree.bbox(&world.grid);
    RootEnuFrame::from_bbox(&crs_from, &root_bbox)
}

fn prepare_input(
    cli: &crate::cli::Cli,
    output_dir: &Path,
) -> Result<PreparedInput, Box<dyn std::error::Error>> {
    let resolved = cityjson_index::resolve_dataset(&cli.input, None)?;
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
    let feature_base_document = derive_base_document(&city_index)?;
    let metadata_dir = output_dir.join("metadata");
    fs::create_dir_all(&metadata_dir)?;
    let metadata_path = metadata_dir.join("cjindex-metadata.city.json");
    fs::write(&metadata_path, &feature_base_document)?;
    Ok(PreparedInput {
        source: parser::InputSource::from_cjindex_resolved(&resolved),
        metadata_path,
        feature_base_document,
    })
}

fn derive_base_document(
    city_index: &cityjson_index::CityIndex,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    // Tyler treats this as a *representative* document, not a canonical one:
    // it supplies the dataset-wide fields tyler needs (CRS, extensions, and
    // the seed for `write_cityjsonseq_auto_transform`). Per-feature
    // `transform` values are applied by cjindex during feature reads, so
    // multi-document datasets (e.g. cityjson layout) decode correctly even
    // when source metadata documents differ. CRS consistency across sources
    // is a dataset-level invariant owned by cityjson-index.
    let metadata = city_index.metadata()?;
    let Some(base_document) = metadata.first() else {
        return Err("cjindex dataset does not contain any source metadata".into());
    };
    Ok(serde_json::to_vec(base_document.as_ref())?)
}

fn collect_tile_feature_ids(
    world: &parser::World,
    qtree_node: &spatial_structs::QuadTree,
) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut feature_ids = Vec::new();
    for cellid in qtree_node.cells() {
        let cell = world.grid.cell(cellid);
        for fid in &cell.feature_ids {
            if seen.insert(*fid) {
                feature_ids.push(*fid);
            }
        }
    }
    feature_ids
}

fn explicit_tile_export_jobs(
    world: &parser::World,
    quadtree: &spatial_structs::QuadTree,
    tileset: &formats::cesium3dtiles::Tileset,
    clip_to_tile_bounds: bool,
) -> Vec<TileExportJob> {
    tileset
        .collect_leaves()
        .into_iter()
        .filter_map(|tile_ref| {
            let tile = tile_ref.clone();
            let qtree_nodeid: spatial_structs::QuadTreeNodeId = (&tile.id).into();
            let qtree_node = quadtree.node(&qtree_nodeid)?;
            let feature_ids = collect_tile_feature_ids(world, qtree_node);
            Some(TileExportJob {
                source_node_id: Some((&qtree_nodeid).into()),
                content_tile_coord: (&tile.id).into(),
                feature_ids,
                clip: if clip_to_tile_bounds {
                    TileClip::SourceBbox
                } else {
                    TileClip::None
                },
            })
        })
        .collect()
}

fn quadtree_leaf_tile_export_jobs(
    world: &parser::World,
    quadtree: &spatial_structs::QuadTree,
) -> Vec<TileExportJob> {
    quadtree
        .collect_leaves()
        .into_iter()
        .map(|qtree_node| {
            let coord = TileCoord::from(&qtree_node.id);
            TileExportJob {
                source_node_id: Some(coord.clone()),
                content_tile_coord: coord,
                feature_ids: collect_tile_feature_ids(world, qtree_node),
                clip: TileClip::None,
            }
        })
        .collect()
}

fn geographic_implicit_tile_export_jobs(
    world: &parser::World,
    quadtree: &spatial_structs::QuadTree,
    tileset: &formats::cesium3dtiles::Tileset,
    root_region: GeographicBounds,
    transformer: &Proj,
    clip_to_tile_bounds: bool,
) -> Result<Vec<TileExportJob>, Box<dyn std::error::Error>> {
    let mut feature_content_tiles: HashMap<usize, HashSet<TileId>> = HashMap::new();
    let mut feature_geographic_bounds: HashMap<usize, GeographicBounds> = HashMap::new();
    let unique_assignment = geographic_implicit_unique_assignment(world);
    let mut raw_feature_tile_assignments = 0usize;

    for tile_ref in tileset.collect_leaves() {
        let source_tile_id = tile_ref.id.clone();
        let qtree_nodeid: spatial_structs::QuadTreeNodeId = (&source_tile_id).into();
        let Some(qtree_node) = quadtree.node(&qtree_nodeid) else {
            continue;
        };
        if qtree_node.nr_items == 0 {
            continue;
        }

        for feature_id in collect_tile_feature_ids(world, qtree_node) {
            let content_level = source_tile_id.level;
            let content_tile_ids = if unique_assignment {
                vec![geographic_tile_id_for_feature_centroid(
                    root_region,
                    &world.features[feature_id],
                    content_level,
                    transformer,
                )?]
            } else {
                let feature_bounds =
                    if let Some(bounds) = feature_geographic_bounds.get(&feature_id) {
                        *bounds
                    } else {
                        let bounds = geographic_bounds_from_source_bbox(
                            &world.features[feature_id].bbox,
                            transformer,
                        )?;
                        feature_geographic_bounds.insert(feature_id, bounds);
                        bounds
                    };
                geographic_tile_ids_for_bounds(root_region, feature_bounds, content_level)
            };
            for content_tile_id in content_tile_ids {
                raw_feature_tile_assignments += 1;
                feature_content_tiles
                    .entry(feature_id)
                    .or_default()
                    .insert(content_tile_id);
            }
        }
    }

    let mut content_tile_features: HashMap<TileId, HashSet<usize>> = HashMap::new();
    let mut feature_tile_assignments = 0usize;
    for (feature_id, content_tile_ids) in feature_content_tiles {
        for content_tile_id in non_overlapping_tile_ids(content_tile_ids) {
            feature_tile_assignments += 1;
            content_tile_features
                .entry(content_tile_id)
                .or_default()
                .insert(feature_id);
        }
    }

    let mut jobs: Vec<TileExportJob> = content_tile_features
        .into_iter()
        .map(|(content_tile_id, feature_ids)| {
            let mut feature_ids: Vec<usize> = feature_ids.into_iter().collect();
            feature_ids.sort_unstable();
            TileExportJob {
                source_node_id: None,
                content_tile_coord: (&content_tile_id).into(),
                feature_ids,
                clip: if clip_to_tile_bounds {
                    TileClip::GeographicRegion(geographic_bounds_for_tile(
                        root_region,
                        &content_tile_id,
                    ))
                } else {
                    TileClip::None
                },
            }
        })
        .collect();
    jobs.sort_by(|lhs, rhs| {
        TileId::from(&lhs.content_tile_coord).cmp(&TileId::from(&rhs.content_tile_coord))
    });
    info!(
        "Geographic implicit tiling assigned {} source features to {} content tiles ({} feature-tile assignments, {} before ancestor deduplication)",
        world.features.len(),
        jobs.len(),
        feature_tile_assignments,
        raw_feature_tile_assignments
    );
    Ok(jobs)
}

fn non_overlapping_tile_ids(tile_ids: HashSet<TileId>) -> Vec<TileId> {
    let mut tile_ids: Vec<TileId> = tile_ids.into_iter().collect();
    tile_ids.sort_unstable();

    let mut retained = Vec::with_capacity(tile_ids.len());
    for candidate in &tile_ids {
        if tile_ids.iter().any(|ancestor| {
            ancestor.level < candidate.level && tile_id_is_ancestor_of(ancestor, candidate)
        }) {
            continue;
        }
        retained.push(candidate.clone());
    }

    retained
}

fn tile_id_is_ancestor_of(ancestor: &TileId, descendant: &TileId) -> bool {
    if ancestor.level >= descendant.level {
        return false;
    }
    let shift = descendant.level - ancestor.level;
    ancestor.x == (descendant.x >> shift) && ancestor.y == (descendant.y >> shift)
}

fn geographic_implicit_unique_assignment(world: &parser::World) -> bool {
    world.cityobject_types.as_ref().is_some_and(|types| {
        types.iter().any(|object_type| {
            matches!(
                object_type,
                parser::CityObjectType::Building | parser::CityObjectType::BuildingPart
            )
        })
    })
}

fn geographic_tile_id_for_feature_centroid(
    root: GeographicBounds,
    feature: &parser::Feature,
    level: u16,
    transformer: &Proj,
) -> Result<TileId, Box<dyn std::error::Error>> {
    let z = f64::midpoint(feature.bbox[2], feature.bbox[5]);
    let (lon, lat, _height) =
        transformer.convert((feature.centroid()[0], feature.centroid()[1], z))?;
    let tiles_per_axis = 1_usize << level;
    let tile_width = (root.east - root.west) / tiles_per_axis as f64;
    let tile_height = (root.north - root.south) / tiles_per_axis as f64;

    Ok(TileId::new(
        geographic_tile_index(lon, root.west, tile_width, tiles_per_axis),
        geographic_tile_index(lat, root.south, tile_height, tiles_per_axis),
        level,
    ))
}

fn geographic_bounds_from_source_bbox(
    bbox: &spatial_structs::Bbox,
    transformer: &Proj,
) -> Result<GeographicBounds, Box<dyn std::error::Error>> {
    let mut west = f64::INFINITY;
    let mut south = f64::INFINITY;
    let mut east = f64::NEG_INFINITY;
    let mut north = f64::NEG_INFINITY;

    for [x, y, z] in bbox_corners(bbox) {
        let (lon, lat, _height) = transformer.convert((x, y, z))?;
        west = west.min(lon);
        south = south.min(lat);
        east = east.max(lon);
        north = north.max(lat);
    }

    Ok(GeographicBounds {
        west,
        south,
        east,
        north,
    })
}

fn geographic_tile_ids_for_bounds(
    root: GeographicBounds,
    bounds: GeographicBounds,
    level: u16,
) -> Vec<TileId> {
    let tiles_per_axis = 1_usize << level;
    let tile_width = (root.east - root.west) / tiles_per_axis as f64;
    let tile_height = (root.north - root.south) / tiles_per_axis as f64;

    if tile_width <= 0.0 || tile_height <= 0.0 {
        return Vec::new();
    }

    let x_min = geographic_tile_index(bounds.west, root.west, tile_width, tiles_per_axis);
    let x_max = geographic_tile_index(bounds.east, root.west, tile_width, tiles_per_axis);
    let y_min = geographic_tile_index(bounds.south, root.south, tile_height, tiles_per_axis);
    let y_max = geographic_tile_index(bounds.north, root.south, tile_height, tiles_per_axis);

    let mut tile_ids = Vec::new();
    for y in y_min..=y_max {
        for x in x_min..=x_max {
            tile_ids.push(TileId::new(x, y, level));
        }
    }
    tile_ids
}

fn geographic_bounds_for_tile(root: GeographicBounds, tile_id: &TileId) -> GeographicBounds {
    let tiles_per_axis = 1_usize << tile_id.level;
    let tile_width = (root.east - root.west) / tiles_per_axis as f64;
    let tile_height = (root.north - root.south) / tiles_per_axis as f64;
    let west = root.west + tile_width * tile_id.x as f64;
    let south = root.south + tile_height * tile_id.y as f64;
    GeographicBounds {
        west,
        south,
        east: west + tile_width,
        north: south + tile_height,
    }
}

fn geographic_tile_index(value: f64, origin: f64, tile_size: f64, tiles_per_axis: usize) -> usize {
    let max_index = tiles_per_axis.saturating_sub(1) as isize;
    (((value - origin) / tile_size).floor() as isize).clamp(0, max_index) as usize
}

fn bbox_corners(bbox: &spatial_structs::Bbox) -> [[f64; 3]; 8] {
    [
        [bbox[0], bbox[1], bbox[2]],
        [bbox[0], bbox[1], bbox[5]],
        [bbox[0], bbox[4], bbox[2]],
        [bbox[0], bbox[4], bbox[5]],
        [bbox[3], bbox[1], bbox[2]],
        [bbox[3], bbox[1], bbox[5]],
        [bbox[3], bbox[4], bbox[2]],
        [bbox[3], bbox[4], bbox[5]],
    ]
}

fn read_tile_feature_models(
    world: &parser::World,
    feature_ids: &[usize],
) -> Result<Vec<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let started = Instant::now();
    let mut models = Vec::with_capacity(feature_ids.len());
    let mut cjindex_refs = Vec::with_capacity(feature_ids.len());
    for fid in feature_ids {
        match &world.features[*fid].reference {
            parser::FeatureReference::CjIndexRef(feature) => {
                cjindex_refs.push(feature.clone());
            }
            parser::FeatureReference::CjIndexId(_) => {
                let city_index = world.input_source.open_index()?;
                for fid in feature_ids {
                    let parser::FeatureReference::CjIndexId(feature_id) =
                        &world.features[*fid].reference
                    else {
                        return Err("cjindex input mixed row references with feature ids".into());
                    };
                    let model = city_index.get(feature_id)?.ok_or_else(|| {
                        format!("feature {feature_id} could not be resolved from cjindex")
                    })?;
                    models.push(model);
                }
                return Ok(models);
            }
        }
    }
    models = parser::World::read_cjindex_features_thread_local(&world.input_source, &cjindex_refs)?;
    debug!(
        "Read {} tile features from cjindex in {:?}",
        models.len(),
        started.elapsed()
    );

    Ok(models)
}

fn deduplicate_tile_feature_ids(world: &parser::World, feature_ids: &[usize]) -> Vec<usize> {
    deduplicate_feature_ids_by_reference(&world.features, feature_ids)
}

fn deduplicate_feature_ids_by_reference(
    features: &[parser::Feature],
    feature_ids: &[usize],
) -> Vec<usize> {
    let mut retained_by_source_id = BTreeMap::<String, usize>::new();

    for feature_id in feature_ids {
        let key = feature_reference_public_id(&features[*feature_id].reference);
        match retained_by_source_id.entry(key) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(*feature_id);
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                let retained = *entry.get();
                if feature_reference_precedes(
                    &features[*feature_id].reference,
                    &features[retained].reference,
                ) {
                    entry.insert(*feature_id);
                }
            }
        }
    }

    let mut retained = retained_by_source_id.into_values().collect::<Vec<_>>();
    retained.sort_unstable();
    if retained.len() != feature_ids.len() {
        debug!(
            "Deduplicated tile feature list from {} to {} source feature IDs",
            feature_ids.len(),
            retained.len()
        );
    }
    retained
}

fn feature_reference_public_id(reference: &parser::FeatureReference) -> String {
    match reference {
        parser::FeatureReference::CjIndexRef(feature) => feature.feature_id.clone(),
        parser::FeatureReference::CjIndexId(feature_id) => feature_id.clone(),
    }
}

fn feature_reference_precedes(
    lhs: &parser::FeatureReference,
    rhs: &parser::FeatureReference,
) -> bool {
    match (lhs, rhs) {
        (parser::FeatureReference::CjIndexRef(lhs), parser::FeatureReference::CjIndexRef(rhs)) => {
            (
                lhs.source_id,
                lhs.row_id,
                lhs.offset,
                lhs.length,
                &lhs.source_path,
            ) < (
                rhs.source_id,
                rhs.row_id,
                rhs.offset,
                rhs.length,
                &rhs.source_path,
            )
        }
        (parser::FeatureReference::CjIndexRef(_), parser::FeatureReference::CjIndexId(_)) => true,
        (parser::FeatureReference::CjIndexId(_), parser::FeatureReference::CjIndexRef(_)) => false,
        (parser::FeatureReference::CjIndexId(lhs), parser::FeatureReference::CjIndexId(rhs)) => {
            lhs < rhs
        }
    }
}

fn build_tile_model_from_feature_ids(
    world: &parser::World,
    feature_ids: &[usize],
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let deduplicated_feature_ids = deduplicate_tile_feature_ids(world, feature_ids);
    let models = prepare_tile_feature_models(
        world,
        &deduplicated_feature_ids,
        object_attribute_types,
        include_parent_attributes,
        false,
    )?;
    if models.is_empty() {
        return Err("tile model preparation removed all CityObjects".into());
    }
    let merge_started = Instant::now();
    let merged = cityjson_lib::ops::merge(models)?;
    debug!(
        "Merged tile model for {} selected features ({} after deduplication) in {:?}",
        feature_ids.len(),
        deduplicated_feature_ids.len(),
        merge_started.elapsed()
    );
    Ok(merged)
}

#[cfg(test)]
fn build_tile_model(
    world: &parser::World,
    qtree_node: &spatial_structs::QuadTree,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let feature_ids = collect_tile_feature_ids(world, qtree_node);
    build_tile_model_from_feature_ids(world, &feature_ids, &BTreeMap::new(), false)
}

fn build_tile_debug_cityjsonseq(
    world: &parser::World,
    feature_ids: &[usize],
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let deduplicated_feature_ids = deduplicate_tile_feature_ids(world, feature_ids);
    let models = prepare_tile_feature_models(
        world,
        &deduplicated_feature_ids,
        object_attribute_types,
        include_parent_attributes,
        true,
    )?;
    let base_root = cityjson_lib::json::from_slice(&world.feature_base_document)?;
    let mut feature_output = Vec::new();
    cityjson_convert::write_cityjsonseq(
        &mut feature_output,
        &base_root,
        &models,
        &cityjson_convert::CityJsonSeqExportOptions::default(),
    )?;
    Ok(feature_output)
}

fn build_tile_cityjsonseq_models(
    world: &parser::World,
    feature_ids: &[usize],
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
) -> Result<Vec<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let deduplicated_feature_ids = deduplicate_tile_feature_ids(world, feature_ids);
    prepare_tile_feature_models(
        world,
        &deduplicated_feature_ids,
        object_attribute_types,
        include_parent_attributes,
        true,
    )
}

fn prepare_tile_feature_models(
    world: &parser::World,
    feature_ids: &[usize],
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
    cleanup_features: bool,
) -> Result<Vec<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let models = read_tile_feature_models(world, feature_ids)?;
    models
        .into_iter()
        .zip(feature_ids.iter().copied())
        .filter_map(|(model, feature_id)| {
            match prepare_feature_model(
                model,
                feature_id,
                world.cityobject_types.as_ref(),
                &world.feature_filter,
                object_attribute_types,
                include_parent_attributes,
                cleanup_features,
            ) {
                Ok(Some(model)) => Some(Ok(model)),
                Ok(None) => None,
                Err(error) => Some(Err(error)),
            }
        })
        .collect()
}

fn prepare_feature_model(
    model: cityjson_lib::CityModel,
    _feature_id: usize,
    _cityobject_types: Option<&Vec<parser::CityObjectType>>,
    feature_filter: &cityjson_index::FeatureFilter,
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
    cleanup_feature: bool,
) -> Result<Option<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let mut model = model;
    if include_parent_attributes {
        inherit_parent_attributes(&mut model)?;
    }
    if !object_attribute_types.is_empty() {
        apply_object_attribute_types(&mut model, object_attribute_types)?;
    }
    let filtered = feature_filter.apply(&model)?;
    model = filtered.model;
    let remove_empty_geometry =
        cleanup_feature || include_parent_attributes || !object_attribute_types.is_empty();
    let model = if remove_empty_geometry {
        remove_empty_geometry_cityobjects(&model)?
    } else {
        model
    };
    if model.cityobjects().is_empty() {
        return Ok(None);
    }
    if cleanup_feature {
        let cleaned = cleanup_and_update_extents(model)?;
        Ok(Some(cleaned))
    } else {
        Ok(Some(model))
    }
}

fn inherit_parent_attributes(
    model: &mut cityjson_lib::CityModel,
) -> Result<(), Box<dyn std::error::Error>> {
    let geometry_bearing_handles = model
        .cityobjects()
        .iter()
        .filter_map(|(handle, cityobject)| {
            cityobject
                .geometry()
                .is_some_and(|geometries| !geometries.is_empty())
                .then_some(handle)
        })
        .collect::<Vec<_>>();

    for handle in geometry_bearing_handles {
        inherit_parent_attributes_for_cityobject(model, handle)?;
    }

    Ok(())
}

fn inherit_parent_attributes_for_cityobject(
    model: &mut cityjson_lib::CityModel,
    child_handle: CityObjectHandle,
) -> Result<(), Box<dyn std::error::Error>> {
    let existing_keys = model
        .cityobjects()
        .get(child_handle)
        .ok_or_else(|| {
            format!("missing CityObject handle {child_handle} during attribute inheritance")
        })?
        .attributes()
        .map(|attributes| attributes.keys().cloned().collect::<HashSet<_>>())
        .unwrap_or_default();

    let parent_handles = model
        .cityobjects()
        .get(child_handle)
        .ok_or_else(|| format!("missing CityObject handle {child_handle} during parent lookup"))?
        .parents()
        .map(<[CityObjectHandle]>::to_vec)
        .unwrap_or_default();

    let mut inherited_keys = existing_keys;
    let mut inherited_attributes = Vec::new();
    let mut visited = HashSet::new();
    collect_parent_attributes(
        model,
        &parent_handles,
        &mut visited,
        &mut inherited_keys,
        &mut inherited_attributes,
    )?;

    if inherited_attributes.is_empty() {
        return Ok(());
    }

    let cityobject = model
        .cityobjects_mut()
        .get_mut(child_handle)
        .ok_or_else(|| {
            format!("missing CityObject handle {child_handle} during attribute update")
        })?;
    let attributes = cityobject.attributes_mut();
    for (key, value) in inherited_attributes {
        attributes.insert(key, value);
    }

    Ok(())
}

fn collect_parent_attributes(
    model: &cityjson_lib::CityModel,
    parent_handles: &[CityObjectHandle],
    visited: &mut HashSet<CityObjectHandle>,
    inherited_keys: &mut HashSet<String>,
    inherited_attributes: &mut Vec<(
        String,
        cityjson_lib::cityjson_types::v2_0::OwnedAttributeValue,
    )>,
) -> Result<(), Box<dyn std::error::Error>> {
    for parent_handle in parent_handles {
        if !visited.insert(*parent_handle) {
            continue;
        }

        let Some(parent) = model.cityobjects().get(*parent_handle) else {
            return Err(format!(
                "missing parent CityObject handle {parent_handle} during attribute inheritance"
            )
            .into());
        };

        if let Some(attributes) = parent.attributes() {
            for (key, value) in attributes.iter() {
                if inherited_keys.insert(key.clone()) {
                    inherited_attributes.push((key.clone(), value.clone()));
                }
            }
        }

        if let Some(grandparents) = parent.parents() {
            collect_parent_attributes(
                model,
                grandparents,
                visited,
                inherited_keys,
                inherited_attributes,
            )?;
        }
    }

    Ok(())
}

fn filter_cityjsonfeature_preserving_root<F>(
    model: &cityjson_lib::CityModel,
    predicate: F,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>>
where
    F: FnMut(cityjson_lib::ops::CityObjectSelectionContext<'_>) -> bool,
{
    let had_feature_root = model.id().is_some();
    let selection = cityjson_lib::ops::select_cityobjects(model, predicate)?;
    let mut filtered = if selection.is_empty() {
        let mut empty = model.clone();
        empty.clear_cityobjects();
        empty.set_id(None);
        empty
    } else {
        cityjson_lib::ops::extract(model, &selection)?
    };

    if !had_feature_root || filtered.id().is_some() || filtered.cityobjects().is_empty() {
        return Ok(filtered);
    }

    let replacement_root = parentless_cityobject_handle(&filtered).ok_or(
        "filtered CityJSONFeature kept CityObjects but has no parentless replacement root",
    )?;
    filtered.set_id(Some(replacement_root));

    Ok(filtered)
}

#[cfg(test)]
pub(crate) fn filter_cityjsonfeature_preserving_root_with_policy(
    model: &cityjson_lib::CityModel,
    cityobject_types: Option<&Vec<parser::CityObjectType>>,
    feature_type_lods: &BTreeMap<String, String>,
    default_highest_lod: bool,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let filter = cityjson_index::FeatureFilter {
        cityobject_types: cityobject_types.map(|types| {
            types
                .iter()
                .map(std::string::ToString::to_string)
                .collect::<BTreeSet<_>>()
        }),
        default_lod: if default_highest_lod {
            cityjson_index::LodSelection::Highest
        } else {
            cityjson_index::LodSelection::All
        },
        lods_by_type: feature_type_lods
            .iter()
            .map(|(feature_type, lod)| {
                (
                    feature_type.clone(),
                    cityjson_index::LodSelection::Exact(lod.clone()),
                )
            })
            .collect(),
    };
    Ok(filter.apply(model)?.model)
}

fn parentless_cityobject_handle(model: &cityjson_lib::CityModel) -> Option<CityObjectHandle> {
    model.cityobjects().iter().find_map(|(handle, cityobject)| {
        let has_surviving_parent = cityobject.parents().is_some_and(|parents| {
            parents
                .iter()
                .any(|parent| model.cityobjects().get(*parent).is_some())
        });
        (!has_surviving_parent).then_some(handle)
    })
}

#[cfg(test)]
fn filter_cityobject_types(
    model: cityjson_lib::CityModel,
    cityobject_types: Option<&Vec<parser::CityObjectType>>,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    filter_cityjsonfeature_preserving_root_with_policy(
        &model,
        cityobject_types,
        &BTreeMap::new(),
        false,
    )
}

fn apply_object_attribute_types(
    model: &mut cityjson_lib::CityModel,
    object_attribute_types: &BTreeMap<String, ObjectAttributeType>,
) -> Result<(), Box<dyn std::error::Error>> {
    let handles = model.cityobjects().ids().collect::<Vec<_>>();
    for handle in handles {
        let Some(cityobject) = model.cityobjects_mut().get_mut(handle) else {
            return Err(format!(
                "missing CityObject handle {handle} during object attribute remapping"
            )
            .into());
        };
        let remapped = object_attribute_types
            .iter()
            .filter_map(|(name, attribute_type)| {
                cityobject
                    .attributes()
                    .and_then(|attributes| attributes.get(name))
                    .and_then(|value| coerce_object_attribute_value(value, *attribute_type))
                    .map(|value| (name.clone(), value))
            })
            .collect::<Vec<_>>();
        let attributes = cityobject.attributes_mut();
        attributes.clear();
        for (name, value) in remapped {
            attributes.insert(name, value);
        }
    }

    Ok(())
}

fn coerce_object_attribute_value(
    value: &cityjson_lib::cityjson_types::v2_0::OwnedAttributeValue,
    attribute_type: ObjectAttributeType,
) -> Option<cityjson_lib::cityjson_types::v2_0::OwnedAttributeValue> {
    use cityjson_lib::cityjson_types::v2_0::OwnedAttributeValue as AttributeValue;

    match attribute_type {
        ObjectAttributeType::String => match value {
            AttributeValue::Bool(value) => Some(AttributeValue::String(value.to_string())),
            AttributeValue::Unsigned(value) => Some(AttributeValue::String(value.to_string())),
            AttributeValue::Integer(value) => Some(AttributeValue::String(value.to_string())),
            AttributeValue::Float(value) if value.is_finite() => {
                Some(AttributeValue::String(value.to_string()))
            }
            AttributeValue::String(value) => Some(AttributeValue::String(value.clone())),
            _ => None,
        },
        ObjectAttributeType::Bool => match value {
            AttributeValue::Bool(value) => Some(AttributeValue::Bool(*value)),
            AttributeValue::Unsigned(value) => Some(AttributeValue::Bool(*value != 0)),
            AttributeValue::Integer(value) => Some(AttributeValue::Bool(*value != 0)),
            AttributeValue::Float(value) if value.is_finite() => {
                Some(AttributeValue::Bool(*value != 0.0))
            }
            AttributeValue::String(value) => parse_bool_attribute(value).map(AttributeValue::Bool),
            _ => None,
        },
        ObjectAttributeType::Int => match value {
            AttributeValue::Bool(value) => Some(AttributeValue::Integer(i64::from(*value))),
            AttributeValue::Unsigned(value) => {
                i64::try_from(*value).ok().map(AttributeValue::Integer)
            }
            AttributeValue::Integer(value) => Some(AttributeValue::Integer(*value)),
            AttributeValue::Float(value) if value.is_finite() =>
            {
                #[allow(clippy::cast_possible_truncation)]
                Some(AttributeValue::Integer(*value as i64))
            }
            AttributeValue::String(value) => value.parse::<i64>().ok().map(AttributeValue::Integer),
            _ => None,
        },
        ObjectAttributeType::Float => match value {
            AttributeValue::Bool(value) => Some(AttributeValue::Float(f64::from(u8::from(*value)))),
            AttributeValue::Unsigned(value) => Some(AttributeValue::Float(*value as f64)),
            AttributeValue::Integer(value) => Some(AttributeValue::Float(*value as f64)),
            AttributeValue::Float(value) if value.is_finite() => {
                Some(AttributeValue::Float(*value))
            }
            AttributeValue::String(value) => value
                .parse::<f64>()
                .ok()
                .filter(|value| value.is_finite())
                .map(AttributeValue::Float),
            _ => None,
        },
    }
}

fn parse_bool_attribute(value: &str) -> Option<bool> {
    match value.trim().to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" => Some(true),
        "false" | "0" | "no" => Some(false),
        _ => None,
    }
}

#[cfg(test)]
fn prune_lod_geometries(
    model: &mut cityjson_lib::CityModel,
    feature_type_lods: &BTreeMap<String, String>,
) -> Result<bool, Box<dyn std::error::Error>> {
    let filtered =
        filter_cityjsonfeature_preserving_root_with_policy(model, None, feature_type_lods, true)?;
    let changed = filtered.geometry_count() != model.geometry_count()
        || filtered.cityobjects().len() != model.cityobjects().len();
    *model = filtered;
    Ok(changed)
}

fn remove_empty_geometry_cityobjects(
    model: &cityjson_lib::CityModel,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    filter_cityjsonfeature_preserving_root(model, |ctx| {
        ctx.cityobject()
            .geometry()
            .is_some_and(|geometries| !geometries.is_empty())
    })
}

fn cleanup_and_update_extents(
    model: cityjson_lib::CityModel,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let mut model = cityjson_lib::ops::cleanup(&model)?;
    let handles = model.cityobjects().ids().collect::<Vec<_>>();
    for handle in handles {
        let extent = model.calculate_cityobject_geographical_extent(handle)?;
        let cityobject = model
            .cityobjects_mut()
            .get_mut(handle)
            .ok_or_else(|| format!("missing CityObject handle {handle} during extent update"))?;
        cityobject.set_geographical_extent(extent);
    }
    if let Some(extent) = model.calculate_geographical_extent()? {
        model.metadata_mut().set_geographical_extent(extent);
    }
    Ok(model)
}

fn write_debug_tile_input(
    path_features_input_dir: &Path,
    file_name: &str,
    cityjsonseq_bytes: &[u8],
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    fs::create_dir_all(path_features_input_dir)?;
    let path_tile_ndjson = path_features_input_dir
        .join(file_name)
        .with_extension("city.jsonl");
    if let Some(parent) = path_tile_ndjson.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&path_tile_ndjson, cityjsonseq_bytes)?;
    Ok(path_tile_ndjson)
}

struct TileWriteContext<'a> {
    output_tiles_dir: &'a Path,
    export_options: cityjson_convert::ExportOptions,
    source_crs: String,
    world: &'a parser::World,
    quadtree: &'a spatial_structs::QuadTree,
    object_attribute_types: BTreeMap<String, ObjectAttributeType>,
    include_parent_attributes: bool,
}

struct Cesium3dTilesPreparedOutput {
    tileset: formats::cesium3dtiles::Tileset,
    export_jobs: Vec<TileExportJob>,
    root_enu_frame: RootEnuFrame,
    source_crs: String,
    tileset_path: PathBuf,
    subtrees_path: PathBuf,
}

struct TileFilesPreparedOutput {
    export_jobs: Vec<TileExportJob>,
    source_crs: String,
}

enum PreparedOutput {
    Cesium3dTiles(Box<Cesium3dTilesPreparedOutput>),
    Cityjson(TileFilesPreparedOutput),
    Cityjsonseq(TileFilesPreparedOutput),
}

impl PreparedOutput {
    fn source_crs(&self) -> &str {
        match self {
            PreparedOutput::Cesium3dTiles(prepared) => &prepared.source_crs,
            PreparedOutput::Cityjson(prepared) | PreparedOutput::Cityjsonseq(prepared) => {
                &prepared.source_crs
            }
        }
    }

    fn root_enu_frame(&self) -> Option<&RootEnuFrame> {
        match self {
            PreparedOutput::Cesium3dTiles(prepared) => Some(&prepared.root_enu_frame),
            PreparedOutput::Cityjson(_) | PreparedOutput::Cityjsonseq(_) => None,
        }
    }
}

trait OutputFormatBackend: Sync {
    fn prepare(
        &self,
        cli: &crate::cli::Cli,
        world: &parser::World,
        quadtree: &spatial_structs::QuadTree,
        grid_cellsize: u32,
        geometric_error_factor: f64,
        debug_data_output_path: &Path,
    ) -> Result<PreparedOutput, Box<dyn std::error::Error>>;

    fn jobs(&self, prepared: &PreparedOutput) -> Vec<TileExportJob>;

    fn write_tile(
        &self,
        job: &TileExportJob,
        model: &cityjson_lib::CityModel,
        context: &TileWriteContext<'_>,
    ) -> Result<(), Box<dyn std::error::Error>>;

    fn finalize(
        &self,
        cli: &crate::cli::Cli,
        quadtree: &spatial_structs::QuadTree,
        successful_jobs: &[TileExportJob],
        failed_jobs: &[TileExportJob],
        prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>>;

    fn write_manifest_only(
        &self,
        cli: &crate::cli::Cli,
        prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>>;
}

struct Cesium3dTilesBackend;

impl OutputFormatBackend for Cesium3dTilesBackend {
    fn prepare(
        &self,
        cli: &crate::cli::Cli,
        world: &parser::World,
        quadtree: &spatial_structs::QuadTree,
        grid_cellsize: u32,
        geometric_error_factor: f64,
        debug_data_output_path: &Path,
    ) -> Result<PreparedOutput, Box<dyn std::error::Error>> {
        let tileset_path = cli.output.join("tileset.json");
        let subtrees_path = cli.output.join("subtrees");
        let tileset_path_unpruned = cli.output.join("tileset_unpruned.json");
        let subtrees_path_unpruned = cli.output.join("subtrees_unpruned");
        info!("Generating 3D Tiles tileset");
        let root_enu_frame = compute_root_enu_frame(world, quadtree)?;
        let tileset = formats::cesium3dtiles::Tileset::from_quadtree(
            quadtree,
            world,
            geometric_error_factor,
            grid_cellsize,
            cli.grid_minz,
            cli.grid_maxz,
            cli.cesium3dtiles_content_bv_from_tile,
            cli.cesium3dtiles_content_add_bv,
            &root_enu_frame,
        );

        if cli.debug_dump_grid {
            info!(
                "Exporting the explicit tileset to TSV files to {:?}",
                debug_data_output_path
            );
            tileset.export(Some(debug_data_output_path))?;
        }

        let source_crs = format!("EPSG:{}", world.crs.to_epsg()?);
        let source_to_geographic = Proj::new_known_crs(&source_crs, "EPSG:4979", None)?;
        let root_geographic_bounds =
            geographic_bounds_from_source_bbox(&quadtree.bbox(&world.grid), &source_to_geographic)?;

        let export_jobs = if cli.cesium3dtiles_implicit {
            let export_jobs = geographic_implicit_tile_export_jobs(
                world,
                quadtree,
                &tileset,
                root_geographic_bounds,
                &source_to_geographic,
                cli.cesium3dtiles_content_clip_to_tile_bounds,
            )?;
            let content_tile_ids = content_tile_ids_from_jobs(&export_jobs);
            let mut tileset_implicit = tileset.clone();
            info!("Converting to geographic implicit tiling");
            let subtrees = make_implicit_subtrees(
                &mut tileset_implicit,
                &content_tile_ids,
                &subtrees_path_unpruned,
            );

            if cli.debug_cesium3dtiles_tileset_only || should_dump_debug_data(cli) {
                info!("Writing unpruned 3D Tiles tileset");
                tileset_implicit.to_file(&tileset_path_unpruned)?;
                info!("Writing unpruned subtrees for implicit tiling");
                write_subtrees(&subtrees_path_unpruned, &subtrees)?;
            }

            export_jobs
        } else {
            let export_jobs = explicit_tile_export_jobs(
                world,
                quadtree,
                &tileset,
                cli.cesium3dtiles_content_clip_to_tile_bounds,
            );
            info!("Writing unpruned 3D Tiles tileset");
            tileset.to_file(&tileset_path_unpruned)?;
            export_jobs
        };

        Ok(PreparedOutput::Cesium3dTiles(Box::new(
            Cesium3dTilesPreparedOutput {
                tileset,
                export_jobs,
                root_enu_frame,
                source_crs,
                tileset_path,
                subtrees_path,
            },
        )))
    }

    fn jobs(&self, prepared: &PreparedOutput) -> Vec<TileExportJob> {
        let PreparedOutput::Cesium3dTiles(prepared) = prepared else {
            return Vec::new();
        };
        prepared.export_jobs.clone()
    }

    fn write_tile(
        &self,
        job: &TileExportJob,
        model: &cityjson_lib::CityModel,
        context: &TileWriteContext<'_>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let content_tile_id = TileId::from(&job.content_tile_coord);
        let file_name = content_tile_id.to_string();
        let output_file = context
            .output_tiles_dir
            .join(&file_name)
            .with_extension("glb");
        let mut export_options = context.export_options.clone();
        match job.clip {
            TileClip::None => {}
            TileClip::SourceBbox => {
                let source_node_id = job
                    .source_node_id
                    .as_ref()
                    .ok_or("source bbox clipping requires a source node ID")?;
                let qtree_nodeid: spatial_structs::QuadTreeNodeId = source_node_id.into();
                let qtree_node = context
                    .quadtree
                    .node(&qtree_nodeid)
                    .ok_or_else(|| format!("did not find tile {} in quadtree", source_node_id))?;
                export_options.clip_bbox = Some(qtree_node.bbox(&context.world.grid));
            }
            TileClip::GeographicRegion(bounds) => {
                export_options.clip_geographic_region =
                    Some(cityjson_convert::GeographicClipRegion {
                        source_crs: context.source_crs.clone(),
                        west: bounds.west,
                        south: bounds.south,
                        east: bounds.east,
                        north: bounds.north,
                    });
            }
        }

        let convert_started = Instant::now();
        cityjson_convert::convert_to_glb(model, &output_file, &export_options)?;
        debug!(
            "Converted tile {} to GLB in {:?}",
            content_tile_id,
            convert_started.elapsed()
        );
        if !output_file.exists() {
            return Err(format!("{} was not created", output_file.display()).into());
        }
        match output_file.metadata() {
            Ok(metadata) if metadata.len() == 0 => {
                debug!(
                    "Tile {} conversion produced empty GLB at {}",
                    content_tile_id,
                    output_file.display()
                );
                if let Err(error) = fs::remove_file(&output_file) {
                    warn!(
                        "Failed to remove empty GLB for tile {} at {}: {}",
                        content_tile_id,
                        output_file.display(),
                        error
                    );
                }
                Err("conversion produced an empty GLB".into())
            }
            Ok(_) => Ok(()),
            Err(error) => {
                Err(format!("could not inspect {}: {}", output_file.display(), error).into())
            }
        }
    }

    fn finalize(
        &self,
        cli: &crate::cli::Cli,
        quadtree: &spatial_structs::QuadTree,
        successful_jobs: &[TileExportJob],
        failed_jobs: &[TileExportJob],
        prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let PreparedOutput::Cesium3dTiles(prepared) = prepared else {
            return Err("3D Tiles backend received incompatible prepared output".into());
        };
        info!("Pruning tileset of {} failed tiles", failed_jobs.len());
        for (i, failed) in failed_jobs.iter().enumerate() {
            debug!(
                "{}, removing failed from the tileset: {}",
                i, failed.content_tile_coord
            );
        }
        if cli.cesium3dtiles_implicit {
            let content_tile_ids = content_tile_ids_from_jobs(successful_jobs);
            let subtrees = make_implicit_subtrees(
                &mut prepared.tileset,
                &content_tile_ids,
                &prepared.subtrees_path,
            );
            info!("Writing subtrees for implicit tiling");
            write_subtrees(&prepared.subtrees_path, &subtrees)?;
        } else {
            let failed_tiles = failed_source_tiles(&prepared.tileset, failed_jobs);
            prepared.tileset.prune(&failed_tiles, quadtree);
            write_external_tilesets_if_needed(cli, &mut prepared.tileset)?;
        }
        info!("Writing 3D Tiles tileset");
        prepared.tileset.to_file(&prepared.tileset_path)?;
        Ok(())
    }

    fn write_manifest_only(
        &self,
        cli: &crate::cli::Cli,
        prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let PreparedOutput::Cesium3dTiles(prepared) = prepared else {
            return Err("3D Tiles backend received incompatible prepared output".into());
        };
        if cli.cesium3dtiles_implicit {
            let content_tile_ids = content_tile_ids_from_jobs(&prepared.export_jobs);
            let subtrees = make_implicit_subtrees(
                &mut prepared.tileset,
                &content_tile_ids,
                &prepared.subtrees_path,
            );
            info!("Writing subtrees for implicit tiling");
            write_subtrees(&prepared.subtrees_path, &subtrees)?;
        }
        info!("Writing 3D Tiles tileset");
        prepared.tileset.to_file(&prepared.tileset_path)?;
        Ok(())
    }
}

struct CityjsonBackend;

impl OutputFormatBackend for CityjsonBackend {
    fn prepare(
        &self,
        _cli: &crate::cli::Cli,
        world: &parser::World,
        quadtree: &spatial_structs::QuadTree,
        _grid_cellsize: u32,
        _geometric_error_factor: f64,
        _debug_data_output_path: &Path,
    ) -> Result<PreparedOutput, Box<dyn std::error::Error>> {
        Ok(PreparedOutput::Cityjson(TileFilesPreparedOutput {
            export_jobs: quadtree_leaf_tile_export_jobs(world, quadtree),
            source_crs: format!("EPSG:{}", world.crs.to_epsg()?),
        }))
    }

    fn jobs(&self, prepared: &PreparedOutput) -> Vec<TileExportJob> {
        let PreparedOutput::Cityjson(prepared) = prepared else {
            return Vec::new();
        };
        prepared.export_jobs.clone()
    }

    fn write_tile(
        &self,
        job: &TileExportJob,
        model: &cityjson_lib::CityModel,
        context: &TileWriteContext<'_>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output_file = context
            .output_tiles_dir
            .join(job.content_tile_coord.to_string())
            .with_extension("city.json");
        cityjson_convert::convert_to_cityjson(
            model,
            &output_file,
            &cityjson_convert::JsonExportOptions::default(),
        )?;
        Ok(())
    }

    fn finalize(
        &self,
        _cli: &crate::cli::Cli,
        _quadtree: &spatial_structs::QuadTree,
        _successful_jobs: &[TileExportJob],
        failed_jobs: &[TileExportJob],
        _prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        info!("Skipped {} failed CityJSON tile outputs", failed_jobs.len());
        Ok(())
    }

    fn write_manifest_only(
        &self,
        _cli: &crate::cli::Cli,
        _prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }
}

struct CityjsonseqBackend;

impl OutputFormatBackend for CityjsonseqBackend {
    fn prepare(
        &self,
        _cli: &crate::cli::Cli,
        world: &parser::World,
        quadtree: &spatial_structs::QuadTree,
        _grid_cellsize: u32,
        _geometric_error_factor: f64,
        _debug_data_output_path: &Path,
    ) -> Result<PreparedOutput, Box<dyn std::error::Error>> {
        Ok(PreparedOutput::Cityjsonseq(TileFilesPreparedOutput {
            export_jobs: quadtree_leaf_tile_export_jobs(world, quadtree),
            source_crs: format!("EPSG:{}", world.crs.to_epsg()?),
        }))
    }

    fn jobs(&self, prepared: &PreparedOutput) -> Vec<TileExportJob> {
        let PreparedOutput::Cityjsonseq(prepared) = prepared else {
            return Vec::new();
        };
        prepared.export_jobs.clone()
    }

    fn write_tile(
        &self,
        job: &TileExportJob,
        _model: &cityjson_lib::CityModel,
        context: &TileWriteContext<'_>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output_file = context
            .output_tiles_dir
            .join(job.content_tile_coord.to_string())
            .with_extension("city.jsonl");
        let models = build_tile_cityjsonseq_models(
            context.world,
            &job.feature_ids,
            &context.object_attribute_types,
            context.include_parent_attributes,
        )?;
        if models.is_empty() {
            return Err("tile CityJSONSeq preparation removed all CityObjects".into());
        }
        let base_root = cityjson_lib::json::from_slice(&context.world.feature_base_document)?;
        cityjson_convert::convert_to_cityjsonseq(
            &base_root,
            &models,
            &output_file,
            &cityjson_convert::CityJsonSeqExportOptions::default(),
        )?;
        Ok(())
    }

    fn finalize(
        &self,
        _cli: &crate::cli::Cli,
        _quadtree: &spatial_structs::QuadTree,
        _successful_jobs: &[TileExportJob],
        failed_jobs: &[TileExportJob],
        _prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        info!(
            "Skipped {} failed CityJSONSeq tile outputs",
            failed_jobs.len()
        );
        Ok(())
    }

    fn write_manifest_only(
        &self,
        _cli: &crate::cli::Cli,
        _prepared: &mut PreparedOutput,
    ) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }
}

fn output_format_backend(kind: OutputFormatKind) -> Box<dyn OutputFormatBackend> {
    match kind {
        OutputFormatKind::Cesium3dTiles => Box::new(Cesium3dTilesBackend),
        OutputFormatKind::Cityjson => Box::new(CityjsonBackend),
        OutputFormatKind::Cityjsonseq => Box::new(CityjsonseqBackend),
    }
}

fn content_tile_ids_from_jobs(jobs: &[TileExportJob]) -> Vec<TileId> {
    jobs.iter()
        .map(|job| TileId::from(&job.content_tile_coord))
        .collect()
}

fn make_implicit_subtrees(
    tileset: &mut formats::cesium3dtiles::Tileset,
    content_tile_ids: &[TileId],
    subtrees_path: &Path,
) -> Vec<(TileId, Vec<u8>)> {
    let components: Vec<_> = subtrees_path
        .components()
        .map(|comp| comp.as_os_str())
        .collect();
    let subtrees_dir_option = components.last().cloned().unwrap().to_str();
    tileset.make_implicit_from_content_tile_ids(content_tile_ids, subtrees_dir_option)
}

fn write_subtrees(
    subtrees_path: &Path,
    subtrees: &[(TileId, Vec<u8>)],
) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(subtrees_path)?;
    for (subtree_id, subtree_bytes) in subtrees {
        fs::create_dir_all(subtrees_path.join(format!("{}/{}", subtree_id.level, subtree_id.x)))?;
        let out_path = subtrees_path
            .join(&subtree_id.to_string())
            .with_extension("subtree");
        let mut subtree_file = File::create(&out_path)
            .unwrap_or_else(|_| panic!("could not create {:?} for writing", &out_path));
        if let Err(_e) = subtree_file.write_all(subtree_bytes) {
            warn!("Failed to write subtree {} content", subtree_id);
        }
    }
    Ok(())
}

fn failed_source_tiles(
    tileset: &formats::cesium3dtiles::Tileset,
    failed_jobs: &[TileExportJob],
) -> Vec<Tile> {
    let failed_source_ids = failed_jobs
        .iter()
        .filter_map(|failed| failed.source_node_id.as_ref().map(TileId::from))
        .collect::<HashSet<_>>();
    tileset
        .collect_leaves()
        .into_iter()
        .filter(|tile| failed_source_ids.contains(&tile.id))
        .cloned()
        .collect()
}

fn write_external_tilesets_if_needed(
    cli: &crate::cli::Cli,
    tileset: &mut formats::cesium3dtiles::Tileset,
) -> Result<(), Box<dyn std::error::Error>> {
    let available_levels = tileset.available_levels();
    if available_levels <= 5 {
        return Ok(());
    }
    let mut split_at_level = 0;
    for level in (0..available_levels).rev() {
        let subtree_depth: u32 = (available_levels - level) as u32;
        let nr_tiles_subtree = (4_usize.pow(subtree_depth) - 1) / 3;
        let ancestor_tree_depth: u32 = (available_levels - (available_levels - level)) as u32;
        let nr_tiles_ancestor = (4_usize.pow(ancestor_tree_depth) - 1) / 3;
        if nr_tiles_ancestor < nr_tiles_subtree {
            split_at_level = level;
            break;
        }
    }
    info!(
        "Splitting the explicit tileset into external tilesets at level {}",
        split_at_level
    );
    let external_tilesets = tileset.split(split_at_level);
    for (filename, child_tileset) in &external_tilesets {
        let tileset_path = cli.output.join(filename);
        child_tileset.to_file(&tileset_path)?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // --- Begin argument parsing
    let cli = crate::cli::Cli::parse();
    debug!("{:?}", &cli);
    info!("tyler version: {}", clap::crate_version!());
    if !cli.output.is_dir() {
        fs::create_dir_all(&cli.output)?;
        info!("Created output directory {:#?}", &cli.output);
    }
    // Since we have a default value, we can safely unwrap.
    let grid_cellsize = cli.grid_cellsize.unwrap();
    let geometric_error_factor = cli.cesium3dtiles_geometric_error_factor;
    if geometric_error_factor < 0.0 {
        return Err("--3dtiles-geometric-error-factor must be non-negative".into());
    }
    // Since we have a default value, it is safe to unwrap
    // let qtree_capacity = 0; // override cli.qtree_capacity
    let qtree_criteria = spatial_structs::QuadTreeCriteria::Vertices; // override --qtree-criteria
    let quadtree_capacity = match qtree_criteria {
        spatial_structs::QuadTreeCriteria::Objects => {
            spatial_structs::QuadTreeCapacity::Objects(cli.qtree_capacity.unwrap())
        }
        spatial_structs::QuadTreeCriteria::Vertices => {
            spatial_structs::QuadTreeCapacity::Vertices(cli.qtree_capacity.unwrap())
        }
    };
    if cli.cesium3dtiles_content_bv_from_tile && !cli.cesium3dtiles_content_add_bv {
        warn!(
            "cesium3dtiles_content_bv_from_tile is true, but cesium3dtiles_content_add_bv is false. The tile content bounding volumes are not going to be added, unless you set --3dtiles-content-add-bv"
        );
    }
    let debug_data = match cli.debug_load_data {
        None => DebugData::default(),
        Some(ref dir_path) => {
            if dir_path.is_dir() {
                let world_path = dir_path.join("world.bincode");
                let quadtree_path = dir_path.join("quadtree.bincode");
                let _tileset_path = dir_path.join("tileset.bincode");
                let tiles_results_path = dir_path.join("tiles_results.bincode");
                DebugData {
                    world: world_path.exists().then_some(world_path),
                    quadtree: quadtree_path.exists().then_some(quadtree_path),
                    tiles_results: tiles_results_path.exists().then_some(tiles_results_path),
                }
            } else {
                warn!(
                    "debug_load_data {dir_path:?} is not a directory, cannot load .bincode files"
                );
                DebugData::default()
            }
        }
    };
    debug!("{:?}", debug_data);
    let debug_data_output_path = cli.output.join("debug");
    if (cli.debug_dump_grid || should_dump_debug_data(&cli)) && !debug_data_output_path.exists() {
        fs::create_dir(&debug_data_output_path)?;
    }
    let object_attribute_types = build_object_attribute_types(&cli)?;
    // --- end of argument parsing

    // Populate the World with features
    // Primitive types that implement Copy are efficiently copied into the function and
    // and it is cleaner to avoid the indirection. However, heap-allocated container
    // types are best passed by reference, because it is "expensive" to Clone them
    // (they don't implement Copy). When we move a value, we explicitly transfer
    // ownership of the value (eg cli.object_type).
    let prepared_input = if debug_data.world.is_none() {
        Some(prepare_input(&cli, &cli.output)?)
    } else {
        None
    };
    let cityobject_types = cli.object_type.clone();
    let feature_type_lods = build_feature_type_lods(&cli);
    let feature_filter = build_feature_filter(cityobject_types.as_ref(), &feature_type_lods);

    let world: parser::World = match debug_data.world {
        None => {
            let prepared_input = prepared_input
                .as_ref()
                .expect("prepared input must exist when world is built from source");
            let mut world = parser::World::from_cjindex(
                prepared_input.source.clone(),
                prepared_input.metadata_path.clone(),
                prepared_input.feature_base_document.clone(),
                grid_cellsize,
                cityobject_types,
                feature_filter,
                cli.grid_minz,
                cli.grid_maxz,
            )?;
            world.index_with_grid()?; // todo input: in general, build a line index
            if let Some(grid_path) = cli.debug_load_grid.as_ref() {
                let features_path = grid_path.parent().map(|parent| parent.join("features.tsv"));
                world.grid = spatial_structs::SquareGrid::from_debug_tsv(
                    grid_path,
                    features_path.as_deref(),
                    world.crs.to_epsg()?,
                    Some(&world.grid),
                )?;
            }
            world
        }
        Some(world_path) => {
            info!("Loading world from bincode {world_path:?}");
            let world_file = File::open(world_path)?;
            bincode::deserialize_from(world_file)?
        }
    };

    info!(
        "Computed grid statistics: {}",
        world.grid.compute_statistics()
    );

    if cli.debug_dump_grid {
        info!("Exporting the grid to TSV to {:?}", &debug_data_output_path);
        world.export_grid(cli.debug_dump_grid_features, Some(&debug_data_output_path))?;
    }
    if should_dump_debug_data(&cli) {
        debug!(
            "Exporting the world instance to bincode to {:?}",
            &debug_data_output_path
        );
        world.export_bincode(Some("world"), Some(&debug_data_output_path))?;
    }

    // Build quadtree
    let quadtree: spatial_structs::QuadTree = match debug_data.quadtree {
        None => {
            info!("Building quadtree");
            spatial_structs::QuadTree::from_world(&world, quadtree_capacity)
        }
        Some(quadtree_path) => {
            info!("Loading quadtree from bincode {quadtree_path:?}");
            let quadtree_file = File::open(quadtree_path)?;
            bincode::deserialize_from(quadtree_file)?
        }
    };

    if cli.debug_dump_grid {
        info!(
            "Exporting the quadtree to TSV to {:?}",
            &debug_data_output_path
        );
        quadtree.export(&world, Some(&debug_data_output_path))?;
    }
    if should_dump_debug_data(&cli) {
        debug!(
            "Exporting the quadtree instance to bincode to {:?}",
            &debug_data_output_path
        );
        quadtree.export_bincode(Some("quadtree"), Some(&debug_data_output_path))?;
    }

    let output_format = cli.format;
    let backend = output_format_backend(output_format);
    let mut prepared_output = backend.prepare(
        &cli,
        &world,
        &quadtree,
        grid_cellsize,
        geometric_error_factor,
        &debug_data_output_path,
    )?;
    let export_jobs = backend.jobs(&prepared_output);

    // Export each tile by merging its selected CityJSONFeature stream in memory.
    let path_output_tiles = cli.output.join("t");
    let path_features_input_dir = debug_data_output_path.join("inputs");
    if output_format != OutputFormatKind::Cesium3dTiles && cli.debug_cesium3dtiles_tileset_only {
        warn!("--debug-3dtiles-tileset-only is ignored for non-3D Tiles output formats");
    }
    if output_format != OutputFormatKind::Cesium3dTiles || !cli.debug_cesium3dtiles_tileset_only {
        fs::create_dir_all(&path_output_tiles)?;
        info!("Created output directory {:#?}", &path_output_tiles);
        if should_dump_debug_data(&cli) {
            fs::create_dir_all(&path_features_input_dir)?;
            info!("Created output directory {:#?}", &path_features_input_dir);
        }

        let geometry_placement = if let Some(root_enu_frame) = prepared_output.root_enu_frame() {
            cityjson_convert::GeometryPlacement::Enu {
                source_crs: prepared_output.source_crs().to_string(),
                ecef_origin: root_enu_frame.ecef_origin,
                east: root_enu_frame.east,
                north: root_enu_frame.north,
                up: root_enu_frame.up,
            }
        } else {
            cityjson_convert::GeometryPlacement::SourceCoordinates
        };
        let export_options = build_glb_export_options(&cli, geometry_placement, None);
        let tile_write_context = TileWriteContext {
            output_tiles_dir: &path_output_tiles,
            export_options,
            source_crs: prepared_output.source_crs().to_string(),
            world: &world,
            quadtree: &quadtree,
            object_attribute_types: object_attribute_types.clone(),
            include_parent_attributes: cli.include_parent_attributes,
        };
        let object_attribute_types = object_attribute_types.clone();
        let tiles_len = export_jobs.len();
        let tiles_failed_iter = export_jobs.par_iter().map(|job| {
            if job.feature_ids.is_empty() {
                debug!(
                    "Tile is empty ({}), skipping conversion",
                    job.content_tile_coord
                );
                return None;
            }
            let file_name = job.content_tile_coord.to_string();
            let model = match build_tile_model_from_feature_ids(
                &world,
                &job.feature_ids,
                &object_attribute_types,
                cli.include_parent_attributes,
            ) {
                Ok(model) => model,
                Err(error) => {
                    warn!(
                        "Failed to build CityJSON model for tile {}: {}",
                        job.content_tile_coord, error
                    );
                    return Some(job.clone());
                }
            };
            if should_dump_debug_data(&cli) {
                let cityjsonseq_bytes = match build_tile_debug_cityjsonseq(
                    &world,
                    &job.feature_ids,
                    &object_attribute_types,
                    cli.include_parent_attributes,
                ) {
                    Ok(bytes) => bytes,
                    Err(error) => {
                        warn!(
                            "Failed to build debug CityJSONFeature stream for tile {}: {}",
                            job.content_tile_coord, error
                        );
                        return Some(job.clone());
                    }
                };
                if let Err(error) = write_debug_tile_input(
                    &path_features_input_dir,
                    file_name.as_str(),
                    &cityjsonseq_bytes,
                ) {
                    warn!(
                        "Failed to write debug CityJSONFeature stream for tile {}: {}",
                        job.content_tile_coord, error
                    );
                    return Some(job.clone());
                }
            }
            debug!(
                "Prepared merged CityJSON model for tile {} with {} CityObjects and {} vertices",
                job.content_tile_coord,
                model.cityobjects().len(),
                model.vertices().len()
            );
            if let Err(error) = backend.write_tile(job, &model, &tile_write_context) {
                warn!(
                    "Tile {} conversion failed: {}",
                    job.content_tile_coord, error
                );
                return Some(job.clone());
            }

            None
        });

        let mut tiles_results: Vec<Option<TileExportJob>> = Vec::with_capacity(tiles_len + 2);
        if let Some(tiles_results_path) = debug_data.tiles_results {
            info!("Loading tiles_results from {tiles_results_path:?}");
            let tiles_results_file = File::open(tiles_results_path)?;
            tiles_results = bincode::deserialize_from(tiles_results_file)?
        } else {
            info!("Converting and optimizing {tiles_len} tiles");
            tiles_failed_iter.collect_into_vec(&mut tiles_results);
            if should_dump_debug_data(&cli) {
                debug!(
                    "Exporting the tiles_results instance to bincode to {:?}",
                    &debug_data_output_path
                );
                let outpath = debug_data_output_path.join("tiles_results.bincode");
                let tiles_results_file = File::create(outpath)?;
                bincode::serialize_into(tiles_results_file, &tiles_results)?;
            }
        }
        let tiles_failed: Vec<TileExportJob> = tiles_results.into_iter().flatten().collect();
        let failed_content_tile_coords: HashSet<TileCoord> = tiles_failed
            .iter()
            .map(|failed| failed.content_tile_coord.clone())
            .collect();
        let tiles_successful: Vec<TileExportJob> = export_jobs
            .iter()
            .filter(|job| !failed_content_tile_coords.contains(&job.content_tile_coord))
            .cloned()
            .collect();
        info!("Done");
        backend.finalize(
            &cli,
            &quadtree,
            &tiles_successful,
            &tiles_failed,
            &mut prepared_output,
        )?;
    } else {
        backend.write_manifest_only(&cli, &mut prepared_output)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_test_dir(prefix: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("tyler-{prefix}-{unique}"));
        fs::create_dir_all(&path).expect("create test dir");
        path
    }

    fn resource_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("data")
            .join(name)
    }

    fn build_quadtree(world: &parser::World) -> spatial_structs::QuadTree {
        spatial_structs::QuadTree::from_world(world, spatial_structs::QuadTreeCapacity::Objects(1))
    }

    fn indexed_resource_world(prefix: &str) -> (parser::World, spatial_structs::QuadTree) {
        let dataset_dir = unique_test_dir(prefix);
        let metadata =
            fs::read_to_string(resource_path("3dbag_x00.city.json")).expect("read metadata");
        let feature = fs::read_to_string(resource_path("3dbag_feature_x71.city.jsonl"))
            .expect("read feature");
        let ndjson_source = dataset_dir.join("source.city.jsonl");
        fs::write(&ndjson_source, format!("{metadata}\n{feature}\n")).expect("write ndjson source");

        let resolved =
            cityjson_index::resolve_dataset(&dataset_dir, None).expect("resolve ndjson dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex ndjson dataset");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let feature_filter = build_feature_filter(None, &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex ndjson world");
        world.index_with_grid().expect("index cjindex ndjson world");
        let quadtree = build_quadtree(&world);
        (world, quadtree)
    }

    fn resource_tileset(
        world: &parser::World,
        quadtree: &spatial_structs::QuadTree,
    ) -> formats::cesium3dtiles::Tileset {
        let root_enu_frame =
            compute_root_enu_frame(world, quadtree).expect("compute root ENU frame");
        formats::cesium3dtiles::Tileset::from_quadtree(
            quadtree,
            world,
            1.0,
            200,
            None,
            None,
            false,
            true,
            &root_enu_frame,
        )
    }

    fn feature_root_id(model: &cityjson_lib::CityModel) -> Option<String> {
        model.id().and_then(|handle| {
            model
                .cityobjects()
                .get(handle)
                .map(|cityobject| cityobject.id().to_owned())
        })
    }

    fn feature_root_repair_fixture() -> cityjson_lib::CityModel {
        cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"root-building",
                "CityObjects":{
                    "root-building":{"type":"Building","children":["building-part-1"]},
                    "building-part-1":{"type":"BuildingPart","parents":["root-building"]},
                    "other-building":{"type":"Building"}
                },
                "vertices":[]
            }"#,
        )
        .expect("feature root repair fixture should parse")
    }

    fn parent_attribute_remapping_fixture_bytes(child_attributes: serde_json::Value) -> Vec<u8> {
        let fixture = serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-parent",
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [0.0, 0.0, 0.0]
            },
            "CityObjects": {
                "building-parent": {
                    "type": "Building",
                    "attributes": {
                        "parent_only": "parent",
                        "shared": "parent",
                        "levels": 7
                    },
                    "children": ["building-part"]
                },
                "building-part": {
                    "type": "BuildingPart",
                    "parents": ["building-parent"],
                    "attributes": child_attributes,
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1",
                        "boundaries": [[[0, 1, 2], [0, 2, 3]]]
                    }]
                }
            },
            "vertices": [
                [0, 0, 0],
                [4, 0, 0],
                [4, 4, 0],
                [0, 4, 0]
            ]
        });

        serde_json::to_vec(&fixture).expect("serialize attribute inheritance fixture")
    }

    fn parent_attribute_remapping_fixture(
        child_attributes: serde_json::Value,
    ) -> cityjson_lib::CityModel {
        cityjson_lib::json::from_feature_slice(&parent_attribute_remapping_fixture_bytes(
            child_attributes,
        ))
        .expect("attribute inheritance fixture should parse")
    }

    fn feature_json(model: &cityjson_lib::CityModel) -> Value {
        let mut feature_output = Vec::new();
        cityjson_lib::json::to_feature_writer(&mut feature_output, model)
            .expect("feature should serialize");
        serde_json::from_slice(&feature_output).expect("feature json should parse")
    }

    fn geometry_relevant_signature(model: &cityjson_lib::CityModel) -> (usize, usize, Vec<String>) {
        let mut cityobject_types = model
            .cityobjects()
            .iter()
            .filter(|(_, cityobject)| {
                cityobject
                    .geometry()
                    .is_some_and(|geometries| !geometries.is_empty())
            })
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        cityobject_types.sort();
        (
            model.vertices().len(),
            model.geometry_count(),
            cityobject_types,
        )
    }

    fn feature_attribute_string(feature: &Value, object_id: &str, key: &str) -> Option<String> {
        feature
            .get("CityObjects")?
            .get(object_id)?
            .get("attributes")?
            .get(key)
            .map(|value| match value {
                Value::String(value) => value.clone(),
                Value::Number(value) => value.to_string(),
                Value::Bool(value) => value.to_string(),
                _ => value.to_string(),
            })
    }

    fn feature_attribute_value<'a>(
        feature: &'a Value,
        object_id: &str,
        key: &str,
    ) -> Option<&'a Value> {
        feature
            .get("CityObjects")?
            .get(object_id)?
            .get("attributes")?
            .get(key)
    }

    fn retained_lods_by_object_id(
        model: &cityjson_lib::CityModel,
    ) -> BTreeMap<String, Vec<String>> {
        model
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| {
                let lods = cityobject
                    .geometry()
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(|geometry_handle| {
                        model
                            .get_geometry(*geometry_handle)
                            .and_then(|geometry| geometry.lod())
                            .map(std::string::ToString::to_string)
                    })
                    .collect::<Vec<_>>();
                (cityobject.id().to_string(), lods)
            })
            .collect()
    }

    fn mixed_object_type_fixture() -> cityjson_lib::CityModel {
        cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"building",
                "CityObjects":{
                    "building":{"type":"Building","geometry":[{"type":"MultiSurface","lod":"1","boundaries":[[[0,1,2]]]}]},
                    "water":{"type":"WaterBody","geometry":[{"type":"MultiSurface","lod":"1","boundaries":[[[3,4,5]]]}]},
                    "plant":{"type":"PlantCover","geometry":[{"type":"MultiSurface","lod":"1","boundaries":[[[6,7,8]]]}]}
                },
                "vertices":[[0,0,0],[1,0,0],[0,1,0],[2,0,0],[3,0,0],[2,1,0],[4,0,0],[5,0,0],[4,1,0]]
            }"#,
        )
        .expect("mixed object type fixture should parse")
    }

    fn multi_type_lod_fixture() -> cityjson_lib::CityModel {
        cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"building",
                "CityObjects":{
                    "building":{
                        "type":"Building",
                        "geometry":[
                            {"type":"MultiSurface","lod":"1","boundaries":[[[0,1,2]]]},
                            {"type":"MultiSurface","lod":"2","boundaries":[[[0,2,3]]]}
                        ]
                    },
                    "building-part":{
                        "type":"BuildingPart",
                        "geometry":[
                            {"type":"MultiSurface","lod":"1","boundaries":[[[4,5,6]]]},
                            {"type":"MultiSurface","lod":"3","boundaries":[[[4,6,7]]]}
                        ]
                    }
                },
                "vertices":[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[2,0,0],[3,0,0],[3,1,0],[2,1,0]]
            }"#,
        )
        .expect("multi type lod fixture should parse")
    }

    fn indexed_feature_ref(feature_id: &str, source_id: i64, row_id: i64) -> parser::Feature {
        parser::Feature {
            centroid: [0.0, 0.0],
            reference: parser::FeatureReference::CjIndexRef(cityjson_index::IndexedFeatureRef {
                row_id,
                feature_id: feature_id.to_string(),
                source_id,
                source_path: PathBuf::from(format!("source-{source_id}.city.json")),
                offset: row_id as u64,
                length: 1,
                vertices_offset: None,
                vertices_length: None,
                member_ranges_json: None,
                bounds: cityjson_index::FeatureBounds {
                    min_x: 0.0,
                    max_x: 1.0,
                    min_y: 0.0,
                    max_y: 1.0,
                    min_z: 0.0,
                    max_z: 1.0,
                },
            }),
            bbox: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            needs_type_filter: false,
        }
    }

    #[test]
    fn deduplicate_feature_ids_keeps_canonical_cjindex_duplicate() {
        let features = vec![
            indexed_feature_ref("duplicate", 2, 20),
            indexed_feature_ref("unique", 2, 21),
            indexed_feature_ref("duplicate", 1, 10),
        ];

        let deduplicated = deduplicate_feature_ids_by_reference(&features, &[0, 1, 2]);

        assert_eq!(deduplicated, vec![1, 2]);
    }

    fn prepare_attribute_inheritance_model(
        model: cityjson_lib::CityModel,
        cityobject_types: Vec<parser::CityObjectType>,
        include_parent_attributes: bool,
    ) -> cityjson_lib::CityModel {
        let mut model = model;
        if include_parent_attributes {
            inherit_parent_attributes(&mut model).expect("attribute inheritance should succeed");
        }
        let model = filter_cityjsonfeature_preserving_root_with_policy(
            &model,
            Some(&cityobject_types),
            &BTreeMap::new(),
            true,
        )
        .expect("type filter should succeed");
        let model =
            remove_empty_geometry_cityobjects(&model).expect("empty object removal should succeed");
        cleanup_and_update_extents(model).expect("cleanup should succeed")
    }

    #[test]
    fn prepare_feature_model_skips_empty_geometry_removal_for_gltf_noop_path() {
        let model = cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"building-part",
                "CityObjects":{
                    "building-part":{
                        "type":"BuildingPart",
                        "geometry":[{
                            "type":"MultiSurface",
                            "lod":"1",
                            "boundaries":[[[0,1,2]]]
                        }]
                    }
                },
                "vertices":[[0,0,0],[1,0,0],[0,1,0]]
            }"#,
        )
        .expect("simple model should parse");
        let original_signature = geometry_relevant_signature(&model);

        let prepared = prepare_feature_model(
            model.clone(),
            0,
            None,
            &cityjson_index::FeatureFilter::default(),
            &BTreeMap::new(),
            false,
            false,
        )
        .expect("feature preparation should succeed")
        .expect("feature should remain");

        assert_eq!(geometry_relevant_signature(&prepared), original_signature);
        assert_eq!(
            prepared.cityobjects().len(),
            model.cityobjects().len(),
            "no-op glTF preparation should not rebuild the feature by removing empty parents"
        );
    }

    #[test]
    fn prepare_feature_model_type_filter_skip_matches_filtered_geometry() {
        let model = cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"building-part",
                "CityObjects":{
                    "building-part":{
                        "type":"BuildingPart",
                        "geometry":[{
                            "type":"MultiSurface",
                            "lod":"1",
                            "boundaries":[[[0,1,2]]]
                        }]
                    }
                },
                "vertices":[[0,0,0],[1,0,0],[0,1,0]]
            }"#,
        )
        .expect("building-part fixture should parse");
        let selected_types = vec![parser::CityObjectType::BuildingPart];
        let feature_filter = build_feature_filter(Some(&selected_types), &BTreeMap::new());

        let skipped = prepare_feature_model(
            model.clone(),
            0,
            Some(&selected_types),
            &feature_filter,
            &BTreeMap::new(),
            false,
            false,
        )
        .expect("skipped filter preparation should succeed")
        .expect("skipped filter feature should remain");
        let filtered = prepare_feature_model(
            model,
            0,
            Some(&selected_types),
            &feature_filter,
            &BTreeMap::new(),
            false,
            false,
        )
        .expect("filtered preparation should succeed")
        .expect("filtered feature should remain");

        assert_eq!(
            geometry_relevant_signature(&skipped),
            geometry_relevant_signature(&filtered)
        );
    }

    /// Test plan 1: omitting `--object-type` keeps all CityObjects, selecting
    /// one or many types filters to that union, and duplicate selections remain
    /// valid while producing the same filtered result.
    #[test]
    fn object_type_filter_supports_all_single_multi_and_duplicate_selection() {
        let model = mixed_object_type_fixture();

        let unfiltered = filter_cityobject_types(model.clone(), None).expect("unfiltered types");
        let mut unfiltered_types = unfiltered
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        unfiltered_types.sort();
        assert_eq!(
            unfiltered_types,
            vec!["Building", "PlantCover", "WaterBody"]
        );

        let building_only =
            filter_cityobject_types(model.clone(), Some(&vec![parser::CityObjectType::Building]))
                .expect("building-only filter");
        let mut building_only_types = building_only
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        building_only_types.sort();
        assert_eq!(building_only_types, vec!["Building"]);

        let union = filter_cityobject_types(
            model.clone(),
            Some(&vec![
                parser::CityObjectType::Building,
                parser::CityObjectType::WaterBody,
            ]),
        )
        .expect("union filter");
        let mut union_types = union
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        union_types.sort();
        assert_eq!(union_types, vec!["Building", "WaterBody"]);

        let duplicates = filter_cityobject_types(
            model,
            Some(&vec![
                parser::CityObjectType::Building,
                parser::CityObjectType::Building,
            ]),
        )
        .expect("duplicate selection");
        let mut duplicate_types = duplicates
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        duplicate_types.sort();
        assert_eq!(duplicate_types, vec!["Building"]);
    }

    #[test]
    fn build_feature_filter_maps_cli_types_lods_and_default_highest_policy() {
        let cityobject_types = vec![
            parser::CityObjectType::Building,
            parser::CityObjectType::BuildingPart,
        ];
        let lods = BTreeMap::from([
            ("Building".to_string(), "2.0".to_string()),
            ("BuildingPart".to_string(), "1.3".to_string()),
        ]);

        let filter = build_feature_filter(Some(&cityobject_types), &lods);

        assert_eq!(
            filter.cityobject_types,
            Some(BTreeSet::from([
                "Building".to_string(),
                "BuildingPart".to_string(),
            ]))
        );
        assert_eq!(filter.default_lod, cityjson_index::LodSelection::Highest);
        assert_eq!(
            filter.lods_by_type.get("Building"),
            Some(&cityjson_index::LodSelection::Exact("2.0".to_string()))
        );
        assert_eq!(
            filter.lods_by_type.get("BuildingPart"),
            Some(&cityjson_index::LodSelection::Exact("1.3".to_string()))
        );
    }

    #[test]
    fn shared_filter_object_type_only_keeps_all_geometries_for_selected_types() {
        let model = multi_type_lod_fixture();
        let filtered = filter_cityjsonfeature_preserving_root_with_policy(
            &model,
            Some(&vec![parser::CityObjectType::BuildingPart]),
            &BTreeMap::new(),
            false,
        )
        .expect("type-only filtering should succeed");

        let retained = retained_lods_by_object_id(&filtered);
        assert_eq!(
            retained.get("building-part"),
            Some(&vec!["1".to_string(), "3".to_string()])
        );
        assert!(!retained.contains_key("building"));
    }

    #[test]
    fn shared_filter_defaults_to_highest_lod_for_every_cityobject() {
        let model = multi_type_lod_fixture();
        let filtered = filter_cityjsonfeature_preserving_root_with_policy(
            &model,
            None,
            &BTreeMap::new(),
            true,
        )
        .expect("default highest filtering should succeed");

        let retained = retained_lods_by_object_id(&filtered);
        assert_eq!(retained.get("building"), Some(&vec!["2".to_string()]));
        assert_eq!(retained.get("building-part"), Some(&vec!["3".to_string()]));
    }

    #[test]
    fn world_extent_uses_highest_lod_geometry_for_selected_types() {
        let dataset_dir = unique_test_dir("highest-lod-extent");
        let feature_base_document =
            fs::read(resource_path("3dbag_x00.city.json")).expect("read metadata");
        let metadata_path = dataset_dir.join("metadata.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");
        let feature_path = dataset_dir.join("source.city.jsonl");
        fs::write(
            &feature_path,
            serde_json::json!({
                "type": "CityJSONFeature",
                "id": "lod-extent",
                "CityObjects": {
                    "lod-extent": {
                        "type": "BuildingPart",
                        "geometry": [
                            {
                                "type": "MultiSurface",
                                "lod": "1.0",
                                "boundaries": [[[0, 1, 2]]]
                            },
                            {
                                "type": "MultiSurface",
                                "lod": "2.0",
                                "boundaries": [[[3, 4, 5]]]
                            }
                        ]
                    }
                },
                "vertices": [
                    [0, 0, 0],
                    [1, 0, 0],
                    [0, 1, 0],
                    [10, 10, 0],
                    [11, 10, 0],
                    [10, 11, 0]
                ]
            })
            .to_string(),
        )
        .expect("write highest-lod extent feature");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve highest-lod extent dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex extent dataset");
        let feature_base_document =
            derive_base_document(&city_index).expect("derive base doc for extent dataset");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let low_types = Some(vec![parser::CityObjectType::BuildingPart]);
        let low_lods = BTreeMap::from([("BuildingPart".to_string(), "1.0".to_string())]);
        let low_filter = build_feature_filter(low_types.as_ref(), &low_lods);
        let world_low = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path.clone(),
            feature_base_document.clone(),
            200,
            low_types,
            low_filter,
            None,
            None,
        )
        .expect("build low-lod cjindex world");

        let high_types = Some(vec![parser::CityObjectType::BuildingPart]);
        let high_lods = BTreeMap::from([("BuildingPart".to_string(), "2.0".to_string())]);
        let high_filter = build_feature_filter(high_types.as_ref(), &high_lods);
        let world_high = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            high_types,
            high_filter,
            None,
            None,
        )
        .expect("build high-lod cjindex world");

        assert!(world_low.grid.bbox[0] < world_high.grid.bbox[0]);
        assert!(world_low.grid.bbox[3] < world_high.grid.bbox[3]);
    }

    #[test]
    fn feature_root_hotfix_keeps_surviving_root() {
        let model = feature_root_repair_fixture();

        let filtered =
            filter_cityjsonfeature_preserving_root(&model, |ctx| ctx.id() == "root-building")
                .expect("root-preserving filter should succeed");

        assert_eq!(
            feature_root_id(&filtered),
            Some("root-building".to_string())
        );
    }

    #[test]
    fn feature_root_hotfix_reroots_to_parentless_survivor() {
        let model = feature_root_repair_fixture();

        let filtered =
            filter_cityjsonfeature_preserving_root(&model, |ctx| ctx.id() == "other-building")
                .expect("root-repairing filter should succeed");

        assert_eq!(
            feature_root_id(&filtered),
            Some("other-building".to_string())
        );

        let mut feature_output = Vec::new();
        cityjson_lib::json::to_feature_writer(&mut feature_output, &filtered)
            .expect("repaired feature should serialize");
        let feature: Value =
            serde_json::from_slice(&feature_output).expect("serialized feature should parse");

        assert_eq!(
            feature.get("id").and_then(Value::as_str),
            Some("other-building")
        );
    }

    #[test]
    fn feature_root_hotfix_allows_empty_filtered_feature() {
        let model = feature_root_repair_fixture();

        let filtered = filter_cityjsonfeature_preserving_root(&model, |_| false)
            .expect("empty feature filter should not fail");

        assert!(filtered.cityobjects().is_empty());
        assert_eq!(feature_root_id(&filtered), None);
    }

    /// Supporting regression for test plan 1: after type filtering and cleanup,
    /// the intermediary tile keeps a correct geographical extent.
    #[test]
    fn prepare_model_filters_cityobject_types_and_updates_extent() {
        let model = cityjson_lib::json::merge_feature_stream_slice(include_bytes!(
            "../cityjson-convert/tests/data/multi_feature_types.city.jsonl"
        ))
        .expect("fixture feature stream should parse");
        let filtered =
            filter_cityobject_types(model, Some(&vec![parser::CityObjectType::Building]))
                .expect("type filter should succeed");
        let filtered = cleanup_and_update_extents(filtered).expect("cleanup should succeed");

        let cityobject_types = filtered
            .cityobjects()
            .iter()
            .map(|(_, cityobject)| cityobject.type_cityobject().to_string())
            .collect::<Vec<_>>();
        assert_eq!(cityobject_types, vec!["Building"]);
        assert_eq!(
            filtered
                .metadata()
                .and_then(|metadata| metadata.geographical_extent())
                .copied(),
            filtered
                .calculate_geographical_extent()
                .expect("extent calculation should succeed")
        );
    }

    /// Test plan 3.1: `--object-attributes` subsets the intermediary tile
    /// attributes to only the requested `name:type` mappings.
    #[test]
    fn object_attributes_subset_the_incoming_cityobject_attributes() {
        let mut model = parent_attribute_remapping_fixture(serde_json::json!({
            "child_only": "child",
            "levels": 3,
            "shared": "child"
        }));
        let object_attribute_types = BTreeMap::from([
            ("child_only".to_string(), ObjectAttributeType::String),
            ("levels".to_string(), ObjectAttributeType::Int),
        ]);

        apply_object_attribute_types(&mut model, &object_attribute_types)
            .expect("attribute subsetting should succeed");
        let feature = feature_json(&model);
        let attributes = feature["CityObjects"]["building-part"]["attributes"]
            .as_object()
            .expect("attributes should exist");

        assert_eq!(attributes.len(), 2);
        assert!(attributes.contains_key("child_only"));
        assert!(attributes.contains_key("levels"));
        assert!(!attributes.contains_key("shared"));
    }

    /// Test plan 3.2: `--object-attributes` coerces values to the requested
    /// glTF metadata scalar types for string, bool, int, and float mappings.
    #[test]
    fn object_attributes_coerce_values_to_the_requested_types() {
        let mut model = cityjson_lib::json::from_feature_slice(
            br#"{
                "type":"CityJSONFeature",
                "id":"building",
                "CityObjects":{
                    "building":{
                        "type":"Building",
                        "attributes":{
                            "as_text":7,
                            "as_bool":"true",
                            "as_int":9.0,
                            "as_float":3
                        },
                        "geometry":[{"type":"MultiSurface","lod":"1","boundaries":[[[0,1,2]]]}]
                    }
                },
                "vertices":[[0,0,0],[1,0,0],[0,1,0]]
            }"#,
        )
        .expect("attribute type fixture should parse");
        let object_attribute_types = BTreeMap::from([
            ("as_text".to_string(), ObjectAttributeType::String),
            ("as_bool".to_string(), ObjectAttributeType::Bool),
            ("as_int".to_string(), ObjectAttributeType::Int),
            ("as_float".to_string(), ObjectAttributeType::Float),
        ]);

        apply_object_attribute_types(&mut model, &object_attribute_types)
            .expect("attribute coercion should succeed");
        let feature = feature_json(&model);

        assert_eq!(
            feature_attribute_value(&feature, "building", "as_text"),
            Some(&Value::String("7".to_string()))
        );
        assert_eq!(
            feature_attribute_value(&feature, "building", "as_bool"),
            Some(&Value::Bool(true))
        );
        assert!(feature_attribute_value(&feature, "building", "as_int")
            .and_then(Value::as_i64)
            .is_some_and(|value| value == 9));
        assert!(feature_attribute_value(&feature, "building", "as_float")
            .and_then(Value::as_f64)
            .is_some_and(|value| (value - 3.0).abs() < f64::EPSILON));
    }

    /// Test plan 4.1: without any explicit `--lod-*` selectors, each
    /// multi-geometry CityObject keeps its highest available LoD by default.
    #[test]
    fn prepare_model_defaults_to_the_highest_lod_per_cityobject() {
        let mut model = cityjson_lib::json::merge_feature_stream_slice(include_bytes!(
            "../cityjson-convert/tests/data/multi_lod_building_part.city.jsonl"
        ))
        .expect("fixture feature stream should parse");

        prune_lod_geometries(&mut model, &BTreeMap::new()).expect("LoD pruning should succeed");
        let model =
            remove_empty_geometry_cityobjects(&model).expect("empty object removal should succeed");
        let model = cleanup_and_update_extents(model).expect("cleanup should succeed");

        let retained_lods = model
            .cityobjects()
            .iter()
            .flat_map(|(_, cityobject)| cityobject.geometry().unwrap_or(&[]))
            .map(|geometry_handle| {
                model
                    .get_geometry(*geometry_handle)
                    .and_then(|geometry| geometry.lod())
                    .map(std::string::ToString::to_string)
            })
            .collect::<Vec<_>>();
        assert_eq!(retained_lods, vec![Some("2.2".to_string())]);
        assert_eq!(
            model.geometry_count(),
            1,
            "cleanup should remove geometries no longer referenced by CityObjects"
        );
    }

    /// Test plan 4.2: a type-specific `--lod-*` selection applies only to that
    /// CityObject type, while unconfigured types still fall back to their
    /// highest available LoD.
    #[test]
    fn prepare_model_uses_type_specific_lod_selectors() {
        let mut model = multi_type_lod_fixture();
        let lods = BTreeMap::from([("Building".to_string(), "1".to_string())]);

        prune_lod_geometries(&mut model, &lods).expect("LoD pruning should succeed");
        let retained = retained_lods_by_object_id(&model);

        assert_eq!(retained.get("building"), Some(&vec!["1".to_string()]));
        assert_eq!(
            retained.get("building-part"),
            Some(&vec!["3".to_string()]),
            "unconfigured types should still default to their highest LoD"
        );
    }

    /// Supporting regression for the `--include-parent-attributes` refactor:
    /// when enabled, parent attributes are copied onto selected children.
    #[test]
    fn prepare_model_copies_parent_attributes_when_enabled() {
        let model = parent_attribute_remapping_fixture(serde_json::json!({}));
        let cityobject_types = vec![
            parser::CityObjectType::Building,
            parser::CityObjectType::BuildingPart,
        ];

        let disabled =
            prepare_attribute_inheritance_model(model.clone(), cityobject_types.clone(), false);
        let disabled_feature = feature_json(&disabled);
        assert_eq!(
            feature_attribute_string(&disabled_feature, "building-part", "parent_only"),
            None
        );
        assert_eq!(
            feature_attribute_string(&disabled_feature, "building-part", "shared"),
            None
        );

        let enabled = prepare_attribute_inheritance_model(model, cityobject_types, true);
        let enabled_feature = feature_json(&enabled);
        assert_eq!(
            feature_attribute_string(&enabled_feature, "building-part", "parent_only"),
            Some("parent".to_string())
        );
        assert_eq!(
            feature_attribute_string(&enabled_feature, "building-part", "shared"),
            Some("parent".to_string())
        );
        assert_eq!(feature_root_id(&enabled), Some("building-part".to_string()));
    }

    /// Supporting regression for the `--include-parent-attributes` refactor:
    /// child attributes win when parent and child define the same key.
    #[test]
    fn prepare_model_keeps_child_attributes_on_conflict() {
        let model = parent_attribute_remapping_fixture(serde_json::json!({
            "child_only": "child",
            "levels": 3,
            "shared": "child"
        }));
        let cityobject_types = vec![
            parser::CityObjectType::Building,
            parser::CityObjectType::BuildingPart,
        ];

        let prepared = prepare_attribute_inheritance_model(model, cityobject_types, true);
        let prepared_feature = feature_json(&prepared);

        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "parent_only"),
            Some("parent".to_string())
        );
        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "child_only"),
            Some("child".to_string())
        );
        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "levels"),
            Some("3".to_string())
        );
        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "shared"),
            Some("child".to_string())
        );
    }

    /// Test plan 5.1: `--include-parent-attributes --object-type BuildingPart`
    /// still inherits the parent `Building` attributes even when the parent
    /// object itself is filtered out of the intermediary tile.
    #[test]
    fn prepare_model_inherits_parent_attributes_when_parent_type_is_not_selected() {
        let model = parent_attribute_remapping_fixture(serde_json::json!({}));

        let prepared = prepare_attribute_inheritance_model(
            model,
            vec![parser::CityObjectType::BuildingPart],
            true,
        );
        let prepared_feature = feature_json(&prepared);

        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "parent_only"),
            Some("parent".to_string())
        );
        assert_eq!(
            feature_attribute_string(&prepared_feature, "building-part", "shared"),
            Some("parent".to_string())
        );
        assert_eq!(
            prepared_feature
                .get("CityObjects")
                .and_then(Value::as_object)
                .map(|objects| objects.len()),
            Some(1)
        );
        assert_eq!(
            feature_root_id(&prepared),
            Some("building-part".to_string())
        );
    }

    #[test]
    fn build_tile_model_remaps_parent_attributes_before_glb_conversion() {
        let dataset_dir = unique_test_dir("attribute-inheritance");
        let metadata: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_x00.city.json")).expect("read metadata"),
        )
        .expect("parse metadata");
        let feature_bytes = parent_attribute_remapping_fixture_bytes(serde_json::json!({}));
        let feature_str =
            String::from_utf8(feature_bytes).expect("attribute fixture should be utf8");
        let metadata_str = serde_json::to_string(&metadata).expect("serialize metadata");
        let ndjson_source = dataset_dir.join("source.city.jsonl");
        fs::write(&ndjson_source, format!("{metadata_str}\n{feature_str}\n"))
            .expect("write ndjson source");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve attribute-inheritance dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let cityobject_types = Some(vec![
            parser::CityObjectType::Building,
            parser::CityObjectType::BuildingPart,
        ]);
        let feature_filter = build_feature_filter(cityobject_types.as_ref(), &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            cityobject_types,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex world");
        world.index_with_grid().expect("index cjindex world");
        let quadtree = build_quadtree(&world);
        let feature_ids = collect_tile_feature_ids(&world, &quadtree);

        let model = build_tile_model_from_feature_ids(&world, &feature_ids, &BTreeMap::new(), true)
            .expect("build tile model with inherited attributes");
        let model_feature = feature_json(&model);

        assert_eq!(
            feature_attribute_string(&model_feature, "building-part", "parent_only"),
            Some("parent".to_string())
        );
        assert_eq!(
            feature_attribute_string(&model_feature, "building-part", "shared"),
            Some("parent".to_string())
        );
        assert_eq!(feature_root_id(&model), Some("building-part".to_string()));
    }

    #[test]
    fn build_tile_model_inherits_parent_attributes_when_only_child_type_is_selected() {
        let dataset_dir = unique_test_dir("attribute-inheritance-child-only");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::write(
            &feature_path,
            parent_attribute_remapping_fixture_bytes(serde_json::json!({})),
        )
        .expect("write feature");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve child-only attribute-inheritance dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let cityobject_types = Some(vec![parser::CityObjectType::BuildingPart]);
        let feature_filter = build_feature_filter(cityobject_types.as_ref(), &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            cityobject_types,
            feature_filter,
            None,
            None,
        )
        .expect("build legacy world");
        world.index_with_grid().expect("index legacy world");
        let quadtree = build_quadtree(&world);
        let feature_ids = collect_tile_feature_ids(&world, &quadtree);

        let model = build_tile_model_from_feature_ids(&world, &feature_ids, &BTreeMap::new(), true)
            .expect("build tile model with inherited attributes");
        let model_feature = feature_json(&model);

        assert_eq!(
            feature_attribute_string(&model_feature, "building-part", "parent_only"),
            Some("parent".to_string())
        );
        assert_eq!(
            feature_attribute_string(&model_feature, "building-part", "shared"),
            Some("parent".to_string())
        );
        assert_eq!(
            model_feature
                .get("CityObjects")
                .and_then(Value::as_object)
                .map(|objects| objects.len()),
            Some(1)
        );
        assert_eq!(feature_root_id(&model), Some("building-part".to_string()));
    }

    #[test]
    fn single_building_type_filter_retains_child_geometry_for_tiles() {
        let dataset_dir = unique_test_dir("building-parent-child-geometry");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::write(
            &feature_path,
            parent_attribute_remapping_fixture_bytes(serde_json::json!({})),
        )
        .expect("write feature");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve building-only parent-child dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let cityobject_types = Some(vec![parser::CityObjectType::Building]);
        let feature_filter = build_feature_filter(cityobject_types.as_ref(), &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            cityobject_types,
            feature_filter,
            None,
            None,
        )
        .expect("build building-only world");
        world.index_with_grid().expect("index building-only world");
        let quadtree = build_quadtree(&world);
        let feature_ids = collect_tile_feature_ids(&world, &quadtree);

        let model =
            build_tile_model_from_feature_ids(&world, &feature_ids, &BTreeMap::new(), false)
                .expect("build tile model");
        let retained = retained_lods_by_object_id(&model);

        assert!(retained
            .get("building-part")
            .is_some_and(|lods| !lods.is_empty() && lods.iter().all(|lod| lod == "1")));
        assert!(model.vertices().len() >= 3);
    }

    #[test]
    fn world_from_cjindex_errors_before_export_when_explicit_lod_is_missing() {
        let dataset_dir = unique_test_dir("missing-explicit-lod");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::write(
            &feature_path,
            parent_attribute_remapping_fixture_bytes(serde_json::json!({})),
        )
        .expect("write feature");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve missing-lod dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let lods = BTreeMap::from([("BuildingPart".to_string(), "99".to_string())]);
        let feature_filter = build_feature_filter(None, &lods);
        let Err(error) = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
            feature_filter,
            None,
            None,
        ) else {
            panic!("missing explicit LoD should fail before grid/export");
        };
        let message = error.to_string();

        assert!(message.contains("requested LoD selector matched no geometry"));
        assert!(message.contains("BuildingPart requested LoD '99'"));
        assert!(message.contains("available LoDs are: 1"));
    }

    #[test]
    fn write_debug_tile_input_writes_cityjsonl() {
        let dataset_dir = unique_test_dir("debug-tile-input");
        let inputs_dir = dataset_dir.join("inputs");
        let path = write_debug_tile_input(&inputs_dir, "tile", b"{\"type\":\"CityJSONFeature\"}\n")
            .expect("write debug tile input");

        assert_eq!(path, inputs_dir.join("tile.city.jsonl"));
        assert_eq!(
            fs::read(&path).expect("read debug tile input"),
            b"{\"type\":\"CityJSONFeature\"}\n"
        );
        assert!(!inputs_dir.join("tile.input").exists());

        let nested_path =
            write_debug_tile_input(&inputs_dir, "1/2/3", b"{\"type\":\"CityJSONFeature\"}\n")
                .expect("write nested debug tile input");
        assert_eq!(nested_path, inputs_dir.join("1/2/3.city.jsonl"));
        assert!(nested_path.exists());
    }

    #[test]
    fn build_tile_model_exports_cjindex_ndjson_directly() {
        let dataset_dir = unique_test_dir("cjindex-ndjson");
        let metadata =
            fs::read_to_string(resource_path("3dbag_x00.city.json")).expect("read metadata");
        let feature = fs::read_to_string(resource_path("3dbag_feature_x71.city.jsonl"))
            .expect("read feature");
        let ndjson_source = dataset_dir.join("source.city.jsonl");
        fs::write(&ndjson_source, format!("{metadata}\n{feature}\n")).expect("write ndjson source");

        let resolved =
            cityjson_index::resolve_dataset(&dataset_dir, None).expect("resolve ndjson dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex ndjson dataset");
        let indexed_bounds = city_index
            .iter_all_bbox_pages(1)
            .expect("build bbox page iterator")
            .next()
            .expect("bbox page should exist")
            .expect("bbox page should load")
            .into_iter()
            .next()
            .expect("indexed feature should exist")
            .bounds;
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let feature_filter = build_feature_filter(None, &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex ndjson world");
        #[allow(clippy::float_cmp)]
        {
            assert_eq!(world.grid.bbox[2], indexed_bounds.min_z);
            assert_eq!(world.grid.bbox[5], indexed_bounds.max_z);
        }
        world.index_with_grid().expect("index cjindex ndjson world");
        assert!(world
            .features
            .iter()
            .all(|feature| matches!(feature.reference, parser::FeatureReference::CjIndexRef(_))));
        let quadtree = build_quadtree(&world);
        let model = build_tile_model(&world, &quadtree).expect("build tile model");

        assert!(!model.cityobjects().is_empty());
        assert!(!model.vertices().is_empty());
    }

    #[test]
    fn build_tile_model_exports_cjindex_ndjson_without_type_filter_directly() {
        let dataset_dir = unique_test_dir("cjindex-ndjson-unfiltered");
        let metadata =
            fs::read_to_string(resource_path("3dbag_x00.city.json")).expect("read metadata");
        let feature = fs::read_to_string(resource_path("3dbag_feature_x71.city.jsonl"))
            .expect("read feature");
        let ndjson_source = dataset_dir.join("source.city.jsonl");
        fs::write(&ndjson_source, format!("{metadata}\n{feature}\n")).expect("write ndjson source");

        let resolved =
            cityjson_index::resolve_dataset(&dataset_dir, None).expect("resolve ndjson dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex ndjson dataset");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let feature_filter = build_feature_filter(None, &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex ndjson world");
        world.index_with_grid().expect("index cjindex ndjson world");
        let quadtree = build_quadtree(&world);
        let model = build_tile_model(&world, &quadtree).expect("build tile model");

        assert!(!model.cityobjects().is_empty());
        assert!(!model.vertices().is_empty());
    }

    #[test]
    fn build_tile_model_exports_cjindex_cityjson_directly() {
        let dataset_dir = unique_test_dir("cjindex-cityjson");
        let metadata: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_x00.city.json")).expect("read metadata"),
        )
        .expect("parse metadata");
        let feature: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_feature_x71.city.jsonl")).expect("read feature"),
        )
        .expect("parse feature");
        let mut cityjson = metadata;
        cityjson["CityObjects"] = feature["CityObjects"].clone();
        cityjson["vertices"] = feature["vertices"].clone();
        let cityjson_path = dataset_dir.join("source.city.json");
        fs::write(
            &cityjson_path,
            serde_json::to_vec(&cityjson).expect("serialize cityjson"),
        )
        .expect("write cityjson source");

        let resolved =
            cityjson_index::resolve_dataset(&dataset_dir, None).expect("resolve cityjson dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex cityjson dataset");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let cityobject_types = Some(vec![parser::CityObjectType::Building]);
        let feature_filter = build_feature_filter(cityobject_types.as_ref(), &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            cityobject_types,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex cityjson world");
        world
            .index_with_grid()
            .expect("index cjindex cityjson world");
        let quadtree = build_quadtree(&world);
        let model = build_tile_model(&world, &quadtree).expect("build tile model");

        assert!(!model.cityobjects().is_empty());
        assert!(!model.vertices().is_empty());
    }

    #[test]
    fn build_tile_model_exports_cjindex_cityjson_multi_document() {
        let dataset_dir = unique_test_dir("cjindex-cityjson-multi");
        let metadata: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_x00.city.json")).expect("read metadata"),
        )
        .expect("parse metadata");
        let feature: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_feature_x71.city.jsonl")).expect("read feature"),
        )
        .expect("parse feature");

        let mut document_a = metadata.clone();
        document_a["CityObjects"] = feature["CityObjects"].clone();
        document_a["vertices"] = feature["vertices"].clone();
        let path_a = dataset_dir.join("source_a.city.json");
        fs::write(
            &path_a,
            serde_json::to_vec(&document_a).expect("serialize cityjson a"),
        )
        .expect("write cityjson a");

        // Second document keeps the same CRS but differs from document_a in
        // bytes (distinct title) and in transform.translate. cityjson-index
        // reconstructs each feature against its own source metadata, so the
        // shifted translate must not break decoding.
        let mut document_b = metadata.clone();
        if let Some(meta) = document_b
            .get_mut("metadata")
            .and_then(Value::as_object_mut)
        {
            meta.insert("title".to_string(), Value::String("variant-b".to_string()));
        }
        if let Some(translate) = document_b
            .get_mut("transform")
            .and_then(|t| t.get_mut("translate"))
            .and_then(Value::as_array_mut)
        {
            for component in translate.iter_mut() {
                if let Some(value) = component.as_f64() {
                    *component = serde_json::json!(value + 1.0);
                }
            }
        }
        document_b["CityObjects"] = feature["CityObjects"].clone();
        document_b["vertices"] = feature["vertices"].clone();
        let path_b = dataset_dir.join("source_b.city.json");
        fs::write(
            &path_b,
            serde_json::to_vec(&document_b).expect("serialize cityjson b"),
        )
        .expect("write cityjson b");

        let resolved = cityjson_index::resolve_dataset(&dataset_dir, None)
            .expect("resolve multi-document cityjson dataset");
        let mut city_index =
            cityjson_index::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index
            .reindex()
            .expect("reindex multi-document cityjson dataset");

        let stored_metadata = city_index.metadata().expect("load source metadata");
        assert!(
            stored_metadata.len() >= 2,
            "expected at least two source metadata documents, got {}",
            stored_metadata.len()
        );

        let feature_base_document =
            derive_base_document(&city_index).expect("derive base doc from multi-document dataset");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let cityobject_types = Some(vec![parser::CityObjectType::Building]);
        let feature_filter = build_feature_filter(cityobject_types.as_ref(), &BTreeMap::new());
        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            cityobject_types,
            feature_filter,
            None,
            None,
        )
        .expect("build cjindex cityjson world");
        world
            .index_with_grid()
            .expect("index cjindex cityjson world");

        // Both source documents must contribute features. Per-source transforms
        // are applied by cjindex during read, not by tyler.
        assert!(
            world.features.len() >= 2,
            "expected features from both source documents, got {}",
            world.features.len()
        );

        let quadtree = build_quadtree(&world);
        let model = build_tile_model(&world, &quadtree).expect("build tile model");
        assert!(!model.cityobjects().is_empty());
        assert!(!model.vertices().is_empty());
    }

    #[test]
    fn tile_coord_converts_to_quadtree_node_and_3d_tile_id() {
        let coord = TileCoord::new(3, 4, 5);

        let qtree_node_id = spatial_structs::QuadTreeNodeId::from(&coord);
        assert_eq!(qtree_node_id.level, 3);
        assert_eq!(qtree_node_id.x, 4);
        assert_eq!(qtree_node_id.y, 5);
        assert_eq!(TileCoord::from(&qtree_node_id), coord);

        let tile_id = TileId::from(&coord);
        assert_eq!(tile_id, TileId::new(4, 5, 3));
        assert_eq!(TileCoord::from(&tile_id), coord);
    }

    #[test]
    fn explicit_quadtree_leaf_job_generation_uses_leaf_nodes() {
        let (world, quadtree) = indexed_resource_world("explicit-leaf-jobs");
        let tileset = resource_tileset(&world, &quadtree);

        let jobs = explicit_tile_export_jobs(&world, &quadtree, &tileset, true);
        let leaves = tileset.collect_leaves();

        assert_eq!(jobs.len(), leaves.len());
        assert!(jobs
            .iter()
            .all(|job| matches!(job.clip, TileClip::SourceBbox)));
        for job in jobs {
            let source_node_id = job
                .source_node_id
                .as_ref()
                .expect("explicit jobs should keep source node IDs");
            let qtree_node_id = spatial_structs::QuadTreeNodeId::from(source_node_id);
            let qtree_node = quadtree
                .node(&qtree_node_id)
                .expect("job source node should exist in quadtree");
            assert_eq!(
                job.feature_ids,
                collect_tile_feature_ids(&world, qtree_node)
            );
            assert_eq!(job.content_tile_coord, *source_node_id);
        }
    }

    #[test]
    fn failed_job_pruning_input_uses_explicit_source_tiles() {
        let (world, quadtree) = indexed_resource_world("failed-pruning-input");
        let tileset = resource_tileset(&world, &quadtree);
        let jobs = explicit_tile_export_jobs(&world, &quadtree, &tileset, false);
        let failed_job = jobs
            .iter()
            .find(|job| job.source_node_id.is_some())
            .expect("expected an explicit tile job")
            .clone();
        let expected_tile_id = TileId::from(failed_job.source_node_id.as_ref().unwrap());

        let failed_tiles = failed_source_tiles(&tileset, &[failed_job]);

        assert_eq!(failed_tiles.len(), 1);
        assert_eq!(failed_tiles[0].id, expected_tile_id);
    }

    #[test]
    fn geographic_implicit_jobs_are_sorted_and_non_overlapping() {
        let (world, quadtree) = indexed_resource_world("geographic-implicit-jobs");
        let tileset = resource_tileset(&world, &quadtree);
        let source_crs = format!("EPSG:{}", world.crs.to_epsg().expect("EPSG code"));
        let source_to_geographic =
            Proj::new_known_crs(&source_crs, "EPSG:4979", None).expect("source CRS transform");
        let root_geographic_bounds =
            geographic_bounds_from_source_bbox(&quadtree.bbox(&world.grid), &source_to_geographic)
                .expect("root geographic bounds");

        let jobs = geographic_implicit_tile_export_jobs(
            &world,
            &quadtree,
            &tileset,
            root_geographic_bounds,
            &source_to_geographic,
            true,
        )
        .expect("geographic implicit jobs");

        assert!(!jobs.is_empty());
        let content_tile_ids = content_tile_ids_from_jobs(&jobs);
        assert!(content_tile_ids.windows(2).all(|pair| pair[0] <= pair[1]));
        for (index, candidate) in content_tile_ids.iter().enumerate() {
            assert!(!content_tile_ids
                .iter()
                .enumerate()
                .any(|(ancestor_index, ancestor)| {
                    ancestor_index != index && tile_id_is_ancestor_of(ancestor, candidate)
                }));
        }
        assert!(jobs
            .iter()
            .all(|job| matches!(job.clip, TileClip::GeographicRegion(_))));
    }

    #[test]
    fn geographic_bounds_map_to_lon_lat_implicit_tile_ids() {
        let root = GeographicBounds {
            west: 0.0,
            south: 50.0,
            east: 4.0,
            north: 54.0,
        };
        let bounds = GeographicBounds {
            west: 2.1,
            south: 51.1,
            east: 3.9,
            north: 52.9,
        };

        let tile_ids = geographic_tile_ids_for_bounds(root, bounds, 2);

        assert_eq!(
            tile_ids,
            vec![
                TileId::new(2, 1, 2),
                TileId::new(3, 1, 2),
                TileId::new(2, 2, 2),
                TileId::new(3, 2, 2),
            ]
        );
    }

    #[test]
    fn geographic_bounds_for_tile_matches_implicit_subdivision() {
        let root = GeographicBounds {
            west: 0.0,
            south: 50.0,
            east: 4.0,
            north: 54.0,
        };

        let bounds = geographic_bounds_for_tile(root, &TileId::new(2, 1, 2));

        assert!((bounds.west - 2.0).abs() < f64::EPSILON);
        assert!((bounds.south - 51.0).abs() < f64::EPSILON);
        assert!((bounds.east - 3.0).abs() < f64::EPSILON);
        assert!((bounds.north - 52.0).abs() < f64::EPSILON);
    }

    #[test]
    fn non_overlapping_tile_ids_removes_descendants_when_ancestor_has_content() {
        let tile_ids = HashSet::from([
            TileId::new(0, 0, 1),
            TileId::new(0, 0, 2),
            TileId::new(1, 0, 2),
            TileId::new(2, 0, 2),
            TileId::new(3, 3, 2),
        ]);

        let retained = non_overlapping_tile_ids(tile_ids);

        assert_eq!(
            retained,
            vec![
                TileId::new(0, 0, 1),
                TileId::new(2, 0, 2),
                TileId::new(3, 3, 2),
            ]
        );
    }
}
