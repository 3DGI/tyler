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

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::coordinates::RootEnuFrame;
use crate::formats::cesium3dtiles::{Tile, TileId};
use crate::proj::Proj;
use cityjson_lib::cityjson::prelude::{CityObjectHandle, GeometryHandle};
use clap::Parser;
use log::{debug, info, log_enabled, warn, Level};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, clap::ValueEnum, Eq, PartialEq)]
#[clap(rename_all = "lower")]
pub enum Formats {
    _3DTiles,
    CityJSON,
}

impl ToString for Formats {
    fn to_string(&self) -> String {
        match self {
            Formats::_3DTiles => "3DTiles".to_string(),
            Formats::CityJSON => "CityJSON".to_string(),
        }
    }
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
    feature_base_document: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct TileExportJob {
    source_tile: Option<Tile>,
    source_tile_id: Option<TileId>,
    content_tile_id: TileId,
    feature_ids: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct GeographicBounds {
    west: f64,
    south: f64,
    east: f64,
    north: f64,
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
        metadata_class_name: cli
            .cesium3dtiles_metadata_class
            .clone()
            .unwrap_or_else(|| "cityobject".to_string()),
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
    match cityjson_index::resolve_dataset(&cli.input, None) {
        Ok(resolved) => {
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
                feature_base_document: Some(feature_base_document),
            })
        }
        Err(_error) => {
            let metadata_path = cli.input.join("metadata.city.json");
            if !metadata_path.is_file() {
                return Err(format!(
                    "{} is neither a cjindex dataset root nor a legacy dataset root containing metadata.city.json",
                    cli.input.display()
                )
                .into());
            }
            Ok(PreparedInput {
                source: parser::InputSource::LegacyFeatureFiles {
                    features_root: cli.input.clone(),
                },
                metadata_path,
                feature_base_document: None,
            })
        }
    }
}

fn derive_base_document(
    city_index: &cityjson_index::CityIndex,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let metadata = city_index.metadata()?;
    let Some(base_document) = metadata.first() else {
        return Err("cjindex dataset does not contain any source metadata".into());
    };
    if metadata
        .iter()
        .skip(1)
        .any(|candidate| candidate.as_ref() != base_document.as_ref())
    {
        return Err(
            "cjindex dataset contains multiple metadata documents; tyler requires one shared base document".into(),
        );
    }
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
                source_tile: Some(tile.clone()),
                source_tile_id: Some(tile.id.clone()),
                content_tile_id: tile.id,
                feature_ids,
            })
        })
        .collect()
}

fn geographic_implicit_tile_export_jobs(
    world: &parser::World,
    quadtree: &spatial_structs::QuadTree,
    tileset: &formats::cesium3dtiles::Tileset,
    root_region: GeographicBounds,
    transformer: &Proj,
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
                source_tile: None,
                source_tile_id: None,
                content_tile_id,
                feature_ids,
            }
        })
        .collect();
    jobs.sort_by(|lhs, rhs| lhs.content_tile_id.cmp(&rhs.content_tile_id));
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

fn tiles_results_successful_content_tile_ids(
    all_content_tile_ids: &[TileId],
    failed_content_tile_ids: &HashSet<TileId>,
) -> Vec<TileId> {
    all_content_tile_ids
        .iter()
        .filter(|tile_id| !failed_content_tile_ids.contains(*tile_id))
        .cloned()
        .collect()
}

fn read_tile_feature_models(
    world: &parser::World,
    feature_ids: &[usize],
) -> Result<Vec<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let mut models = Vec::with_capacity(feature_ids.len());
    match &world.input_source {
        parser::InputSource::LegacyFeatureFiles { features_root } => {
            for fid in feature_ids {
                let parser::FeatureReference::LegacyPath(relative_path) =
                    &world.features[*fid].reference
                else {
                    return Err("legacy input unexpectedly referenced a cjindex feature".into());
                };
                let feature_path = features_root.join(relative_path);
                models.push(cityjson_lib::json::staged::from_feature_file_with_base(
                    &feature_path,
                    &world.feature_base_document,
                )?);
            }
        }
        parser::InputSource::CjIndexDataset { .. } => {
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
                                return Err(
                                    "cjindex input mixed row references with legacy feature ids"
                                        .into(),
                                );
                            };
                            let model = city_index.get(feature_id)?.ok_or_else(|| {
                                format!("feature {feature_id} could not be resolved from cjindex")
                            })?;
                            models.push(model);
                        }
                        return Ok(models);
                    }
                    parser::FeatureReference::LegacyPath(_) => {
                        return Err(
                            "cjindex input unexpectedly referenced a legacy feature path".into(),
                        );
                    }
                }
            }
            models = parser::World::read_cjindex_features_thread_local(
                &world.input_source,
                &cjindex_refs,
            )?;
        }
    }

    Ok(models)
}

fn build_tile_model_from_feature_ids(
    world: &parser::World,
    feature_ids: &[usize],
    feature_type_lods: &BTreeMap<String, String>,
    include_parent_attributes: bool,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let models = prepare_tile_feature_models(
        world,
        feature_ids,
        feature_type_lods,
        include_parent_attributes,
        false,
    )?;
    if models.is_empty() {
        return Err("tile model preparation removed all CityObjects".into());
    }
    let merged = cityjson_lib::ops::merge(models)?;
    cleanup_and_update_extents(merged)
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
    feature_type_lods: &BTreeMap<String, String>,
    include_parent_attributes: bool,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let models = prepare_tile_feature_models(
        world,
        feature_ids,
        feature_type_lods,
        include_parent_attributes,
        true,
    )?;
    let base_root = cityjson_lib::json::from_slice(&world.feature_base_document)?;
    let mut feature_output = Vec::new();
    cityjson_lib::json::write_cityjsonseq_auto_transform(
        &mut feature_output,
        &base_root,
        models,
        [0.001, 0.001, 0.001],
    )?;
    Ok(feature_output)
}

fn prepare_tile_feature_models(
    world: &parser::World,
    feature_ids: &[usize],
    feature_type_lods: &BTreeMap<String, String>,
    include_parent_attributes: bool,
    cleanup_features: bool,
) -> Result<Vec<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    read_tile_feature_models(world, feature_ids)?
        .into_iter()
        .filter_map(|model| {
            match prepare_feature_model(
                model,
                world,
                feature_type_lods,
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
    world: &parser::World,
    feature_type_lods: &BTreeMap<String, String>,
    include_parent_attributes: bool,
    cleanup_feature: bool,
) -> Result<Option<cityjson_lib::CityModel>, Box<dyn std::error::Error>> {
    let mut model = filter_cityobject_types(model, world.cityobject_types.as_ref())?;
    prune_lod_geometries(&mut model, feature_type_lods)?;
    if include_parent_attributes {
        inherit_parent_attributes(&mut model)?;
    }
    let model = remove_empty_geometry_cityobjects(&model)?;
    if model.cityobjects().is_empty() {
        return Ok(None);
    }
    if cleanup_feature {
        cleanup_and_update_extents(model).map(Some)
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
    inherited_attributes: &mut Vec<(String, cityjson_lib::cityjson::v2_0::OwnedAttributeValue)>,
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

fn filter_cityobject_types(
    model: cityjson_lib::CityModel,
    cityobject_types: Option<&Vec<parser::CityObjectType>>,
) -> Result<cityjson_lib::CityModel, Box<dyn std::error::Error>> {
    let Some(cityobject_types) = cityobject_types else {
        return Ok(model);
    };
    let selected = cityobject_types
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<HashSet<_>>();
    filter_cityjsonfeature_preserving_root(&model, |ctx| {
        selected.contains(&ctx.cityobject().type_cityobject().to_string())
    })
}

fn prune_lod_geometries(
    model: &mut cityjson_lib::CityModel,
    feature_type_lods: &BTreeMap<String, String>,
) -> Result<(), Box<dyn std::error::Error>> {
    if feature_type_lods.is_empty() {
        return Ok(());
    }

    let retained_by_object = model
        .cityobjects()
        .iter()
        .map(|(handle, cityobject)| {
            let feature_type = cityobject.type_cityobject().to_string();
            let retained = cityobject
                .geometry()
                .unwrap_or(&[])
                .iter()
                .copied()
                .filter(|geometry_handle| {
                    geometry_matches_lod(
                        model,
                        *geometry_handle,
                        feature_type_lods.get(&feature_type),
                    )
                })
                .collect::<Vec<_>>();
            (handle, retained)
        })
        .collect::<Vec<_>>();

    for (handle, retained) in retained_by_object {
        let cityobject = model
            .cityobjects_mut()
            .get_mut(handle)
            .ok_or_else(|| format!("missing CityObject handle {handle} during LoD pruning"))?;
        cityobject.clear_geometry();
        for geometry_handle in retained {
            cityobject.add_geometry(geometry_handle);
        }
    }

    Ok(())
}

fn geometry_matches_lod(
    model: &cityjson_lib::CityModel,
    geometry_handle: GeometryHandle,
    selected_lod: Option<&String>,
) -> bool {
    let Some(selected_lod) = selected_lod else {
        return true;
    };
    model
        .get_geometry(geometry_handle)
        .and_then(|geometry| geometry.lod())
        .is_some_and(|lod| lod.to_string() == *selected_lod)
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
    let geometric_error_factor = cli.geometric_error_factor.unwrap();
    if geometric_error_factor < 0.0 {
        return Err("--geometric-error-factor must be non-negative".into());
    }
    let format = Formats::_3DTiles; // override --format
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
    #[allow(unused)]
    let metadata_class: String = match format {
        Formats::_3DTiles => {
            if cli.cesium3dtiles_tileset_only {
                String::new()
            } else if cli.cesium3dtiles_metadata_class.is_none() {
                panic!("metadata_class must be set for writing 3D Tiles")
            } else {
                cli.cesium3dtiles_metadata_class.clone().unwrap()
            }
        }
        Formats::CityJSON => "".to_string(),
    };
    if cli.cesium3dtiles_content_bv_from_tile && !cli.cesium3dtiles_content_add_bv {
        warn!("cesium3dtiles_content_bv_from_tile is true, but cesium3dtiles_content_add_bv is false. The tile content bounding volumes are not going to be added, unless you set --3dtiles-content-add-bv");
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
    if (cli.grid_export || log_enabled!(Level::Debug)) && !debug_data_output_path.exists() {
        fs::create_dir(&debug_data_output_path)?;
    }
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

    let world: parser::World = match debug_data.world {
        None => {
            let prepared_input = prepared_input
                .as_ref()
                .expect("prepared input must exist when world is built from source");
            let mut world = match &prepared_input.feature_base_document {
                Some(feature_base_document) => parser::World::from_cjindex(
                    prepared_input.source.clone(),
                    prepared_input.metadata_path.clone(),
                    feature_base_document.clone(),
                    grid_cellsize,
                    cityobject_types,
                    cli.grid_minz,
                    cli.grid_maxz,
                )?,
                None => parser::World::new(
                    &prepared_input.metadata_path,
                    &cli.input,
                    grid_cellsize,
                    cityobject_types,
                    cli.grid_minz,
                    cli.grid_maxz,
                )?,
            };
            world.index_with_grid()?; // todo input: in general, build a line index
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

    if cli.grid_export {
        info!("Exporting the grid to TSV to {:?}", &debug_data_output_path);
        world.export_grid(cli.grid_export_features, Some(&debug_data_output_path))?;
    }
    if log_enabled!(Level::Debug) {
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

    if cli.grid_export {
        info!(
            "Exporting the quadtree to TSV to {:?}",
            &debug_data_output_path
        );
        quadtree.export(&world, Some(&debug_data_output_path))?;
    }
    if log_enabled!(Level::Debug) {
        debug!(
            "Exporting the quadtree instance to bincode to {:?}",
            &debug_data_output_path
        );
        quadtree.export_bincode(Some("quadtree"), Some(&debug_data_output_path))?;
    }

    // 3D Tiles

    let tileset_path = cli.output.join("tileset.json");
    let subtrees_path = cli.output.join("subtrees");
    let tileset_path_unpruned = cli.output.join("tileset_unpruned.json");
    let subtrees_path_unpruned = cli.output.join("subtrees_unpruned");
    info!("Generating 3D Tiles tileset");
    let root_enu_frame = compute_root_enu_frame(&world, &quadtree)?;
    let mut tileset = formats::cesium3dtiles::Tileset::from_quadtree(
        &quadtree,
        &world,
        geometric_error_factor,
        grid_cellsize,
        cli.grid_minz,
        cli.grid_maxz,
        cli.cesium3dtiles_content_bv_from_tile,
        cli.cesium3dtiles_content_add_bv,
        &root_enu_frame,
    );

    if cli.grid_export {
        info!(
            "Exporting the explicit tileset to TSV files to {:?}",
            &debug_data_output_path
        );
        tileset.export(Some(&debug_data_output_path))?;
    }

    let source_crs = format!("EPSG:{}", world.crs.to_epsg()?);
    let source_to_geographic = Proj::new_known_crs(&source_crs, "EPSG:4979", None)?;
    let root_geographic_bounds =
        geographic_bounds_from_source_bbox(&quadtree.bbox(&world.grid), &source_to_geographic)?;

    let export_jobs = match cli.cesium3dtiles_implicit {
        true => {
            let export_jobs = geographic_implicit_tile_export_jobs(
                &world,
                &quadtree,
                &tileset,
                root_geographic_bounds,
                &source_to_geographic,
            )?;
            let content_tile_ids: Vec<TileId> = export_jobs
                .iter()
                .map(|job| job.content_tile_id.clone())
                .collect();
            let mut tileset_implicit = tileset.clone();
            info!("Converting to geographic implicit tiling");
            let components: Vec<_> = subtrees_path_unpruned
                .components()
                .map(|comp| comp.as_os_str())
                .collect();
            let subtrees_dir_option = components.last().cloned().unwrap().to_str();
            let subtrees = tileset_implicit
                .make_implicit_from_content_tile_ids(&content_tile_ids, subtrees_dir_option);

            if cli.cesium3dtiles_tileset_only || log_enabled!(Level::Debug) {
                info!("Writing unpruned 3D Tiles tileset");
                tileset_implicit.to_file(&tileset_path_unpruned)?;

                info!("Writing unpruned subtrees for implicit tiling");
                fs::create_dir_all(&subtrees_path_unpruned)?;
                for (subtree_id, subtree_bytes) in &subtrees {
                    fs::create_dir_all(
                        subtrees_path_unpruned
                            .join(format!("{}/{}", subtree_id.level, subtree_id.x)),
                    )
                    .unwrap();
                    let out_path = subtrees_path_unpruned
                        .join(&subtree_id.to_string())
                        .with_extension("subtree");
                    let mut subtree_file = File::create(&out_path)
                        .unwrap_or_else(|_| panic!("could not create {:?} for writing", &out_path));
                    if let Err(_e) = subtree_file.write_all(subtree_bytes) {
                        warn!("Failed to write subtree {} content", subtree_id);
                    }
                }
            }

            export_jobs
        }
        false => {
            let export_jobs = explicit_tile_export_jobs(&world, &quadtree, &tileset);

            info!("Writing unpruned 3D Tiles tileset");
            tileset.to_file(&tileset_path_unpruned)?;

            export_jobs
        }
    };

    // Export each tile by merging its selected CityJSONFeature stream in memory.
    let path_output_tiles = cli.output.join("t");
    let path_features_input_dir = cli.output.join("inputs");
    // TODO: need to refactor this parallel loop somehow that it does not only read the
    //  3d tiles tiles, but also works with cityjson output
    if !cli.cesium3dtiles_tileset_only {
        fs::create_dir_all(&path_output_tiles)?;
        info!("Created output directory {:#?}", &path_output_tiles);
        if cli.debug_tile_inputs {
            fs::create_dir_all(&path_features_input_dir)?;
            info!("Created output directory {:#?}", &path_features_input_dir);
        }

        let geometry_placement = cityjson_convert::GeometryPlacement::Enu {
            source_crs: source_crs.clone(),
            ecef_origin: root_enu_frame.ecef_origin,
            east: root_enu_frame.east,
            north: root_enu_frame.north,
            up: root_enu_frame.up,
        };
        let export_options = build_glb_export_options(&cli, geometry_placement, None);
        let feature_type_lods = build_feature_type_lods(&cli);
        let tiles_len = export_jobs.len();
        let all_content_tile_ids: Vec<TileId> = export_jobs
            .iter()
            .map(|job| job.content_tile_id.clone())
            .collect();
        let tiles_failed_iter = export_jobs.into_par_iter().map(|job| {
            if job.feature_ids.is_empty() {
                // The Tileset.prune() method removes the empty tiles from the tileset,
                //  so skipping the tile conversion without failure is ok if it's empty.
                debug!(
                    "Tile is empty ({}), skipping conversion",
                    job.content_tile_id
                );
                return None;
            }
            let tileid_string = job.content_tile_id.to_string();
            let file_name = tileid_string;
            let output_file = path_output_tiles.join(&file_name).with_extension("glb");
            let model = match build_tile_model_from_feature_ids(
                &world,
                &job.feature_ids,
                &feature_type_lods,
                cli.include_parent_attributes,
            ) {
                Ok(model) => model,
                Err(error) => {
                    warn!(
                        "Failed to build CityJSON model for tile {}: {}",
                        job.content_tile_id, error
                    );
                    return Some(job);
                }
            };
            if cli.debug_tile_inputs {
                let cityjsonseq_bytes = match build_tile_debug_cityjsonseq(
                    &world,
                    &job.feature_ids,
                    &feature_type_lods,
                    cli.include_parent_attributes,
                ) {
                    Ok(bytes) => bytes,
                    Err(error) => {
                        warn!(
                            "Failed to build debug CityJSONFeature stream for tile {}: {}",
                            job.content_tile_id, error
                        );
                        return Some(job);
                    }
                };
                if let Err(error) = write_debug_tile_input(
                    &path_features_input_dir,
                    file_name.as_str(),
                    &cityjsonseq_bytes,
                ) {
                    warn!(
                        "Failed to write debug CityJSONFeature stream for tile {}: {}",
                        job.content_tile_id, error
                    );
                    return Some(job);
                }
            }
            let mut tile_export_options = export_options.clone();
            if cli.cesium3dtiles_content_clip_to_tile_bounds {
                if cli.cesium3dtiles_implicit {
                    let tile_bounds =
                        geographic_bounds_for_tile(root_geographic_bounds, &job.content_tile_id);
                    tile_export_options.clip_geographic_region =
                        Some(cityjson_convert::GeographicClipRegion {
                            source_crs: source_crs.clone(),
                            west: tile_bounds.west,
                            south: tile_bounds.south,
                            east: tile_bounds.east,
                            north: tile_bounds.north,
                        });
                } else if let Some(source_tile_id) = &job.source_tile_id {
                    let qtree_nodeid: spatial_structs::QuadTreeNodeId = source_tile_id.into();
                    let qtree_node = quadtree.node(&qtree_nodeid).unwrap_or_else(|| {
                        panic!("did not find tile {} in quadtree", source_tile_id)
                    });
                    tile_export_options.clip_bbox = Some(qtree_node.bbox(&world.grid));
                }
            }
            if let Err(error) =
                cityjson_convert::convert_to_glb(&model, &output_file, &tile_export_options)
            {
                warn!("Tile {} conversion failed: {}", job.content_tile_id, error);
                return Some(job);
            }
            if !output_file.exists() {
                warn!(
                    "Tile {} conversion failed: {} was not created",
                    job.content_tile_id,
                    output_file.display()
                );
                return Some(job);
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
            if log_enabled!(Level::Debug) {
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
        info!("Done");

        info!("Pruning tileset of {} failed tiles", tiles_failed.len());
        for (i, failed) in tiles_failed.iter().enumerate() {
            debug!(
                "{}, removing failed from the tileset: {}",
                i, failed.content_tile_id
            );
        }
        if cli.cesium3dtiles_implicit {
            let failed_content_tile_ids: HashSet<TileId> = tiles_failed
                .iter()
                .map(|failed| failed.content_tile_id.clone())
                .collect();
            let content_tile_ids: Vec<TileId> = tiles_results_successful_content_tile_ids(
                &all_content_tile_ids,
                &failed_content_tile_ids,
            );
            let components: Vec<_> = subtrees_path
                .components()
                .map(|comp| comp.as_os_str())
                .collect();
            let subtrees_dir_option = components.last().cloned().unwrap().to_str();
            let subtrees =
                tileset.make_implicit_from_content_tile_ids(&content_tile_ids, subtrees_dir_option);
            info!("Writing subtrees for implicit tiling");
            fs::create_dir_all(&subtrees_path)?;
            for (subtree_id, subtree_bytes) in subtrees {
                fs::create_dir_all(
                    subtrees_path.join(format!("{}/{}", subtree_id.level, subtree_id.x)),
                )
                .unwrap();
                let out_path = subtrees_path
                    .join(&subtree_id.to_string())
                    .with_extension("subtree");
                let mut subtree_file = File::create(&out_path)
                    .unwrap_or_else(|_| panic!("could not create {:?} for writing", &out_path));
                if let Err(_e) = subtree_file.write_all(&subtree_bytes) {
                    warn!("Failed to write subtree {} content", subtree_id);
                }
            }
        } else {
            let failed_tiles: Vec<Tile> = tiles_failed
                .into_iter()
                .filter_map(|failed| failed.source_tile)
                .collect();
            // Remove tiles that failed the gltf conversion
            tileset.prune(&failed_tiles, &quadtree);
            let available_levels = tileset.available_levels();
            // A five level deep tree is still managable in size.
            if available_levels > 5 {
                // Try to find the split where each child tileset starts to have more tiles in their
                // tree, than the ancestor tree. This way, the main tileset is smaller in size than
                // the child tilesets, so it loads faster. This method is not very accurate, because
                // it doesn't account for the actual number of tiles on each level, it only
                // calculates with the theoretical maximum.
                let mut split_at_level = 0;
                for level in (0..available_levels).rev() {
                    let subtree_depth: u32 = (available_levels - level) as u32;
                    let nr_tiles_subtree = (4_usize.pow(subtree_depth) - 1) / 3;
                    let ancestor_tree_depth: u32 =
                        (available_levels - (available_levels - level)) as u32;
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
            }
        }
        info!("Writing 3D Tiles tileset");
        tileset.to_file(&tileset_path)?;
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

    fn prepare_attribute_inheritance_model(
        model: cityjson_lib::CityModel,
        include_parent_attributes: bool,
    ) -> cityjson_lib::CityModel {
        let cityobject_types = vec![
            parser::CityObjectType::Building,
            parser::CityObjectType::BuildingPart,
        ];
        let mut model = filter_cityobject_types(model, Some(&cityobject_types))
            .expect("type filter should succeed");
        prune_lod_geometries(&mut model, &BTreeMap::new()).expect("LoD pruning should succeed");
        if include_parent_attributes {
            inherit_parent_attributes(&mut model).expect("attribute inheritance should succeed");
        }
        let model =
            remove_empty_geometry_cityobjects(&model).expect("empty object removal should succeed");
        cleanup_and_update_extents(model).expect("cleanup should succeed")
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

    #[test]
    fn prepare_model_prunes_lod_geometry_before_gltf_writer() {
        let mut model = cityjson_lib::json::merge_feature_stream_slice(include_bytes!(
            "../cityjson-convert/tests/data/multi_lod_building_part.city.jsonl"
        ))
        .expect("fixture feature stream should parse");
        let lods = BTreeMap::from([("BuildingPart".to_string(), "2.2".to_string())]);

        prune_lod_geometries(&mut model, &lods).expect("LoD pruning should succeed");
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

    #[test]
    fn prepare_model_copies_parent_attributes_when_enabled() {
        let model = parent_attribute_remapping_fixture(serde_json::json!({}));

        let disabled = prepare_attribute_inheritance_model(model.clone(), false);
        let disabled_feature = feature_json(&disabled);
        assert_eq!(
            feature_attribute_string(&disabled_feature, "building-part", "parent_only"),
            None
        );
        assert_eq!(
            feature_attribute_string(&disabled_feature, "building-part", "shared"),
            None
        );

        let enabled = prepare_attribute_inheritance_model(model, true);
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

    #[test]
    fn prepare_model_keeps_child_attributes_on_conflict() {
        let model = parent_attribute_remapping_fixture(serde_json::json!({
            "child_only": "child",
            "levels": 3,
            "shared": "child"
        }));

        let prepared = prepare_attribute_inheritance_model(model, true);
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

    #[test]
    fn build_tile_model_remaps_parent_attributes_before_glb_conversion() {
        let dataset_dir = unique_test_dir("attribute-inheritance");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.city.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::write(
            &feature_path,
            parent_attribute_remapping_fixture_bytes(serde_json::json!({})),
        )
        .expect("write feature");

        let mut world = parser::World::new(
            &metadata_path,
            &features_dir,
            200,
            Some(vec![
                parser::CityObjectType::Building,
                parser::CityObjectType::BuildingPart,
            ]),
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
        assert_eq!(feature_root_id(&model), Some("building-part".to_string()));
    }

    #[test]
    fn build_tile_model_exports_legacy_features() {
        let dataset_dir = unique_test_dir("legacy");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.city.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::copy(resource_path("3dbag_feature_x71.city.jsonl"), &feature_path)
            .expect("copy feature");

        let mut world = parser::World::new(
            &metadata_path,
            &features_dir,
            200,
            Some(vec![parser::CityObjectType::Building]),
            None,
            None,
        )
        .expect("build legacy world");
        world.index_with_grid().expect("index legacy world");
        let quadtree = build_quadtree(&world);
        let model = build_tile_model(&world, &quadtree).expect("build tile model");
        let feature_ids = collect_tile_feature_ids(&world, &quadtree);
        let ndjson = String::from_utf8(
            build_tile_debug_cityjsonseq(&world, &feature_ids, &BTreeMap::new(), false)
                .expect("build debug cityjsonseq"),
        )
        .expect("debug cityjsonseq utf8");

        assert!(!model.cityobjects().is_empty());
        let mut lines = ndjson.lines();
        let header: Value =
            serde_json::from_str(lines.next().expect("CityJSONSeq header should exist"))
                .expect("CityJSONSeq header should parse");
        let feature: Value =
            serde_json::from_str(lines.next().expect("CityJSONSeq feature should exist"))
                .expect("CityJSONSeq feature should parse");
        assert_eq!(header["type"], "CityJSON");
        assert_eq!(feature["type"], "CityJSONFeature");
        assert_eq!(lines.count(), 0);
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

        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
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

        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            None,
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

        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            Some(vec![parser::CityObjectType::Building]),
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
