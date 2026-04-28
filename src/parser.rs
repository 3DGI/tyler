//! CityJSON parsing and feature indexing built directly on `cityjson-lib`.
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

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::Instant;

use cityjson::v2_0::vertex::VertexIndex as GeometryVertexIndex;
use cityjson_index::{CityIndex, StorageLayout};
use cityjson_lib::{cityjson, json};
use log::{debug, info};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::spatial_structs::{Bbox, Cell, CellId, GridLayout};

const CJINDEX_PAGE_SIZE: usize = 65_536;
const CJINDEX_PARALLEL_CHUNK_SIZE: usize = 2_048;
const LARGE_FEATURE_VERTEX_COUNT_THRESHOLD: usize = 50_000;

thread_local! {
    static CJINDEX_THREAD_LOCAL: RefCell<Option<(PathBuf, CityIndex)>> = const { RefCell::new(None) };
    static UNIQUE_ASSIGNMENT_VERTEX_COUNTS: RefCell<Vec<(CellId, usize)>> = const { RefCell::new(Vec::new()) };
}

#[derive(Serialize, Deserialize)]
pub struct World {
    pub cityobject_types: Option<Vec<CityObjectType>>,
    pub crs: Crs,
    pub features: FeatureSet,
    pub feature_base_document: Vec<u8>,
    pub grid: crate::spatial_structs::SquareGrid,
    pub path_metadata: PathBuf,
    pub input_source: InputSource,
}

struct FeatureInGridCells {
    feature: Feature,
    cells: Vec<(CellId, Cell)>,
}

#[derive(Default)]
struct SelectedGeometryStats {
    bbox: Option<Bbox>,
    centroid: Option<[f64; 2]>,
    selected_vertices: Vec<GeometryVertexIndex<u32>>,
    selected_object_types: Vec<CityObjectType>,
    ignored_object_types: Vec<CityObjectType>,
}

#[derive(Default)]
struct ExtentStats {
    extent: Option<Bbox>,
    nr_features: usize,
    nr_features_ignored: usize,
    cityobject_types_ignored: Vec<CityObjectType>,
}

impl ExtentStats {
    fn add_selected_geometry_stats(&mut self, stats: SelectedGeometryStats) {
        if let Some(model_bbox) = stats.bbox {
            if let Some(current) = self.extent.as_mut() {
                merge_bbox(current, &model_bbox);
            } else {
                self.extent = Some(model_bbox);
            }
            self.nr_features += 1;
        } else {
            self.nr_features_ignored += 1;
            self.extend_ignored_types(stats.ignored_object_types);
        }
    }

    fn merge(&mut self, other: Self) {
        if let Some(other_extent) = other.extent {
            if let Some(current) = self.extent.as_mut() {
                merge_bbox(current, &other_extent);
            } else {
                self.extent = Some(other_extent);
            }
        }
        self.nr_features += other.nr_features;
        self.nr_features_ignored += other.nr_features_ignored;
        self.extend_ignored_types(other.cityobject_types_ignored);
    }

    fn extend_ignored_types(&mut self, ignored_types: Vec<CityObjectType>) {
        for cotype in ignored_types {
            if !self.cityobject_types_ignored.contains(&cotype) {
                self.cityobject_types_ignored.push(cotype);
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputSource {
    pub dataset_root: PathBuf,
    pub index_path: PathBuf,
    pub layout: CjIndexLayout,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CjIndexLayout {
    Ndjson,
    CityJson,
    FeatureFiles,
}

impl CjIndexLayout {
    fn storage_layout(self, dataset_root: &Path) -> StorageLayout {
        match self {
            Self::Ndjson => StorageLayout::Ndjson {
                paths: vec![dataset_root.to_path_buf()],
            },
            Self::CityJson => StorageLayout::CityJson {
                paths: vec![dataset_root.to_path_buf()],
            },
            Self::FeatureFiles => StorageLayout::FeatureFiles {
                root: dataset_root.to_path_buf(),
                metadata_glob: "**/metadata.json".to_owned(),
                feature_glob: "**/*.city.jsonl".to_owned(),
            },
        }
    }
}

impl InputSource {
    pub fn from_cjindex_resolved(resolved: &cityjson_index::ResolvedDataset) -> Self {
        let layout = match resolved.layout {
            cityjson_index::DatasetLayoutKind::Ndjson => CjIndexLayout::Ndjson,
            cityjson_index::DatasetLayoutKind::CityJson => CjIndexLayout::CityJson,
            cityjson_index::DatasetLayoutKind::FeatureFiles => CjIndexLayout::FeatureFiles,
        };
        Self {
            dataset_root: resolved.dataset_root.clone(),
            index_path: resolved.index_path.clone(),
            layout,
        }
    }

    pub fn open_index(&self) -> Result<CityIndex, Box<dyn std::error::Error>> {
        Ok(CityIndex::open(
            self.layout.storage_layout(&self.dataset_root),
            &self.index_path,
        )?)
    }
}

impl World {
    pub fn from_cjindex(
        input_source: InputSource,
        path_metadata: PathBuf,
        feature_base_document: Vec<u8>,
        cellsize: u32,
        cityobject_types: Option<Vec<CityObjectType>>,
        arg_minz: Option<i32>,
        arg_maxz: Option<i32>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let metadata = json::from_slice(&feature_base_document)?;
        let crs = Crs::from_model(&metadata)?;

        info!(
            "Computing extent from the features of type {:?}",
            cityobject_types
        );

        let (extent, nr_features, nr_features_ignored, cityobject_types_ignored) =
            if cityobject_types.is_none() {
                let city_index = input_source.open_index()?;
                Self::extent_from_cjindex_bbox_pages(&city_index)?
            } else {
                Self::extent_from_cjindex_features(&input_source, cityobject_types.as_ref())?
            };

        let mut extent = extent.ok_or_else(|| {
            format!(
                "Did not find any CityJSONFeatures of type {:?} in cjindex dataset",
                cityobject_types
            )
        })?;

        if let Some(minz) = arg_minz {
            if extent[2] < minz as f64 {
                extent[2] = minz as f64;
            }
        }
        if let Some(maxz) = arg_maxz {
            if extent[5] > maxz as f64 {
                extent[5] = maxz as f64;
            }
        }

        info!(
            "Found {} features of type {:?}",
            nr_features, &cityobject_types
        );
        info!(
            "Ignored {} features of type {:?}",
            nr_features_ignored, &cityobject_types_ignored
        );
        debug!("extent: {:?}", &extent);
        info!(
            "Computed extent from features: {}",
            crate::spatial_structs::bbox_to_wkt(&extent)
        );

        let grid = crate::spatial_structs::SquareGrid::new(&extent, cellsize, crs.to_epsg()?);
        debug!("{}", grid);

        Ok(Self {
            features: Vec::with_capacity(nr_features),
            crs,
            feature_base_document,
            grid,
            cityobject_types,
            path_metadata,
            input_source,
        })
    }

    fn extent_from_cjindex_bbox_pages(
        city_index: &CityIndex,
    ) -> Result<(Option<Bbox>, usize, usize, Vec<CityObjectType>), Box<dyn std::error::Error>> {
        let Some(summary) = city_index.feature_bounds_summary()? else {
            return Ok((None, 0, 0, Vec::new()));
        };
        Ok((
            Some(Self::cjindex_bounds_to_world_bbox(&summary.bounds)),
            summary.feature_count,
            0,
            Vec::new(),
        ))
    }

    fn extent_from_cjindex_features(
        input_source: &InputSource,
        cityobject_types: Option<&Vec<CityObjectType>>,
    ) -> Result<(Option<Bbox>, usize, usize, Vec<CityObjectType>), Box<dyn std::error::Error>> {
        let city_index = input_source.open_index()?;
        let (chunk_tx, chunk_rx) =
            std::sync::mpsc::sync_channel::<Vec<cityjson_index::IndexedFeatureRef>>(64);

        // Same 2-stage pipeline as `index_with_grid`: a dedicated page-loader
        // thread streams owned chunks to `chunk_tx` so workers can start
        // processing before all pages have been read.
        let total = std::thread::scope(|s| -> Result<ExtentStats, std::io::Error> {
            let page_loader = s.spawn(move || -> Result<(), std::io::Error> {
                let pages_iter = city_index
                    .iter_all_feature_ref_pages(CJINDEX_PAGE_SIZE)
                    .map_err(|e| std::io::Error::other(e.to_string()))?;
                for page_result in pages_iter {
                    let page = page_result.map_err(|e| std::io::Error::other(e.to_string()))?;
                    for chunk in page.chunks(CJINDEX_PARALLEL_CHUNK_SIZE) {
                        if chunk_tx.send(chunk.to_vec()).is_err() {
                            return Ok(());
                        }
                    }
                }
                Ok(())
            });

            let total = chunk_rx
                .into_iter()
                .par_bridge()
                .map(|chunk| {
                    Self::extent_from_cjindex_feature_refs_chunk(
                        input_source,
                        &chunk,
                        cityobject_types,
                    )
                })
                .try_reduce(ExtentStats::default, |mut a, b| {
                    a.merge(b);
                    Ok(a)
                })?;

            page_loader.join().expect("page loader thread panicked")?;
            Ok(total)
        })?;

        Ok((
            total.extent,
            total.nr_features,
            total.nr_features_ignored,
            total.cityobject_types_ignored,
        ))
    }

    fn extent_from_cjindex_feature_refs_chunk(
        input_source: &InputSource,
        feature_refs: &[cityjson_index::IndexedFeatureRef],
        cityobject_types: Option<&Vec<CityObjectType>>,
    ) -> Result<ExtentStats, std::io::Error> {
        let models = Self::read_cjindex_features_thread_local(input_source, feature_refs)
            .map_err(|error| std::io::Error::other(error.to_string()))?;
        let mut chunk = ExtentStats::default();
        for model in models {
            chunk.add_selected_geometry_stats(selected_geometry_stats(&model, cityobject_types));
        }
        Ok(chunk)
    }

    fn cjindex_bounds_to_world_bbox(bounds: &cityjson_index::FeatureBounds) -> Bbox {
        [
            bounds.min_x,
            bounds.min_y,
            bounds.min_z,
            bounds.max_x,
            bounds.max_y,
            bounds.max_z,
        ]
    }

    pub(crate) fn read_cjindex_features_thread_local(
        input_source: &InputSource,
        features: &[cityjson_index::IndexedFeatureRef],
    ) -> cityjson_lib::Result<Vec<cityjson_lib::CityModel>> {
        let index_path = &input_source.index_path;

        CJINDEX_THREAD_LOCAL.with(|cell| {
            let needs_open = {
                let slot = cell.borrow();
                match slot.as_ref() {
                    Some((cached_index_path, _)) => cached_index_path != index_path,
                    None => true,
                }
            };

            if needs_open {
                let city_index = input_source.open_index().map_err(|error| {
                    cityjson_lib::Error::Io(std::io::Error::other(error.to_string()))
                })?;
                *cell.borrow_mut() = Some((index_path.clone(), city_index));
            }

            let slot = cell.borrow();
            let Some((_, city_index)) = slot.as_ref() else {
                return Err(cityjson_lib::Error::Io(std::io::Error::other(
                    "cjindex thread-local index cache was not initialized",
                )));
            };
            city_index.read_features(features)
        })
    }

    pub fn index_with_grid(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        info!("Counting vertices in grid cells");

        self.features.clear();
        let started = Instant::now();
        let city_index = self.input_source.open_index()?;

        // Split the &mut self borrow into disjoint field borrows so workers can
        // hold immutable views (input_source, grid_layout, cityobject_types)
        // while the integrator mutates `features` and `grid`.
        let features: &mut FeatureSet = &mut self.features;
        let grid: &mut crate::spatial_structs::SquareGrid = &mut self.grid;
        let input_source: &InputSource = &self.input_source;
        let cityobject_types: Option<&Vec<CityObjectType>> = self.cityobject_types.as_ref();
        let grid_layout = grid.layout();

        let (chunk_tx, chunk_rx) =
            std::sync::mpsc::sync_channel::<Vec<cityjson_index::IndexedFeatureRef>>(64);
        let (result_tx, result_rx) = std::sync::mpsc::sync_channel::<
            Result<Vec<Option<FeatureInGridCells>>, std::io::Error>,
        >(64);

        let mut indexed_features = 0usize;
        let mut consumer_err: Option<std::io::Error> = None;

        // 3-stage pipeline:
        //   Stage 1 (page_loader): streams pages from cjindex and pushes owned
        //     chunks to chunk_tx. Runs on its own OS thread so I/O overlaps with
        //     CPU work in the rayon pool.
        //   Stage 2 (workers): par_bridge consumes chunk_rx and dispatches each
        //     chunk to a rayon worker, which parses the features, counts vertices
        //     in cells, and pushes the result to result_tx.
        //   Stage 3 (integrator, main thread): drains result_rx and applies the
        //     mutations to `features` and `grid`.
        let loader_outcome = std::thread::scope(|s| -> Result<(usize, usize), std::io::Error> {
            let page_loader = s.spawn(move || -> Result<(usize, usize), std::io::Error> {
                let pages_iter = city_index
                    .iter_all_feature_ref_pages(CJINDEX_PAGE_SIZE)
                    .map_err(|e| std::io::Error::other(e.to_string()))?;
                let mut page_count = 0usize;
                let mut scanned_features = 0usize;
                for page_result in pages_iter {
                    let page = page_result.map_err(|e| std::io::Error::other(e.to_string()))?;
                    page_count += 1;
                    scanned_features += page.len();
                    for chunk in page.chunks(CJINDEX_PARALLEL_CHUNK_SIZE) {
                        if chunk_tx.send(chunk.to_vec()).is_err() {
                            return Ok((page_count, scanned_features));
                        }
                    }
                }
                Ok((page_count, scanned_features))
            });

            s.spawn(move || {
                chunk_rx
                    .into_iter()
                    .par_bridge()
                    .for_each_with(result_tx, |tx, chunk| {
                        let result = Self::index_cjindex_feature_refs_chunk(
                            input_source,
                            &grid_layout,
                            cityobject_types,
                            &chunk,
                        );
                        let _ = tx.send(result);
                    });
            });

            while let Ok(result) = result_rx.recv() {
                match result {
                    Ok(fics) => {
                        for fic in fics.into_iter().flatten() {
                            indexed_features += 1;
                            integrate_feature_in_cells(features, grid, fic);
                        }
                    }
                    Err(e) => {
                        if consumer_err.is_none() {
                            consumer_err = Some(e);
                        }
                    }
                }
            }

            page_loader.join().expect("page loader thread panicked")
        });

        if let Some(e) = consumer_err {
            return Err(Box::new(e));
        }
        let (page_count, scanned_features) = loader_outcome?;

        info!(
            "Indexed {} of {} scanned features into grid cells across {} pages in {:?}",
            features.len(),
            scanned_features,
            page_count,
            started.elapsed()
        );
        Ok(())
    }

    fn index_cjindex_feature_refs_chunk(
        input_source: &InputSource,
        grid_layout: &GridLayout,
        cityobject_types: Option<&Vec<CityObjectType>>,
        feature_refs: &[cityjson_index::IndexedFeatureRef],
    ) -> Result<Vec<Option<FeatureInGridCells>>, std::io::Error> {
        let models = Self::read_cjindex_features_thread_local(input_source, feature_refs)
            .map_err(|error| std::io::Error::other(error.to_string()))?;
        Ok(feature_refs
            .iter()
            .cloned()
            .zip(models.iter())
            .map(|(feature_ref, model)| {
                Self::index_feature_model(
                    grid_layout,
                    cityobject_types,
                    FeatureReference::CjIndexRef(feature_ref),
                    model,
                )
            })
            .collect())
    }

    fn index_feature_model(
        grid_layout: &GridLayout,
        cityobject_types: Option<&Vec<CityObjectType>>,
        feature_reference: FeatureReference,
        model: &cityjson_lib::CityModel,
    ) -> Option<FeatureInGridCells> {
        let stats_started = Instant::now();
        let stats = selected_geometry_stats(model, cityobject_types);
        debug!(
            "selected_geometry_stats collected {} vertices for {:?} in {:?}",
            stats.selected_vertices.len(),
            feature_reference,
            stats_started.elapsed()
        );
        let bbox = stats.bbox?;
        let centroid = stats.centroid?;

        let uses_unique_assignment =
            selected_types_use_unique_assignment(&stats.selected_object_types);
        let feature = Feature {
            centroid,
            reference: feature_reference,
            bbox,
            needs_type_filter: !stats.ignored_object_types.is_empty(),
        };

        if uses_unique_assignment {
            let count_started = Instant::now();
            let (cellid, nr_vertices) =
                count_unique_assignment_cell(model, &stats.selected_vertices, grid_layout, &bbox)?;
            debug!(
                "Selected unique assignment grid cell {:?} with {} vertices in {:?}",
                cellid,
                nr_vertices,
                count_started.elapsed()
            );
            return Some(Self::feature_to_unique_cell(feature, cellid, nr_vertices));
        }

        let count_started = Instant::now();
        let cell_vtx_cnt =
            count_vertices_in_grid(model, &stats.selected_vertices, grid_layout, &bbox);
        debug!(
            "Counted {} grid cells for {:?} in {:?}",
            cell_vtx_cnt.len(),
            feature.reference,
            count_started.elapsed()
        );
        if cell_vtx_cnt.is_empty() {
            return None;
        }
        Self::feature_to_cells(feature, cell_vtx_cnt)
    }

    fn feature_to_cells(feature: Feature, cell_vtx_cnt: CellCounts) -> Option<FeatureInGridCells> {
        let mut cells: Vec<(CellId, Cell)> = Vec::with_capacity(cell_vtx_cnt.len());

        for (cellid, nr_vertices) in cell_vtx_cnt.iter() {
            cells.push((
                *cellid,
                Cell {
                    feature_ids: Vec::new(),
                    nr_vertices: *nr_vertices,
                },
            ));
        }
        Some(FeatureInGridCells { feature, cells })
    }

    fn feature_to_unique_cell(
        feature: Feature,
        cellid: CellId,
        nr_vertices: usize,
    ) -> FeatureInGridCells {
        FeatureInGridCells {
            feature,
            cells: vec![(
                cellid,
                Cell {
                    feature_ids: Vec::new(),
                    nr_vertices,
                },
            )],
        }
    }

    pub fn export_grid(
        &self,
        export_features: bool,
        output_dir: Option<&Path>,
    ) -> std::io::Result<()> {
        if export_features {
            self.grid.export(Some(&self.features), output_dir)
        } else {
            self.grid.export(None, output_dir)
        }
    }

    pub fn export_bincode(
        &self,
        name: Option<&str>,
        output_dir: Option<&Path>,
    ) -> bincode::Result<()> {
        let file_name: &str = name.unwrap_or("world");
        let file = match output_dir {
            None => File::create(format!("{file_name}.bincode"))?,
            Some(outdir) => File::create(outdir.join(format!("{file_name}.bincode")))?,
        };
        bincode::serialize_into(file, self)
    }
}

fn integrate_feature_in_cells(
    features: &mut FeatureSet,
    grid: &mut crate::spatial_structs::SquareGrid,
    feature_in_cells: FeatureInGridCells,
) {
    let fid = features.len();
    features.push(feature_in_cells.feature);
    for (cellid, cell) in feature_in_cells.cells {
        let grid_cell = grid.cell_mut(&cellid);
        grid_cell.nr_vertices += cell.nr_vertices;
        grid_cell.feature_ids.push(fid);
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Crs(String);

impl Crs {
    fn from_model(
        model: &cityjson::v2_0::OwnedCityModel,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let metadata = model
            .metadata()
            .ok_or_else(|| "CityJSON metadata is missing".to_string())?;
        let reference_system = metadata
            .reference_system()
            .ok_or_else(|| "CityJSON metadata.referenceSystem is missing".to_string())?;
        Ok(Self(reference_system.to_string()))
    }

    pub fn to_epsg(&self) -> Result<u16, Box<dyn std::error::Error>> {
        let parts: Vec<&str> = self.0.split('/').collect();
        if let Some(authority) = parts.get(parts.len() - 3) {
            if *authority != "EPSG" {
                return Err(Box::try_from(format!(
                    "the CRS definition should be EPSG: {}",
                    self.0
                ))
                .unwrap());
            }
        }
        if let Some(code) = parts.last() {
            Ok(code.parse::<u16>()?)
        } else {
            Err(Box::try_from(format!(
                "the CRS definition should contain the EPSG code as its last element: {}",
                self.0
            ))
            .unwrap())
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct Feature {
    pub(crate) centroid: [f64; 2],
    pub reference: FeatureReference,
    pub bbox: Bbox,
    #[serde(default = "default_feature_needs_type_filter")]
    pub needs_type_filter: bool,
}

impl Feature {
    pub fn centroid(&self) -> [f64; 2] {
        self.centroid
    }
}

fn default_feature_needs_type_filter() -> bool {
    true
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum FeatureReference {
    CjIndexRef(cityjson_index::IndexedFeatureRef),
    CjIndexId(String),
}

impl Default for FeatureReference {
    fn default() -> Self {
        Self::CjIndexId(String::new())
    }
}

#[derive(
    Debug, Serialize, Deserialize, clap::ValueEnum, Clone, Copy, Ord, PartialOrd, Eq, PartialEq,
)]
#[clap(rename_all = "PascalCase")]
pub enum CityObjectType {
    Bridge,
    BridgePart,
    BridgeInstallation,
    BridgeConstructiveElement,
    BridgeRoom,
    BridgeFurniture,
    Building,
    BuildingPart,
    BuildingInstallation,
    BuildingConstructiveElement,
    BuildingFurniture,
    BuildingStorey,
    BuildingRoom,
    BuildingUnit,
    CityFurniture,
    CityObjectGroup,
    GenericCityObject,
    LandUse,
    OtherConstruction,
    PlantCover,
    SolitaryVegetationObject,
    TINRelief,
    WaterBody,
    Road,
    Railway,
    Waterway,
    TransportSquare,
    Tunnel,
    TunnelPart,
    TunnelInstallation,
    TunnelConstructiveElement,
    TunnelHollowSpace,
    TunnelFurniture,
}

impl fmt::Display for CityObjectType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub type FeatureSet = Vec<Feature>;

fn selected_geometry_stats(
    model: &cityjson::v2_0::OwnedCityModel,
    filter: Option<&Vec<CityObjectType>>,
) -> SelectedGeometryStats {
    let vertex_scratch_capacity = model.vertices().len().min(4096);
    let mut selected_vertices = Vec::with_capacity(vertex_scratch_capacity);
    let mut geometry_scratch = Vec::with_capacity(vertex_scratch_capacity);
    let mut bbox: Option<Bbox> = None;
    let mut x_sum = 0.0;
    let mut y_sum = 0.0;
    let mut nr_vertices = 0usize;
    let mut selected_object_types = Vec::new();
    let mut ignored_object_types = Vec::new();

    collect_selected_vertex_indices(
        model,
        filter,
        &mut selected_vertices,
        &mut geometry_scratch,
        &mut selected_object_types,
        &mut ignored_object_types,
    );

    for vertex_ref in &selected_vertices {
        let Some(vertex) = model.vertices().get(*vertex_ref) else {
            continue;
        };
        let coordinate = vertex.to_array();
        update_bbox(&mut bbox, coordinate);
        x_sum += coordinate[0];
        y_sum += coordinate[1];
        nr_vertices += 1;
    }

    SelectedGeometryStats {
        bbox,
        centroid: (nr_vertices > 0)
            .then_some([x_sum / nr_vertices as f64, y_sum / nr_vertices as f64]),
        selected_vertices,
        selected_object_types,
        ignored_object_types,
    }
}

fn selected_types_use_unique_assignment(selected_object_types: &[CityObjectType]) -> bool {
    selected_object_types.iter().any(|cotype| {
        // Building-like features are assigned once, to the cell with the
        // highest feature score. Terrain and similar surface features keep
        // their bbox-spanning duplicate assignment.
        matches!(
            cotype,
            CityObjectType::Building | CityObjectType::BuildingPart
        )
    })
}

fn collect_selected_vertex_indices(
    model: &cityjson::v2_0::OwnedCityModel,
    filter: Option<&Vec<CityObjectType>>,
    selected_vertices: &mut Vec<GeometryVertexIndex<u32>>,
    geometry_scratch: &mut Vec<GeometryVertexIndex<u32>>,
    selected_object_types: &mut Vec<CityObjectType>,
    ignored_object_types: &mut Vec<CityObjectType>,
) {
    selected_vertices.clear();
    selected_object_types.clear();
    ignored_object_types.clear();

    for (_id, cityobject) in model.cityobjects().iter() {
        let Some(object_type) = map_cityobject_type(cityobject.type_cityobject()) else {
            continue;
        };
        if is_selected_type(filter, object_type) {
            if !selected_object_types.contains(&object_type) {
                selected_object_types.push(object_type);
            }
            let geometry_handles = cityobject.geometry().unwrap_or(&[]);
            for geometry_handle in geometry_handles {
                let Some(geometry) = model.get_geometry(*geometry_handle) else {
                    continue;
                };
                let Some(indices) = geometry.unique_vertex_indices(geometry_scratch) else {
                    continue;
                };
                selected_vertices.extend_from_slice(indices);
            }
        } else if !ignored_object_types.contains(&object_type) {
            ignored_object_types.push(object_type);
        }
    }

    selected_vertices.sort_unstable();
    selected_vertices.dedup();
}

fn count_vertices_in_grid(
    model: &cityjson::v2_0::OwnedCityModel,
    selected_vertices: &[GeometryVertexIndex<u32>],
    grid: &GridLayout,
    bbox: &Bbox,
) -> CellCounts {
    if selected_vertices.is_empty() {
        return CellCounts::default();
    }

    let vertex_counts = if selected_vertices.len() >= LARGE_FEATURE_VERTEX_COUNT_THRESHOLD {
        count_vertex_cells_parallel(model, selected_vertices, grid)
    } else {
        count_vertex_cells(model, selected_vertices, grid)
    };

    if vertex_counts.is_empty() {
        return CellCounts::default();
    }

    CellCounts::merge_vertex_counts_with_bbox(vertex_counts, grid, bbox)
}

fn count_unique_assignment_cell(
    model: &cityjson::v2_0::OwnedCityModel,
    selected_vertices: &[GeometryVertexIndex<u32>],
    grid: &GridLayout,
    bbox: &Bbox,
) -> Option<(CellId, usize)> {
    if selected_vertices.is_empty() {
        return None;
    }

    if selected_vertices.len() >= LARGE_FEATURE_VERTEX_COUNT_THRESHOLD {
        let vertex_counts = count_vertex_cells_parallel(model, selected_vertices, grid);
        if vertex_counts.is_empty() {
            return None;
        }
        return max_merged_vertex_count_with_bbox(&vertex_counts, grid, bbox);
    }

    UNIQUE_ASSIGNMENT_VERTEX_COUNTS.with_borrow_mut(|vertex_counts| {
        count_vertex_cells_into(model, selected_vertices, grid, vertex_counts);
        if vertex_counts.is_empty() {
            return None;
        }
        max_merged_vertex_count_with_bbox(vertex_counts, grid, bbox)
    })
}

fn count_vertex_cells(
    model: &cityjson::v2_0::OwnedCityModel,
    selected_vertices: &[GeometryVertexIndex<u32>],
    grid: &GridLayout,
) -> Vec<(CellId, usize)> {
    let mut vertex_counts = Vec::with_capacity(selected_vertices.len());
    count_vertex_cells_into(model, selected_vertices, grid, &mut vertex_counts);
    vertex_counts
}

fn count_vertex_cells_into(
    model: &cityjson::v2_0::OwnedCityModel,
    selected_vertices: &[GeometryVertexIndex<u32>],
    grid: &GridLayout,
    vertex_counts: &mut Vec<(CellId, usize)>,
) {
    vertex_counts.clear();
    let missing_capacity = selected_vertices
        .len()
        .saturating_sub(vertex_counts.capacity());
    if missing_capacity > 0 {
        vertex_counts.reserve(missing_capacity);
    }

    let vertices = model.vertices();
    for vertex_ref in selected_vertices {
        let Some(vertex) = vertices.get(*vertex_ref) else {
            continue;
        };
        let point = [vertex.x(), vertex.y()];
        vertex_counts.push((grid.locate_point(&point), 0));
    }
    compress_vertex_counts(vertex_counts);
}

fn compress_vertex_counts(vertex_counts: &mut Vec<(CellId, usize)>) {
    if vertex_counts.is_empty() {
        return;
    }

    vertex_counts.sort_unstable_by_key(|(cellid, _)| *cellid);

    let mut read = 0;
    let mut write = 0;
    while read < vertex_counts.len() {
        let cellid = vertex_counts[read].0;
        let mut nr_vertices = 1usize;
        read += 1;

        while read < vertex_counts.len() && vertex_counts[read].0 == cellid {
            nr_vertices += 1;
            read += 1;
        }

        // Keep the historical scoring semantics: the first vertex in a cell
        // starts from 2, matching BTreeMap::entry(...).or_insert(1) += 1.
        vertex_counts[write] = (cellid, nr_vertices + 1);
        write += 1;
    }

    vertex_counts.truncate(write);
}

fn count_vertex_cells_parallel(
    model: &cityjson::v2_0::OwnedCityModel,
    selected_vertices: &[GeometryVertexIndex<u32>],
    grid: &GridLayout,
) -> Vec<(CellId, usize)> {
    let vertices = model.vertices();
    selected_vertices
        .par_iter()
        .fold(BTreeMap::new, |mut vertex_counts, vertex_ref| {
            let Some(vertex) = vertices.get(*vertex_ref) else {
                return vertex_counts;
            };
            let point = [vertex.x(), vertex.y()];
            increment_vertex_count(&mut vertex_counts, grid.locate_point(&point));
            vertex_counts
        })
        .reduce(BTreeMap::new, merge_vertex_count_maps)
        .into_iter()
        .collect()
}

fn max_merged_vertex_count_with_bbox(
    vertex_counts: &[(CellId, usize)],
    grid: &GridLayout,
    bbox: &Bbox,
) -> Option<(CellId, usize)> {
    let (columns, rows) = grid.intersect_bbox_ranges(bbox);
    let min_column = *columns.start();
    let max_column = *columns.end();
    let min_row = *rows.start();
    let max_row = *rows.end();
    let mut vertex_counts = vertex_counts.iter().copied().peekable();
    let mut best: Option<(CellId, usize)> = None;

    for row in min_row..=max_row {
        for column in min_column..=max_column {
            let bbox_cellid = CellId { row, column };
            while let Some((cellid, count)) =
                vertex_counts.next_if(|(cellid, _)| *cellid < bbox_cellid)
            {
                update_best_cell(&mut best, cellid, count);
            }

            if let Some((_, count)) = vertex_counts.next_if(|(cellid, _)| *cellid == bbox_cellid) {
                update_best_cell(&mut best, bbox_cellid, count + 1);
            } else {
                update_best_cell(&mut best, bbox_cellid, 2);
            }
        }
    }

    for (cellid, count) in vertex_counts {
        update_best_cell(&mut best, cellid, count);
    }

    best
}

fn update_best_cell(best: &mut Option<(CellId, usize)>, cellid: CellId, count: usize) {
    if best.as_ref().is_none_or(|(best_cellid, best_count)| {
        count > *best_count || (count == *best_count && cellid > *best_cellid)
    }) {
        *best = Some((cellid, count));
    }
}

fn increment_vertex_count(vertex_counts: &mut BTreeMap<CellId, usize>, cellid: CellId) {
    *vertex_counts.entry(cellid).or_insert(1) += 1;
}

fn merge_vertex_count_maps(
    mut left: BTreeMap<CellId, usize>,
    right: BTreeMap<CellId, usize>,
) -> BTreeMap<CellId, usize> {
    for (cellid, count) in right {
        match left.get_mut(&cellid) {
            Some(left_count) => *left_count += count - 1,
            None => {
                left.insert(cellid, count);
            }
        }
    }
    left
}

#[derive(Default)]
struct CellCounts {
    entries: Vec<(CellId, usize)>,
}

impl CellCounts {
    fn len(&self) -> usize {
        self.entries.len()
    }

    fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn iter(&self) -> impl Iterator<Item = (&CellId, &usize)> {
        self.entries.iter().map(|(cellid, count)| (cellid, count))
    }

    fn merge_vertex_counts_with_bbox(
        vertex_counts: Vec<(CellId, usize)>,
        grid: &GridLayout,
        bbox: &Bbox,
    ) -> Self {
        let (columns, rows) = grid.intersect_bbox_ranges(bbox);
        let min_column = *columns.start();
        let max_column = *columns.end();
        let min_row = *rows.start();
        let max_row = *rows.end();
        let bbox_len = (max_column - min_column + 1) * (max_row - min_row + 1);
        let mut entries = Vec::with_capacity(vertex_counts.len().max(bbox_len));
        let mut vertex_counts = vertex_counts.into_iter().peekable();

        for row in min_row..=max_row {
            for column in min_column..=max_column {
                let bbox_cellid = CellId { row, column };
                while let Some((cellid, count)) =
                    vertex_counts.next_if(|(cellid, _)| *cellid < bbox_cellid)
                {
                    entries.push((cellid, count));
                }

                if let Some((_, count)) =
                    vertex_counts.next_if(|(cellid, _)| *cellid == bbox_cellid)
                {
                    entries.push((bbox_cellid, count + 1));
                } else {
                    entries.push((bbox_cellid, 2));
                }
            }
        }

        entries.extend(vertex_counts);
        Self { entries }
    }
}

fn is_selected_type(filter: Option<&Vec<CityObjectType>>, object_type: CityObjectType) -> bool {
    match filter {
        Some(types) => types.contains(&object_type),
        None => true,
    }
}

fn map_cityobject_type<SS: cityjson::resources::storage::StringStorage>(
    object_type: &cityjson::v2_0::CityObjectType<SS>,
) -> Option<CityObjectType> {
    use cityjson::v2_0::CityObjectType as CjType;

    Some(match object_type {
        CjType::Bridge => CityObjectType::Bridge,
        CjType::BridgePart => CityObjectType::BridgePart,
        CjType::BridgeInstallation => CityObjectType::BridgeInstallation,
        CjType::BridgeConstructiveElement => CityObjectType::BridgeConstructiveElement,
        CjType::BridgeRoom => CityObjectType::BridgeRoom,
        CjType::BridgeFurniture => CityObjectType::BridgeFurniture,
        CjType::Building => CityObjectType::Building,
        CjType::BuildingPart => CityObjectType::BuildingPart,
        CjType::BuildingInstallation => CityObjectType::BuildingInstallation,
        CjType::BuildingConstructiveElement => CityObjectType::BuildingConstructiveElement,
        CjType::BuildingFurniture => CityObjectType::BuildingFurniture,
        CjType::BuildingStorey => CityObjectType::BuildingStorey,
        CjType::BuildingRoom => CityObjectType::BuildingRoom,
        CjType::BuildingUnit => CityObjectType::BuildingUnit,
        CjType::CityFurniture => CityObjectType::CityFurniture,
        CjType::CityObjectGroup => CityObjectType::CityObjectGroup,
        CjType::GenericCityObject => CityObjectType::GenericCityObject,
        CjType::LandUse => CityObjectType::LandUse,
        CjType::OtherConstruction => CityObjectType::OtherConstruction,
        CjType::PlantCover => CityObjectType::PlantCover,
        CjType::SolitaryVegetationObject => CityObjectType::SolitaryVegetationObject,
        CjType::TINRelief => CityObjectType::TINRelief,
        CjType::WaterBody => CityObjectType::WaterBody,
        CjType::Road => CityObjectType::Road,
        CjType::Railway => CityObjectType::Railway,
        CjType::Waterway => CityObjectType::Waterway,
        CjType::TransportSquare => CityObjectType::TransportSquare,
        CjType::Tunnel => CityObjectType::Tunnel,
        CjType::TunnelPart => CityObjectType::TunnelPart,
        CjType::TunnelInstallation => CityObjectType::TunnelInstallation,
        CjType::TunnelConstructiveElement => CityObjectType::TunnelConstructiveElement,
        CjType::TunnelHollowSpace => CityObjectType::TunnelHollowSpace,
        CjType::TunnelFurniture => CityObjectType::TunnelFurniture,
        CjType::Default | CjType::Extension(_) => return None,
        _ => return None,
    })
}

fn update_bbox(bbox: &mut Option<Bbox>, coordinate: [f64; 3]) {
    match bbox {
        Some(current) => {
            if coordinate[0] < current[0] {
                current[0] = coordinate[0];
            }
            if coordinate[1] < current[1] {
                current[1] = coordinate[1];
            }
            if coordinate[2] < current[2] {
                current[2] = coordinate[2];
            }
            if coordinate[0] > current[3] {
                current[3] = coordinate[0];
            }
            if coordinate[1] > current[4] {
                current[4] = coordinate[1];
            }
            if coordinate[2] > current[5] {
                current[5] = coordinate[2];
            }
        }
        None => {
            *bbox = Some([
                coordinate[0],
                coordinate[1],
                coordinate[2],
                coordinate[0],
                coordinate[1],
                coordinate[2],
            ])
        }
    }
}

fn merge_bbox(target: &mut Bbox, other: &Bbox) {
    if other[0] < target[0] {
        target[0] = other[0];
    }
    if other[1] < target[1] {
        target[1] = other[1];
    }
    if other[2] < target[2] {
        target[2] = other[2];
    }
    if other[3] > target[3] {
        target[3] = other[3];
    }
    if other[4] > target[4] {
        target[4] = other[4];
    }
    if other[5] > target[5] {
        target[5] = other[5];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cityjson_lib::json::staged::from_feature_file_with_base;

    fn resource_path(name: &str) -> PathBuf {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let direct = manifest.join("resources").join("data").join(name);
        if direct.exists() {
            direct
        } else {
            manifest
                .join("..")
                .join("resources")
                .join("data")
                .join(name)
        }
    }

    #[test]
    fn test_crs_to_epsg() {
        let crs = Crs("https://www.opengis.net/def/crs/EPSG/0/7415".to_string());
        let epsg_code = crs.to_epsg().unwrap();
        assert_eq!(7415_u16, epsg_code);
    }

    #[test]
    fn test_feature_file_loads_through_cjlib() {
        let pb = resource_path("3dbag_feature_x71.city.jsonl");
        let base = std::fs::read(resource_path("3dbag_x00.city.json")).unwrap();
        let model = from_feature_file_with_base(pb, &base).unwrap();
        let stats = selected_geometry_stats(&model, Some(&vec![CityObjectType::Building]));
        assert!(stats.bbox.is_some());
        assert!(!stats.selected_vertices.is_empty());
        let bbox = stats.bbox.unwrap();
        assert!(bbox[0] > 0.0);
        assert!(bbox[1] > 0.0);
    }

    fn parse_feature(feature: serde_json::Value) -> cityjson::v2_0::OwnedCityModel {
        cityjson_lib::json::from_feature_slice(&serde_json::to_vec(&feature).unwrap()).unwrap()
    }

    fn unique_assignment_winner_from_full_counts(counts: &CellCounts) -> Option<(CellId, usize)> {
        counts
            .iter()
            .max_by(|left, right| left.1.cmp(right.1))
            .map(|(cellid, count)| (*cellid, *count))
    }

    fn assert_unique_assignment_matches_full_counts(
        model: &cityjson::v2_0::OwnedCityModel,
        object_type: CityObjectType,
        grid: &crate::spatial_structs::SquareGrid,
    ) -> (CellId, usize) {
        let stats = selected_geometry_stats(model, Some(&vec![object_type]));
        let bbox = stats.bbox.unwrap();
        let layout = grid.layout();
        let full_counts = count_vertices_in_grid(model, &stats.selected_vertices, &layout, &bbox);
        let expected = unique_assignment_winner_from_full_counts(&full_counts);
        let optimized =
            count_unique_assignment_cell(model, &stats.selected_vertices, &layout, &bbox);

        assert_eq!(optimized, expected);
        optimized.expect("unique assignment should produce a winning cell")
    }

    #[test]
    fn unique_assignment_single_cell_matches_full_cell_counts() {
        let model = parse_feature(serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2]]]
                    }]
                }
            },
            "vertices": [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        }));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 100.0, 100.0, 10.0],
            100,
            7415,
        );

        let (cellid, count) =
            assert_unique_assignment_matches_full_counts(&model, CityObjectType::Building, &grid);

        assert_eq!(cellid, grid.locate_point(&[0.0, 0.0]));
        assert_eq!(count, 5);
    }

    #[test]
    fn unique_assignment_multi_cell_counts_match_full_cell_counts() {
        let model = parse_feature(serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2, 3, 4]]]
                    }]
                }
            },
            "vertices": [
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [150, 0, 0],
                [250, 0, 0]
            ]
        }));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 300.0, 100.0, 10.0],
            100,
            7415,
        );

        let (cellid, count) =
            assert_unique_assignment_matches_full_counts(&model, CityObjectType::Building, &grid);

        assert_eq!(cellid, grid.locate_point(&[0.0, 0.0]));
        assert_eq!(count, 5);
    }

    #[test]
    fn unique_assignment_equal_count_tie_matches_full_cell_counts() {
        let model = parse_feature(serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2, 3]]]
                    }]
                }
            },
            "vertices": [
                [0, 0, 0],
                [1, 0, 0],
                [150, 0, 0],
                [151, 0, 0]
            ]
        }));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 200.0, 100.0, 10.0],
            100,
            7415,
        );

        let (cellid, count) =
            assert_unique_assignment_matches_full_counts(&model, CityObjectType::Building, &grid);

        assert_eq!(cellid, grid.locate_point(&[150.0, 0.0]));
        assert_eq!(count, 4);
    }

    #[test]
    fn unique_assignment_large_bbox_keeps_bbox_only_cells_lower_scored() {
        let model = parse_feature(serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2]]]
                    }]
                }
            },
            "vertices": [[0, 0, 0], [250, 250, 0], [0, 250, 0]]
        }));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 300.0, 300.0, 10.0],
            100,
            7415,
        );

        let (cellid, count) =
            assert_unique_assignment_matches_full_counts(&model, CityObjectType::Building, &grid);

        assert_eq!(cellid, grid.locate_point(&[250.0, 250.0]));
        assert_eq!(count, 3);
    }

    #[test]
    fn unique_assignment_building_part_matches_full_cell_counts() {
        let model = parse_feature(serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-part-1",
            "CityObjects": {
                "building-part-1": {
                    "type": "BuildingPart",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2]]]
                    }]
                }
            },
            "vertices": [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        }));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 100.0, 100.0, 10.0],
            100,
            7415,
        );

        let (cellid, count) = assert_unique_assignment_matches_full_counts(
            &model,
            CityObjectType::BuildingPart,
            &grid,
        );

        assert_eq!(cellid, grid.locate_point(&[0.0, 0.0]));
        assert_eq!(count, 5);
    }

    #[test]
    fn count_vertices_in_grid_only_counts_selected_vertices() {
        let base = serde_json::to_vec(&serde_json::json!({
            "type": "CityJSON",
            "version": "2.0",
            "CityObjects": {},
            "vertices": [],
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [0.0, 0.0, 0.0]
            }
        }))
        .unwrap();
        let feature = serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2]]]
                    }]
                },
                "road-1": {
                    "type": "Road",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[3, 4, 5]]]
                    }]
                }
            },
            "vertices": [
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [100, 100, 0],
                [101, 100, 0],
                [100, 101, 0]
            ]
        });
        let feature_path = std::env::temp_dir().join(format!(
            "tyler-parser-selected-vertices-{}.city.jsonl",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&feature_path, serde_json::to_vec(&feature).unwrap()).unwrap();

        let model = from_feature_file_with_base(&feature_path, &base).unwrap();
        let stats = selected_geometry_stats(&model, Some(&vec![CityObjectType::Building]));
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 200.0, 200.0, 10.0],
            100,
            7415,
        );
        let counts = count_vertices_in_grid(
            &model,
            &stats.selected_vertices,
            &grid.layout(),
            &stats.bbox.unwrap(),
        );

        let building_cell = grid.locate_point(&[0.0, 0.0]);
        let road_cell = grid.locate_point(&[100.0, 100.0]);

        assert!(counts.iter().any(|(cellid, _)| *cellid == building_cell));
        assert!(!counts.iter().any(|(cellid, _)| *cellid == road_cell));

        let _ = std::fs::remove_file(feature_path);
    }

    #[test]
    fn count_vertices_in_grid_matches_reference_hashmap_counts() {
        let base = serde_json::to_vec(&serde_json::json!({
            "type": "CityJSON",
            "version": "2.0",
            "CityObjects": {},
            "vertices": [],
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [0.0, 0.0, 0.0]
            }
        }))
        .unwrap();
        let feature = serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building",
                    "geometry": [{
                        "type": "MultiSurface",
                        "lod": "1.0",
                        "boundaries": [[[0, 1, 2, 3]]]
                    }]
                }
            },
            "vertices": [
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [250, 250, 0]
            ]
        });
        let feature_path = std::env::temp_dir().join(format!(
            "tyler-parser-grid-counts-{}.city.jsonl",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&feature_path, serde_json::to_vec(&feature).unwrap()).unwrap();

        let model = from_feature_file_with_base(&feature_path, &base).unwrap();
        let stats = selected_geometry_stats(&model, Some(&vec![CityObjectType::Building]));
        let bbox = stats.bbox.unwrap();
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 300.0, 300.0, 10.0],
            100,
            7415,
        );

        let layout = grid.layout();
        let optimized = count_vertices_in_grid(&model, &stats.selected_vertices, &layout, &bbox);
        let reference =
            reference_count_vertices_in_grid(&model, &stats.selected_vertices, &grid, &bbox);

        assert_eq!(
            optimized
                .iter()
                .map(|(cellid, count)| (*cellid, *count))
                .collect::<Vec<_>>(),
            reference
        );

        let optimized = optimized
            .iter()
            .map(|(cellid, count)| (*cellid, *count))
            .collect::<BTreeMap<_, _>>();
        let same_cell = grid.locate_point(&[0.0, 0.0]);
        let bbox_only_cell = grid.locate_point(&[200.0, 0.0]);
        assert_eq!(optimized.get(&same_cell), Some(&5));
        assert_eq!(optimized.get(&bbox_only_cell), Some(&2));

        let _ = std::fs::remove_file(feature_path);
    }

    #[test]
    fn count_vertices_in_grid_large_feature_path_matches_reference_hashmap_counts() {
        let base = serde_json::to_vec(&serde_json::json!({
            "type": "CityJSON",
            "version": "2.0",
            "CityObjects": {},
            "vertices": [],
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [0.0, 0.0, 0.0]
            }
        }))
        .unwrap();
        let feature = serde_json::json!({
            "type": "CityJSONFeature",
            "id": "building-1",
            "CityObjects": {
                "building-1": {
                    "type": "Building"
                }
            },
            "vertices": [
                [0, 0, 0],
                [1, 0, 0],
                [250, 250, 0]
            ]
        });
        let feature_path = std::env::temp_dir().join(format!(
            "tyler-parser-large-grid-counts-{}.city.jsonl",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&feature_path, serde_json::to_vec(&feature).unwrap()).unwrap();

        let model = from_feature_file_with_base(&feature_path, &base).unwrap();
        let bbox: Bbox = [0.0, 0.0, 0.0, 250.0, 250.0, 0.0];
        let grid = crate::spatial_structs::SquareGrid::new(
            &[0.0, 0.0, 0.0, 300.0, 300.0, 10.0],
            100,
            7415,
        );
        let selected_vertices = (0..LARGE_FEATURE_VERTEX_COUNT_THRESHOLD)
            .map(|index| {
                GeometryVertexIndex::new(match index % 4 {
                    0 | 3 => 0,
                    1 => 1,
                    _ => 2,
                })
            })
            .collect::<Vec<_>>();

        let layout = grid.layout();
        let optimized = count_vertices_in_grid(&model, &selected_vertices, &layout, &bbox);
        let reference = reference_count_vertices_in_grid(&model, &selected_vertices, &grid, &bbox);

        assert_eq!(
            optimized
                .iter()
                .map(|(cellid, count)| (*cellid, *count))
                .collect::<Vec<_>>(),
            reference
        );

        let _ = std::fs::remove_file(feature_path);
    }

    fn reference_count_vertices_in_grid(
        model: &cityjson::v2_0::OwnedCityModel,
        selected_vertices: &[GeometryVertexIndex<u32>],
        grid: &crate::spatial_structs::SquareGrid,
        bbox: &Bbox,
    ) -> Vec<(CellId, usize)> {
        let mut cell_vtx_cnt = std::collections::BTreeMap::new();

        for vertex_ref in selected_vertices {
            let Some(vertex) = model.vertices().get(*vertex_ref) else {
                continue;
            };
            let point = [vertex.x(), vertex.y()];
            let cellid = grid.locate_point(&point);
            *cell_vtx_cnt.entry(cellid).or_insert(1) += 1;
        }

        if !selected_vertices.is_empty() {
            for cellid in grid.intersect_bbox(bbox) {
                *cell_vtx_cnt.entry(cellid).or_insert(1) += 1;
            }
        }

        cell_vtx_cnt.into_iter().collect()
    }
}
