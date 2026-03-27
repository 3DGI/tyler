//! CityJSON parsing and feature indexing built directly on `cjlib`.
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

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};

use cjlib::{cityjson, json};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use walkdir::WalkDir;

use crate::spatial_structs::{Bbox, Cell, CellId};

#[derive(Serialize, Deserialize)]
pub struct World {
    pub cityobject_types: Option<Vec<CityObjectType>>,
    pub crs: Crs,
    pub features: FeatureSet,
    pub feature_base_document: Vec<u8>,
    pub grid: crate::spatial_structs::SquareGrid,
    pub path_features_root: PathBuf,
    pub path_metadata: PathBuf,
}

struct ExtentResult {
    extent: Bbox,
    nr_features: usize,
    cityobject_types_ignored: Vec<CityObjectType>,
    nr_features_ignored: usize,
}

struct FeatureDirsFiles {
    feature_dirs: Vec<PathBuf>,
    feature_files: Vec<PathBuf>,
}

struct FeatureInGridCells {
    feature: Feature,
    cells: Vec<(CellId, Cell)>,
}

#[derive(Default)]
struct SelectedGeometryStats {
    bbox: Option<Bbox>,
    centroid: Option<[f64; 2]>,
    nr_vertices: usize,
    selected_object_types: Vec<CityObjectType>,
    ignored_object_types: Vec<CityObjectType>,
}

impl World {
    pub fn new<P: AsRef<Path>>(
        path_metadata: P,
        path_features_root: P,
        cellsize: u32,
        cityobject_types: Option<Vec<CityObjectType>>,
        arg_minz: Option<i32>,
        arg_maxz: Option<i32>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let path_features_root = path_features_root.as_ref().to_path_buf();
        let path_metadata = path_metadata.as_ref().to_path_buf();
        let feature_base_document = std::fs::read(&path_metadata)?;
        let metadata = json::from_slice(&feature_base_document)?;
        let crs = Crs::from_model(metadata.as_inner())?;

        info!(
            "Computing extent from the features of type {:?}",
            cityobject_types
        );

        let feature_dirs_files = Self::find_feature_dirs_and_files(&path_features_root);
        debug!(
            "Found {} subdirectories and {} CityJSONFeature files at the root directory",
            feature_dirs_files.feature_dirs.len(),
            feature_dirs_files.feature_files.len()
        );

        let extents: Vec<ExtentResult> = feature_dirs_files
            .feature_dirs
            .into_par_iter()
            .filter_map(|dir| {
                Self::extent(dir, cityobject_types.as_ref(), feature_base_document.as_slice())
            })
            .collect();

        let mut nr_features = 0;
        let mut nr_features_ignored = 0;
        let mut extent = Self::extent_init(
            &path_features_root,
            cityobject_types.as_ref(),
            feature_base_document.as_slice(),
        )
            .unwrap_or_else(|| {
                panic!(
                    "Did not find any CityJSONFeature of type {:?} in {}",
                    cityobject_types,
                    path_features_root.display()
                )
            });
        let mut cityobject_types_ignored: Vec<CityObjectType> = Vec::new();

        for (i, extent_result) in extents.iter().enumerate() {
            nr_features += extent_result.nr_features;
            nr_features_ignored += extent_result.nr_features_ignored;
            if i == 0 {
                extent = extent_result.extent;
            } else {
                merge_bbox(&mut extent, &extent_result.extent);
            }
            for cotype in &extent_result.cityobject_types_ignored {
                if !cityobject_types_ignored.contains(cotype) {
                    cityobject_types_ignored.push(*cotype);
                }
            }
        }

        for feature_path in &feature_dirs_files.feature_files {
            Self::extent_file(
                cityobject_types.as_ref(),
                &mut extent,
                &mut nr_features,
                &mut nr_features_ignored,
                &mut cityobject_types_ignored,
                feature_path,
                feature_base_document.as_slice(),
            );
        }

        if nr_features == 0 {
            panic!(
                "Did not find any CityJSONFeatures of type {:?}",
                cityobject_types
            );
        }

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
            path_features_root,
            path_metadata,
        })
    }

    fn find_feature_dirs_and_files(path_features_root: &PathBuf) -> FeatureDirsFiles {
        let mut path_features_root_dirs: Vec<PathBuf> = Vec::new();
        let mut path_features_root_files: Vec<PathBuf> = Vec::new();
        for entry_res in WalkDir::new(path_features_root).min_depth(1).max_depth(1) {
            if let Ok(entry) = entry_res {
                if entry.file_type().is_dir() {
                    path_features_root_dirs.push(entry.path().to_path_buf());
                } else if entry.file_type().is_file() {
                    if let Some(jsonl_path) = Self::direntry_to_jsonl(entry) {
                        path_features_root_files.push(jsonl_path)
                    }
                }
            } else {
                error!(
                    "Error in walking the directory {}, error: {}",
                    &path_features_root.display(),
                    entry_res.unwrap_err()
                )
            }
        }
        FeatureDirsFiles {
            feature_dirs: path_features_root_dirs,
            feature_files: path_features_root_files,
        }
    }

    fn extent<P: AsRef<Path> + std::fmt::Debug>(
        path_features: P,
        cityobject_types: Option<&Vec<CityObjectType>>,
        feature_base_document: &[u8],
    ) -> Option<ExtentResult> {
        let features_enum_iter = WalkDir::new(&path_features)
            .into_iter()
            .filter_map(Self::jsonl_path);

        let mut extent: Option<Bbox> = None;
        let mut nr_features = 0;
        let mut nr_features_ignored = 0;
        let mut cityobject_types_ignored: Vec<CityObjectType> = Vec::new();

        for feature_path in features_enum_iter {
            match json::from_feature_file_with_base(&feature_path, feature_base_document) {
                Ok(feature) => {
                    let stats = selected_geometry_stats(feature.as_inner(), cityobject_types);
                    if let Some(bbox) = stats.bbox {
                        if let Some(accumulated) = extent.as_mut() {
                            merge_bbox(accumulated, &bbox);
                        } else {
                            extent = Some(bbox);
                        }
                        nr_features += 1;
                    } else {
                        nr_features_ignored += 1;
                        for cotype in stats.ignored_object_types {
                            if !cityobject_types_ignored.contains(&cotype) {
                                cityobject_types_ignored.push(cotype);
                            }
                        }
                    }
                }
                Err(e) => warn!("Failed to parse {:?} with {:?}", &feature_path, e),
            }
        }

        extent.map(|extent| ExtentResult {
            extent,
            nr_features,
            cityobject_types_ignored,
            nr_features_ignored,
        })
    }

    fn extent_init<P: AsRef<Path> + std::fmt::Debug>(
        path_features: P,
        cityobject_types: Option<&Vec<CityObjectType>>,
        feature_base_document: &[u8],
    ) -> Option<Bbox> {
        let features_enum_iter = WalkDir::new(&path_features)
            .into_iter()
            .filter_map(Self::jsonl_path);

        for feature_path in features_enum_iter {
            match json::from_feature_file_with_base(&feature_path, feature_base_document) {
                Ok(feature) => {
                    let stats = selected_geometry_stats(feature.as_inner(), cityobject_types);
                    if let Some(bbox) = stats.bbox {
                        return Some(bbox);
                    }
                }
                Err(e) => warn!("Failed to parse {:?} with {:?}", &feature_path, e),
            }
        }
        None
    }

    fn extent_file(
        cityobject_types: Option<&Vec<CityObjectType>>,
        extent: &mut Bbox,
        nr_features: &mut usize,
        nr_features_ignored: &mut usize,
        cityobject_types_ignored: &mut Vec<CityObjectType>,
        feature_path: &PathBuf,
        feature_base_document: &[u8],
    ) {
        if let Ok(feature) = json::from_feature_file_with_base(feature_path, feature_base_document)
        {
            let stats = selected_geometry_stats(feature.as_inner(), cityobject_types);
            if let Some(bbox) = stats.bbox {
                merge_bbox(extent, &bbox);
                *nr_features += 1;
            } else {
                *nr_features_ignored += 1;
                for cotype in stats.ignored_object_types {
                    if !cityobject_types_ignored.contains(&cotype) {
                        cityobject_types_ignored.push(cotype);
                    }
                }
            }
        } else {
            error!("Failed to parse {:?}", &feature_path);
        }
    }

    pub fn jsonl_path(walkdir_res: Result<walkdir::DirEntry, walkdir::Error>) -> Option<PathBuf> {
        if let Ok(entry) = walkdir_res {
            Self::direntry_to_jsonl(entry)
        } else {
            None
        }
    }

    fn direntry_to_jsonl(entry: walkdir::DirEntry) -> Option<PathBuf> {
        if let Some(ext) = entry.path().extension() {
            if ext == "jsonl" {
                Some(entry.path().to_path_buf())
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn index_with_grid(&mut self) {
        let feature_dirs_files = Self::find_feature_dirs_and_files(&self.path_features_root);
        info!("Counting vertices in grid cells");

        self.features.clear();
        for dir in feature_dirs_files.feature_dirs {
            for feature_path in WalkDir::new(dir).into_iter().filter_map(Self::jsonl_path) {
                if let Some(feature_in_cells) = self.index_feature_path(&feature_path) {
                    self.integrate_feature_in_cells(feature_in_cells);
                }
            }
        }
        for feature_path in feature_dirs_files.feature_files {
            if let Some(feature_in_cells) = self.index_feature_path(&feature_path) {
                self.integrate_feature_in_cells(feature_in_cells);
            }
        }
        debug!("indexed {} features", self.features.len());
    }

    fn index_feature_path(&self, feature_path: &PathBuf) -> Option<FeatureInGridCells> {
        let feature = json::from_feature_file_with_base(feature_path, &self.feature_base_document);
        if let Ok(model) = feature {
            let stats = selected_geometry_stats(model.as_inner(), self.cityobject_types.as_ref());
            let bbox = stats.bbox?;
            let centroid = stats.centroid?;
            let cell_vtx_cnt = count_vertices_in_grid(
                model.as_inner(),
                self.cityobject_types.as_ref(),
                &self.grid,
                &bbox,
            );
            if cell_vtx_cnt.is_empty() {
                return None;
            }
            self.feature_to_cells(
                feature_path,
                Feature {
                    centroid,
                    path_jsonl: feature_path
                        .strip_prefix(&self.path_features_root)
                        .unwrap_or(feature_path)
                        .to_path_buf(),
                    bbox,
                },
                cell_vtx_cnt,
                &stats.selected_object_types,
            )
        } else {
            error!("Failed to parse the feature {:?}", &feature_path);
            None
        }
    }

    fn integrate_feature_in_cells(&mut self, feature_in_cells: FeatureInGridCells) {
        let fid = self.features.len();
        self.features.push(feature_in_cells.feature);
        for (cellid, cell) in feature_in_cells.cells {
            let grid_cell = self.grid.cell_mut(&cellid);
            grid_cell.nr_vertices += cell.nr_vertices;
            if !grid_cell.feature_ids.contains(&fid) {
                grid_cell.feature_ids.push(fid)
            }
        }
    }

    fn feature_to_cells(
        &self,
        feature_path: &PathBuf,
        feature: Feature,
        cell_vtx_cnt: HashMap<CellId, usize>,
        selected_object_types: &[CityObjectType],
    ) -> Option<FeatureInGridCells> {
        let _ = feature_path;
        let unique_assignment = selected_object_types.iter().any(|cotype| {
            matches!(cotype, CityObjectType::Building | CityObjectType::BuildingPart)
        });
        let mut cells: Vec<(CellId, Cell)> = Vec::with_capacity(cell_vtx_cnt.len());

        if unique_assignment {
            let (cellid, nr_vertices) = cell_vtx_cnt
                .iter()
                .max_by(|a, b| a.1.cmp(b.1))
                .map(|(k, v)| (k, v))
                .unwrap();
            cells.push((
                *cellid,
                Cell {
                    feature_ids: Vec::new(),
                    nr_vertices: *nr_vertices,
                },
            ));
        } else {
            for (cellid, nr_vertices) in cell_vtx_cnt.iter() {
                cells.push((
                    *cellid,
                    Cell {
                        feature_ids: Vec::new(),
                        nr_vertices: *nr_vertices,
                    },
                ));
            }
        }
        Some(FeatureInGridCells { feature, cells })
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

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Crs(String);

impl Crs {
    fn from_model(model: &cityjson::v2_0::OwnedCityModel) -> Result<Self, Box<dyn std::error::Error>> {
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
    pub path_jsonl: PathBuf,
    pub bbox: Bbox,
}

impl Feature {
    pub fn centroid(&self) -> [f64; 2] {
        self.centroid
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
    let mut bbox: Option<Bbox> = None;
    let mut x_sum = 0.0;
    let mut y_sum = 0.0;
    let mut nr_vertices = 0usize;
    let mut selected_object_types = Vec::new();
    let mut ignored_object_types = Vec::new();

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
                let Some(boundaries) = geometry.boundaries() else {
                    continue;
                };
                for vertex_ref in boundaries.vertices() {
                    let Some(vertex) = model.vertices().get(*vertex_ref) else {
                        continue;
                    };
                    let coordinate = vertex.to_array();
                    update_bbox(&mut bbox, coordinate);
                    x_sum += coordinate[0];
                    y_sum += coordinate[1];
                    nr_vertices += 1;
                }
            }
        } else if !ignored_object_types.contains(&object_type) {
            ignored_object_types.push(object_type);
        }
    }

    SelectedGeometryStats {
        bbox,
        centroid: (nr_vertices > 0).then_some([x_sum / nr_vertices as f64, y_sum / nr_vertices as f64]),
        nr_vertices,
        selected_object_types,
        ignored_object_types,
    }
}

fn count_vertices_in_grid(
    model: &cityjson::v2_0::OwnedCityModel,
    filter: Option<&Vec<CityObjectType>>,
    grid: &crate::spatial_structs::SquareGrid,
    bbox: &Bbox,
) -> HashMap<CellId, usize> {
    let mut cell_vtx_cnt: HashMap<CellId, usize> = HashMap::new();

    for (_id, cityobject) in model.cityobjects().iter() {
        let Some(object_type) = map_cityobject_type(cityobject.type_cityobject()) else {
            continue;
        };
        if !is_selected_type(filter, object_type) {
            continue;
        }
        let geometry_handles = cityobject.geometry().unwrap_or(&[]);
        for geometry_handle in geometry_handles {
            let Some(geometry) = model.get_geometry(*geometry_handle) else {
                continue;
            };
            let Some(boundaries) = geometry.boundaries() else {
                continue;
            };
            for vertex_ref in boundaries.vertices() {
                let Some(vertex) = model.vertices().get(*vertex_ref) else {
                    continue;
                };
                let point = [vertex.x(), vertex.y()];
                let cellid = grid.locate_point(&point);
                *cell_vtx_cnt.entry(cellid).or_insert(0) += 1;
            }
        }
    }

    for cellid in grid.intersect_bbox(bbox) {
        *cell_vtx_cnt.entry(cellid).or_insert(0) += 1;
    }

    cell_vtx_cnt
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
        None => *bbox = Some([coordinate[0], coordinate[1], coordinate[2], coordinate[0], coordinate[1], coordinate[2]]),
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

    #[test]
    fn test_crs_to_epsg() {
        let crs = Crs("https://www.opengis.net/def/crs/EPSG/0/7415".to_string());
        let epsg_code = crs.to_epsg().unwrap();
        assert_eq!(7415_u16, epsg_code);
    }

    #[test]
    fn test_feature_file_loads_through_cjlib() {
        let pb = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("data")
            .join("3dbag_feature_x71.city.jsonl");
        let base = std::fs::read(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("resources")
                .join("data")
                .join("3dbag_x00.city.json"),
        )
        .unwrap();
        let model = json::from_feature_file_with_base(pb, &base).unwrap();
        let stats = selected_geometry_stats(model.as_inner(), Some(&vec![CityObjectType::Building]));
        assert!(stats.bbox.is_some());
        assert!(stats.nr_vertices > 0);
        let bbox = stats.bbox.unwrap();
        assert!(bbox[0] > 0.0);
        assert!(bbox[1] > 0.0);
    }
}
