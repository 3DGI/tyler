//! CityJSON parser.
//! The module is responsible for parsing CityJSON data and populating the World.
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
use std::fs::{read_to_string, File};
use std::path::{Path, PathBuf};

use log::{debug, error, info, warn};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::from_str;
use walkdir::WalkDir;

use crate::spatial_structs::BboxQc;

/// Represents the "world" that contains some features and needs to be partitioned into
/// tiles.
///
/// # Members
///
/// `path_features_root` - The path to the root directory containing all features.
///
/// `path_metadata` - The path to the JSON file that stores the
/// [CityJSON object](https://www.cityjson.org/specs/1.1.3/#cityjson-object)
/// (also called CityJSON metadata in *tyler*).
///
/// `cityobject_types` - The World only contains features of these types.
#[derive(Serialize, Deserialize)]
pub struct World {
    pub cityobject_types: Option<Vec<CityObjectType>>,
    pub crs: Crs,
    pub features: FeatureSet,
    pub grid: crate::spatial_structs::SquareGrid,
    pub path_features_root: PathBuf,
    pub path_metadata: PathBuf,
    pub transform: Transform,
}

struct ExtentQcResult {
    extent_qc: BboxQc,
    nr_features: usize,
    cityobject_types_ignored: Vec<CityObjectType>,
    nr_features_ignored: usize,
}

struct FeatureDirsFiles {
    feature_dirs: Vec<PathBuf>,
    feature_files: Vec<PathBuf>,
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
        let cm = CityJSONMetadata::from_file(&path_metadata)?;
        let crs = cm.metadata.reference_system;
        let transform = cm.transform;

        info!(
            "Computing extent from the features of type {:?}",
            cityobject_types
        );
        // FIXME: if cityobject_types is None, then all cityobject are ignored, instead of included
        // Compute the extent of the features and the number of features.
        // We don't store the computed extent explicitly, because the grid contains that info.
        let feature_dirs_files = Self::find_feature_dirs_and_files(&path_features_root);
        // Walk the subdirectories of the root
        debug!(
            "Found {} subdirectories and {} CityJSONFeature files at the root directory",
            feature_dirs_files.feature_dirs.len(),
            feature_dirs_files.feature_files.len()
        );
        let extents: Vec<ExtentQcResult> = feature_dirs_files
            .feature_dirs
            .into_par_iter()
            .filter_map(|dir| Self::extent_qc(dir, cityobject_types.as_ref()))
            .collect();
        let mut nr_features = 0;
        let mut nr_features_ignored = 0;
        let mut extent_qc: BboxQc = BboxQc::default();
        let mut cityobject_types_ignored: Vec<CityObjectType> = Vec::new();
        for (i, extent) in extents.iter().enumerate() {
            nr_features += extent.nr_features;
            nr_features_ignored += extent.nr_features_ignored;
            if i == 0 {
                extent_qc = extent.extent_qc.clone();
            } else {
                extent_qc.update_with(&extent.extent_qc);
            }
            for cotype in &extent.cityobject_types_ignored {
                if !cityobject_types_ignored.contains(cotype) {
                    cityobject_types_ignored.push(*cotype);
                }
            }
        }
        // Walk the files at the root and update the counters
        for feature_path in &feature_dirs_files.feature_files {
            Self::extent_qc_file(
                cityobject_types.as_ref(),
                &mut extent_qc,
                &mut nr_features,
                &mut nr_features_ignored,
                &mut cityobject_types_ignored,
                feature_path,
            );
        }
        if nr_features == 0 {
            panic!(
                "Did not find any CityJSONFeatures of type {:?}",
                cityobject_types
            );
        }
        info!(
            "Found {} features of type {:?}",
            nr_features, &cityobject_types
        );
        info!(
            "Ignored {} features of type {:?}",
            nr_features_ignored, &cityobject_types_ignored
        );
        debug!("extent_qc: {:?}", &extent_qc);
        let extent_rw = extent_qc.to_bbox(&transform, arg_minz, arg_maxz);
        info!(
            "Computed extent from features: {}",
            crate::spatial_structs::bbox_to_wkt(&extent_rw)
        );

        // Allocate the grid, but at this point it is still empty
        let epsg = crs.to_epsg()?;
        let grid = crate::spatial_structs::SquareGrid::new(&extent_rw, cellsize, epsg, Some(10.0));
        debug!("{}", grid);

        // Allocate the features container, but at this point it is still empty
        let mut features: FeatureSet = Vec::with_capacity(nr_features + 1);
        features.resize(nr_features + 1, Feature::default());

        Ok(Self {
            features,
            crs,
            transform,
            grid,
            cityobject_types,
            path_features_root,
            path_metadata,
        })
    }

    /// Find the direct subdirectories and CityJSONFeature files in the directory.
    /// Returns a Vec of the subdirectory paths and a Vec of CityJSONFeature paths.
    /// Note that it is not guaranteed that the returned directories contain any CityJSONFeatures.
    fn find_feature_dirs_and_files(path_features_root: &PathBuf) -> FeatureDirsFiles {
        let mut path_features_root_dirs: Vec<PathBuf> = Vec::new();
        let mut path_features_root_files: Vec<PathBuf> = Vec::new();
        for entry_res in WalkDir::new(&path_features_root).min_depth(1).max_depth(1) {
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

    /// Compute the extent (in quantized coordinates), the number of features and the
    /// CityObject types that are present in the data but not selected.
    fn extent_qc<P: AsRef<Path> + std::fmt::Debug>(
        path_features: P,
        cityobject_types: Option<&Vec<CityObjectType>>,
    ) -> Option<ExtentQcResult> {
        // Do a first loop over the features to calculate their extent and their number.
        // Need a mutable iterator, because .next() consumes the next value and advances the iterator.
        let mut features_enum_iter = WalkDir::new(&path_features)
            .into_iter()
            .filter_map(Self::jsonl_path);
        // Init the extent with from the first feature of the requested types
        let mut extent_qc = BboxQc([0, 0, 0, 0, 0, 0]);
        let mut found_feature_type = false;
        let mut nr_features = 0;
        let mut nr_features_ignored = 0;
        let mut cityobject_types_ignored: Vec<CityObjectType> = Vec::new();
        // Iterate only until the first feature is found
        #[allow(clippy::while_let_on_iterator)]
        while let Some(feature_path) = features_enum_iter.next() {
            if let Ok(cf) = CityJSONFeatureVertices::from_file(&feature_path) {
                if let Some(eqc) = cf.bbox_of_types(cityobject_types) {
                    extent_qc = eqc;
                    found_feature_type = true;
                    nr_features += 1;
                    break;
                } else {
                    for (_, co) in cf.cityobjects.iter() {
                        if !cityobject_types_ignored.contains(&co.cotype) {
                            cityobject_types_ignored.push(co.cotype);
                        }
                        nr_features_ignored += 1;
                    }
                }
            } else {
                warn!("Failed to parse {:?}", &feature_path)
            }
        }
        if !found_feature_type {
            return None;
        }
        for feature_path in features_enum_iter {
            Self::extent_qc_file(
                cityobject_types,
                &mut extent_qc,
                &mut nr_features,
                &mut nr_features_ignored,
                &mut cityobject_types_ignored,
                &feature_path,
            );
        }
        Some(ExtentQcResult {
            extent_qc,
            nr_features,
            cityobject_types_ignored,
            nr_features_ignored,
        })
    }

    fn extent_qc_file(
        cityobject_types: Option<&Vec<CityObjectType>>,
        extent_qc: &mut BboxQc,
        nr_features: &mut usize,
        nr_features_ignored: &mut usize,
        cityobject_types_ignored: &mut Vec<CityObjectType>,
        feature_path: &PathBuf,
    ) {
        if let Ok(cf) = CityJSONFeatureVertices::from_file(feature_path) {
            if let Some(bbox_qc) = cf.bbox_of_types(cityobject_types) {
                let [x_min, y_min, z_min, x_max, y_max, z_max] = bbox_qc.0;
                if x_min < extent_qc.0[0] {
                    extent_qc.0[0] = x_min
                } else if x_max > extent_qc.0[3] {
                    extent_qc.0[3] = x_max
                }
                if y_min < extent_qc.0[1] {
                    extent_qc.0[1] = y_min
                } else if y_max > extent_qc.0[4] {
                    extent_qc.0[4] = y_max
                }
                if z_min < extent_qc.0[2] {
                    extent_qc.0[2] = z_min
                } else if z_max > extent_qc.0[5] {
                    extent_qc.0[5] = z_max
                }
                *nr_features += 1;
            } else {
                for (_, co) in cf.cityobjects.iter() {
                    if !cityobject_types_ignored.contains(&co.cotype) {
                        cityobject_types_ignored.push(co.cotype);
                    }
                    *nr_features_ignored += 1;
                }
            }
        } else {
            error!("Failed to parse {:?}", &feature_path);
        }
    }

    /// Return the file path if the 'DirEntry' is a .jsonl file (eg. .city.jsonl).
    pub fn jsonl_path(walkdir_res: Result<walkdir::DirEntry, walkdir::Error>) -> Option<PathBuf> {
        if let Ok(entry) = walkdir_res {
            Self::direntry_to_jsonl(entry)
        } else {
            // TODO: notify the user if some path cannot be accessed (eg. permission), https://docs.rs/walkdir/latest/walkdir/struct.Error.html
            None
        }
    }

    /// Convert a [walkdir::DirEntry] to [PathBuf] if the file is CityJSONFeature file.
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

    // Loop through the features and assign the features to the grid cells.
    pub fn index_with_grid(&mut self) {
        let feature_set_paths_iter = WalkDir::new(&self.path_features_root)
            .into_iter()
            .filter_map(Self::jsonl_path)
            .enumerate();
        // For each feature_path (parallel) -- but we would need to mutate a variable from a parallel loop, creating a data race condition, we'll fix this later
        //  parse the feature
        //  for each vertex of the feature
        //      cellid <- locate vertex in grid
        //      cell <- get mutable cell reference from cellid
        //      increment vertex count in cell
        //      add feature id to cell
        info!("Counting vertices in grid cells");
        let mut fid: usize = 0;
        for (_, feature_path) in feature_set_paths_iter {
            let cf = CityJSONFeatureVertices::from_file(&feature_path);
            if let Ok(featurevertices) = cf {
                // We make a (cellid, vertex count) map and assign the feature to the cell that
                // contains the most of the feature's vertices.
                // But maybe a HashMap is not the most performant solution here? A Vec of tuples?
                let mut cell_vtx_cnt: HashMap<crate::spatial_structs::CellId, usize> =
                    HashMap::new();
                for (_, co) in featurevertices.cityobjects.iter() {
                    // If the object_type argument was not passed, that means that we need all
                    // CityObject types. If it was passed, then we filter with its values.
                    // Doing this condition-tree would be much simpler if Option.is_some_and()
                    // was stable feature already.
                    let mut do_compute = self.cityobject_types.is_none();
                    if let Some(ref cotypes) = self.cityobject_types {
                        do_compute = cotypes.contains(&co.cotype);
                    }
                    if do_compute {
                        // Just counting vertices here
                        for vtx_qc in featurevertices.vertices.iter() {
                            let vtx_rw = [
                                (vtx_qc[0] as f64 * self.transform.scale[0])
                                    + self.transform.translate[0],
                                (vtx_qc[1] as f64 * self.transform.scale[1])
                                    + self.transform.translate[1],
                            ];
                            let cellid = self.grid.locate_point(&vtx_rw);
                            *cell_vtx_cnt.entry(cellid).or_insert(1) += 1;
                        }
                    }
                }
                // After counting the object vertices in the cells, we need to
                // assign the object to the cells that intersect with its bbox,
                // because of https://github.com/3DGI/tyler/issues/28
                if let Some(bbox_qc) = featurevertices.bbox_of_types(self.cityobject_types.as_ref())
                {
                    let bbox = bbox_qc.to_bbox(&self.transform, None, None);
                    let intersecting_cellids = self.grid.intersect_bbox(&bbox);
                    for cellid in intersecting_cellids {
                        // Just add a new entry with the intersecting cell to the map, but no not
                        // increase the vertex count, because the vertices have been counted
                        // already, these might be cells where the object does not actually have a
                        // vertex.
                        // REVIEW: actually, let's just increase the vertex count
                        *cell_vtx_cnt.entry(cellid).or_insert(1) += 1;
                    }
                }

                if !cell_vtx_cnt.is_empty() {
                    // We found at least one CityObject of the required type
                    self.features[fid] = featurevertices.to_feature(&feature_path);
                    // TODO: what other cityobject types need to have 1-1 cell assignment?
                    if let Some(ref cotypes) = self.cityobject_types {
                        if cotypes.contains(&CityObjectType::Building)
                            || cotypes.contains(&CityObjectType::BuildingPart)
                        {
                            // In this case we have a 1-1 feature-to-cell assignment, we only retain the vertex
                            // count in the cell that gets the feature.
                            // The cell that receives the feature is the one with the highest vertex count
                            // of the feature.
                            // However, with this method it is not possible to combine cityobject types that
                            // require different cell-assignment methods into the same tileset.
                            // E.g. terrain features need to be duplicated across cells, buildings need to
                            // unique. The tileset for them must be generated separately.
                            let (cellid, nr_vertices) = cell_vtx_cnt
                                .iter()
                                .max_by(|a, b| a.1.cmp(b.1))
                                .map(|(k, v)| (k, v))
                                .unwrap();
                            let cell = self.grid.cell_mut(cellid);
                            cell.nr_vertices += nr_vertices;
                            if !cell.feature_ids.contains(&fid) {
                                cell.feature_ids.push(fid)
                            }
                        } else {
                            for (cellid, nr_vertices) in cell_vtx_cnt.iter() {
                                let cell = self.grid.cell_mut(cellid);
                                cell.nr_vertices += nr_vertices;
                                if !cell.feature_ids.contains(&fid) {
                                    cell.feature_ids.push(fid)
                                }
                            }
                        }
                        fid += 1;
                    }
                }
            } else {
                error!("Failed to parse the feature {:?}", &feature_path);
            }
        }
    }

    /// Export the grid of the World into the working directory.
    pub fn export_grid(
        &self,
        export_features: bool,
        output_dir: Option<&Path>,
    ) -> std::io::Result<()> {
        if export_features {
            self.grid
                .export(Some(&self.features), Some(&self.transform), output_dir)
        } else {
            self.grid.export(None, None, output_dir)
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

/// A partial [CityJSON object](https://www.cityjson.org/specs/1.1.3/#cityjson-object).
/// It is partial, because we only store the metadata that is necessary for parsing the
/// CityJSONFeatures.
#[derive(Deserialize, Debug)]
pub struct CityJSONMetadata {
    pub transform: Transform,
    pub metadata: Metadata,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Transform {
    pub scale: [f64; 3],
    pub translate: [f64; 3],
}

#[derive(Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct Metadata {
    pub reference_system: Crs,
}

/// Coordinate Reference System as defined by the
/// [referenceSystem](https://www.cityjson.org/specs/1.1.3/#referencesystem-crs) CityJSON object.
#[derive(Serialize, Deserialize, Debug)]
pub struct Crs(String);

impl Crs {
    /// Return the EPSG code from the CRS definition, if the CRS definition is indeed an EPSG.
    ///
    /// ## Examples
    /// ```
    /// let crs = CRS("https://www.opengis.net/def/crs/EPSG/0/7415");
    /// let epsg_code = crs.to_epsg().unwrap();
    /// assert_eq!(7415_u16, epsg_code);
    /// ```
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
        return if let Some(c) = parts.last() {
            let code: u16 = c.parse::<u16>().unwrap();
            Ok(code)
        } else {
            Err(Box::try_from(format!(
                "the CRS definition should contain the EPSG code as its last element: {}",
                self.0
            ))
            .unwrap())
        };
    }
}

/// Container for storing the CityJSONFeature vertices.
///
/// CityJSONFeature coordinates are supposed to be within the range of an `i32`,
/// `[-2147483648, 2147483647]`.
/// It allocates for the vertex container. I tried zero-copy (zero-allocation) deserialization
/// from the JSON string with the [zerovec](https://crates.io/crates/zerovec) crate
/// (see [video](https://youtu.be/DM2DI3ZI_BQ) for details), but I was getting an error of
/// "Attempted to build VarZeroVec out of elements that cumulatively are larger than a u32 in size"
/// from the zerovec crate, and I didn't investigate further.
#[derive(Deserialize, Debug)]
pub struct CityJSONFeatureVertices {
    #[serde(rename = "CityObjects")]
    pub cityobjects: HashMap<String, CityObject>,
    pub vertices: Vec<[i64; 3]>,
}

impl CityJSONMetadata {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let cm_str = read_to_string(path.as_ref())?;
        let cm: CityJSONMetadata = from_str(&cm_str)?;
        Ok(cm)
    }
}

impl CityJSONFeatureVertices {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let cf_str = read_to_string(path.as_ref())?;
        let cf: CityJSONFeatureVertices = from_str(&cf_str)?;
        Ok(cf)
    }

    /// Return the number of vertices of the feature.
    /// We assume that the number of vertices in a feature does not exceed 65535 (thus `u16`).
    fn vertex_count(&self) -> u16 {
        self.vertices.len() as u16
    }

    /// Feature centroid (2D) computed as the average coordinate.
    /// The centroid coordinates are quantized, so they need to be transformed back to real-world
    /// coordinates.
    /// It is more efficient to apply the transformation once, when the centroid is computed, than
    /// applying it to each vertex in the loop of computing the average coordinate.
    fn centroid_qc(&self) -> [i64; 2] {
        let mut x_sum: i64 = 0;
        let mut y_sum: i64 = 0;
        for [x, y, _z] in self.vertices.iter() {
            x_sum += *x;
            y_sum += *y;
        }
        // Yes, we divide an integer with an integer and we discard the decimals, but that's ok,
        // because the quantized coordinates (integers) already include the decimals of the
        // real-world coordinates. Thus, when the quantized centroid is scaled to the real-world
        // coordinate with a factor `< 0` (eg. 0.001), we will get accurate-enough coordinates
        // for the centroid.
        [
            (x_sum / self.vertices.len() as i64),
            (y_sum / self.vertices.len() as i64),
        ]
    }

    /// Feature centroid (2D) computed as the average coordinate.
    /// The centroid coordinates are real-world coordinates (thus they are transformed back to
    /// real-world coordinates from the quantized coordinates).
    #[allow(dead_code)]
    fn centroid(&self, transform: &Transform) -> [f64; 2] {
        let [ctr_x, ctr_y] = self.centroid_qc();
        [
            (ctr_x as f64 * transform.scale[0]) + transform.translate[0],
            (ctr_y as f64 * transform.scale[1]) + transform.translate[1],
        ]
    }

    /// Compute the 3D bounding box of the feature.
    /// Returns quantized coordinates.
    #[allow(dead_code)]
    pub fn bbox_qc(&self) -> crate::spatial_structs::BboxQc {
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0];
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0];
        for [x, y, z] in self.vertices.iter() {
            if *x < x_min {
                x_min = *x
            } else if *x > x_max {
                x_max = *x
            }
            if *y < y_min {
                y_min = *y
            } else if *y > y_max {
                y_max = *y
            }
            if *z < z_min {
                z_min = *z
            } else if *z > z_max {
                z_max = *z
            }
        }
        BboxQc([x_min, y_min, z_min, x_max, y_max, z_max])
    }

    /// Compute the 3D bounding box of only the provided CityObject types in the feature.
    /// Returns quantized coordinates.
    pub fn bbox_of_types(&self, cityobject_types: Option<&Vec<CityObjectType>>) -> Option<BboxQc> {
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0];
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0];
        let mut found_co_geometry = false;
        for (_, co) in self.cityobjects.iter() {
            // If the object_type argument was not passed, that means that we need all
            // CityObject types. If it was passed, then we filter with its values.
            // Doing this condition-tree would be much simpler if Option.is_some_and()
            // was stable feature already.
            let mut do_compute = cityobject_types.is_none();
            if let Some(cotypes) = cityobject_types {
                do_compute = cotypes.contains(&co.cotype);
            }
            if do_compute {
                for geom in co.geometry.iter() {
                    match geom {
                        Geometry::MultiSurface { boundaries, .. } => {
                            for srf in boundaries {
                                for ring in srf {
                                    for vtx in ring {
                                        let [x, y, z] = &self.vertices[*vtx];
                                        if *x < x_min {
                                            x_min = *x
                                        } else if *x > x_max {
                                            x_max = *x
                                        }
                                        if *y < y_min {
                                            y_min = *y
                                        } else if *y > y_max {
                                            y_max = *y
                                        }
                                        if *z < z_min {
                                            z_min = *z
                                        } else if *z > z_max {
                                            z_max = *z
                                        }
                                    }
                                }
                            }
                            found_co_geometry = true;
                        }
                        Geometry::Solid { boundaries, .. } => {
                            for shell in boundaries {
                                for srf in shell {
                                    for ring in srf {
                                        for vtx in ring {
                                            let [x, y, z] = &self.vertices[*vtx];
                                            if *x < x_min {
                                                x_min = *x
                                            } else if *x > x_max {
                                                x_max = *x
                                            }
                                            if *y < y_min {
                                                y_min = *y
                                            } else if *y > y_max {
                                                y_max = *y
                                            }
                                            if *z < z_min {
                                                z_min = *z
                                            } else if *z > z_max {
                                                z_max = *z
                                            }
                                        }
                                    }
                                }
                            }
                            found_co_geometry = true;
                        }
                    }
                }
            }
        }
        if found_co_geometry {
            Some(BboxQc([x_min, y_min, z_min, x_max, y_max, z_max]))
        } else {
            None
        }
    }

    /// Compute the 2D quantized centroid and the 3D bounding box in one loop.
    ///
    /// Combines the [centroid_quantized] and [bbox] methods to compute the values in a single
    /// loop over the vertices.
    fn centroid_bbox_qc(&self) -> [i64; 8] {
        let mut x_sum: i64 = 0;
        let mut y_sum: i64 = 0;
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0];
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0];
        for [x, y, z] in self.vertices.iter() {
            x_sum += x;
            y_sum += y;
            if *x < x_min {
                x_min = *x
            } else if *x > x_max {
                x_max = *x
            }
            if *y < y_min {
                y_min = *y
            } else if *y > y_max {
                y_max = *y
            }
            if *z < z_min {
                z_min = *z
            } else if *z > z_max {
                z_max = *z
            }
        }
        let x_ctr = x_sum / self.vertices.len() as i64;
        let y_ctr = y_sum / self.vertices.len() as i64;
        [x_ctr, y_ctr, x_min, y_min, z_min, x_max, y_max, z_max]
    }

    /// Sets the 'path_jsonl' to default.
    pub fn to_feature<P: AsRef<Path>>(&self, path: P) -> Feature {
        let ctr_bbox = self.centroid_bbox_qc();
        Feature {
            centroid_qc: [ctr_bbox[0], ctr_bbox[1]],
            nr_vertices: self.vertex_count(),
            path_jsonl: path.as_ref().to_path_buf(),
            bbox_qc: BboxQc([
                ctr_bbox[2],
                ctr_bbox[3],
                ctr_bbox[4],
                ctr_bbox[5],
                ctr_bbox[6],
                ctr_bbox[7],
            ]),
        }
    }
}

/// Stores the information that is computed from a CityJSONFeature.
#[derive(Debug, Default, Clone, Ord, PartialOrd, Eq, PartialEq, Serialize, Deserialize)]
pub struct Feature {
    pub(crate) centroid_qc: [i64; 2],
    pub(crate) nr_vertices: u16,
    pub path_jsonl: PathBuf,
    pub bbox_qc: BboxQc,
}

impl Feature {
    pub fn centroid(&self, transform: &Transform) -> [f64; 2] {
        let [ctr_x, ctr_y] = self.centroid_qc;
        [
            (ctr_x as f64 * transform.scale[0]) + transform.translate[0],
            (ctr_y as f64 * transform.scale[1]) + transform.translate[1],
        ]
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
    #[serde(rename = "+GenericCityObject")]
    GenericCityObject,
}

impl fmt::Display for CityObjectType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

// Indexed geometry
type Vertex = usize;
type Ring = Vec<Vertex>;
type Surface = Vec<Ring>;
type Shell = Vec<Surface>;
type MultiSurface = Vec<Surface>;
type Solid = Vec<Shell>;

#[derive(Deserialize, Debug)]
#[serde(tag = "type")]
enum Geometry {
    MultiSurface { boundaries: MultiSurface },
    Solid { boundaries: Solid },
}

#[derive(Deserialize, Debug)]
pub struct CityObject {
    #[serde(rename = "type")]
    pub cotype: CityObjectType,
    geometry: Vec<Geometry>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::from_str;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("data")
    }

    #[test]
    fn test_crs_to_epsg() {
        let crs = Crs("https://www.opengis.net/def/crs/EPSG/0/7415".to_string());
        let epsg_code = crs.to_epsg().unwrap();
        assert_eq!(7415_u16, epsg_code);
    }

    #[test]
    fn test_cityjsonmetadata() -> serde_json::Result<()> {
        let cityjson_str = r#"{
            "type": "CityJSON",
            "version": "1.1",
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [0.0, 0.0, 0.0]
            },
            "metadata": {
                "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/7415",
                "title": "MyTitle"
            },
            "CityObjects": {},
            "vertices": []
        }"#;
        let cm: CityJSONMetadata = from_str(cityjson_str)?;
        println!("{:#?}", cm.metadata.reference_system);
        println!("{:#?}, {:#?}", cm.transform.scale, cm.transform.translate);
        Ok(())
    }

    #[test]
    fn test_cityjsonfeaturevertices() -> serde_json::Result<()> {
        let cityjsonfeature_str = r#"{"type":"CityJSONFeature","CityObjects":{"b70a1e56f-debe-11e7-8ec4-89be260623ee":{"type":"Road","geometry":[{"type":"MultiSurface","lod":"1","boundaries":[[[0,1,2]],[[1,3,4]],[[1,0,3]],[[2,5,0]],[[2,6,5]],[[7,8,6]],[[9,10,11]],[[10,12,13]],[[14,15,16]],[[17,14,16]],[[18,19,20]],[[21,22,23]],[[24,25,26]],[[20,27,25]],[[20,19,27]],[[28,29,30]],[[9,23,10]],[[31,32,28]],[[31,33,32]],[[34,31,28]],[[35,34,28]],[[35,28,30]],[[36,22,37]],[[30,29,18]],[[36,38,39]],[[18,29,19]],[[40,26,41]],[[42,40,41]],[[24,20,25]],[[17,43,42]],[[40,42,43]],[[26,40,24]],[[43,17,16]],[[15,14,39]],[[39,38,15]],[[37,38,36]],[[21,37,22]],[[9,21,23]],[[11,10,44]],[[44,10,13]],[[13,12,7]],[[7,12,8]],[[2,7,6]],[[45,46,4]],[[46,1,4]],[[47,46,45]],[[48,47,45]]]}],"attributes":{"3df_id":"G0200.42b3d391aef50268e0530a0a28492340"}}},"vertices":[[23241731,-6740287,16980],[23243271,-6737886,17050],[23241947,-6737751,17030],[23243688,-6740239,16990],[23244961,-6739729,16990],[23241021,-6740116,16970],[23240334,-6739867,16960],[23240760,-6737152,17020],[23239680,-6739542,16950],[23207572,-6713437,17050],[23206398,-6715354,17010],[23211403,-6716175,17030],[23239066,-6739146,16950],[23224416,-6725473,17000],[23154567,-6713216,17160],[23200871,-6711570,17040],[23153430,-6710683,17210],[23152498,-6713168,17190],[23148683,-6700000,17400],[23145589,-6704251,17390],[23148683,-6706399,17330],[23205998,-6712657,17050],[23204080,-6714161,17010],[23205285,-6714668,17010],[23149399,-6707907,17300],[23146208,-6708310,17330],[23147259,-6710093,17300],[23145640,-6706320,17370],[23146890,-6619484,17810],[23145656,-6700000,17440],[23149034,-6662558,17710],[23140404,-6619323,17890],[23139266,-6623569,17820],[23139266,-6619957,17790],[23149466,-6614336,17770],[23149281,-6634334,17790],[23202811,-6713844,17010],[23204339,-6712080,17050],[23202621,-6711716,17040],[23201509,-6713726,17010],[23150482,-6709178,17270],[23148723,-6711555,17260],[23150508,-6712602,17220],[23151857,-6710125,17240],[23219449,-6721924,17010],[23246174,-6738901,16990],[23244554,-6737539,17080],[23245629,-6736755,17120],[23246913,-6738228,17040]],"id":"b70a1e56f-debe-11e7-8ec4-89be260623ee"}"#;
        let cf: CityJSONFeatureVertices = from_str(cityjsonfeature_str)?;
        for v in cf.vertices.iter() {
            println!("{:#?}", v.first());
        }
        Ok(())
    }

    #[test]
    fn test_centroid() -> serde_json::Result<()> {
        let pb: PathBuf = test_data_dir().join("3dbag_feature_x71.city.jsonl");
        let cf: CityJSONFeatureVertices = CityJSONFeatureVertices::from_file(&pb).unwrap();
        let ctr_quantized = cf.centroid_qc();
        println!("quantized centroid: {:#?}", ctr_quantized);

        let pb: PathBuf = test_data_dir().join("3dbag_x00.city.json");
        let cm: CityJSONMetadata = CityJSONMetadata::from_file(&pb).unwrap();

        let ctr_real_world: (f64, f64) = (
            (ctr_quantized[0] as f64 * cm.transform.scale[0]) + cm.transform.translate[0],
            (ctr_quantized[1] as f64 * cm.transform.scale[1]) + cm.transform.translate[1],
        );
        println!("real-world centroid: {:#?}", ctr_real_world);

        Ok(())
    }
}

pub type FeatureSet = Vec<Feature>;
