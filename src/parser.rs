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
use std::fs::read_to_string;
use std::iter;
use std::path::{Path, PathBuf};

use log::{debug, error, info};
use postgres::fallible_iterator::FallibleIterator;
use postgres::types::FromSql;
use postgres::{Client, Error, NoTls};
use serde::Deserialize;
use serde_json::from_str;
use walkdir::WalkDir;

use crate::spatial_structs::{Bbox, BboxQc, SquareGrid};

pub struct WorldDb {
    pub uri: String,
    table: String,
    geometry_column: String,
    primary_key: String,
    pub grid: SquareGrid,
    pub features: FeatureSetDb,
}

impl WorldDb {
    pub fn new(
        uri: &str,
        table: &str,
        geometry_column: &str,
        primary_key: &str,
        cellsize: u16,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let mut client = Client::connect(uri, NoTls)?;
        let extent_rw = Self::extent_rw(&mut client, table, geometry_column)?;
        let epsg: u16 = Self::epsg(&mut client, table, geometry_column)?;
        let q = format!("SELECT count(*) FROM {table};");
        let res = client.query_one(&q, &[])?;
        let nr_features: i64 = res.get(0);
        let mut features: FeatureSetDb = Vec::with_capacity(nr_features as usize + 1);
        features.resize(nr_features as usize + 1, FeatureDb::default());

        let grid = SquareGrid::new(&extent_rw, cellsize, epsg, Some(10.0));
        Ok(Self {
            uri: uri.to_string(),
            table: table.to_string(),
            geometry_column: geometry_column.to_string(),
            primary_key: primary_key.to_string(),
            grid,
            features,
        })
    }

    fn extent_rw(
        client: &mut Client,
        table: &str,
        geometry_column: &str,
    ) -> Result<Bbox, Box<dyn std::error::Error>> {
        let q = format!(
            "SELECT st_xmin(extent) xmin
             , st_ymin(extent) ymin
             , st_zmin(extent) zmin
             , st_xmax(extent) xmax
             , st_ymax(extent) ymax
             , st_zmax(extent) zmax
        FROM (SELECT st_3dextent({geometry_column}) extent FROM {table}) AS e;"
        );
        let row = client.query_one(&q, &[])?;
        let b: Bbox = [
            row.get(0),
            row.get(1),
            row.get(2),
            row.get(3),
            row.get(4),
            row.get(5),
        ];
        Ok(b)
    }

    fn epsg(
        client: &mut Client,
        table: &str,
        geometry_column: &str,
    ) -> Result<u16, Box<dyn std::error::Error>> {
        let q = format!("SELECT st_srid({geometry_column}) FROM {table} LIMIT 1;");
        let row = client.query_one(&q, &[])?;
        let epsg: i32 = row.get(0);
        Ok(epsg as u16)
    }

    pub fn index_with_grid(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let transaction_size: i32 = 100;
        info!("Predicting runtime in grid cells");
        let mut fid: usize = 0;
        let mut client = Client::connect(&self.uri, NoTls)?;

        let q = format!(
            "SELECT {pk}, st_astext({geom}), st_area({geom}) FROM {tbl};",
            pk = &self.primary_key,
            geom = &self.geometry_column,
            tbl = &self.table
        );
        let params: Vec<String> = vec![];
        let mut res_it = client.query_raw(&q, params)?;
        while let Some(feature_row) = res_it.next()? {
            let mut cell_vtx_cnt: HashMap<crate::spatial_structs::CellId, usize> = HashMap::new();
            let fid_pk: i32 = feature_row.get(0);
            let feature = FeatureDb {
                primary_key: fid_pk,
            };
            let wkt: &str = feature_row.get(1);
            let area = feature_row.get(2);
            let vertices = parse_wkt_polygon(wkt);
            let centroid = Self::centroid(vertices);
            let cellid = self.grid.locate_point(&centroid);
            *cell_vtx_cnt.entry(cellid).or_insert(1) += 1;
            if !cell_vtx_cnt.is_empty() {
                self.features[fid] = feature;
                let (cellid, nr_vertices) = cell_vtx_cnt
                    .iter()
                    .max_by(|a, b| a.1.cmp(b.1))
                    .map(|(k, v)| (k, v))
                    .unwrap();
                let cell = self.grid.cell_mut(cellid);
                // Yes, in fact we don't store the number of vertices, but we store the
                // total predicted runtime for the cell. This just a hack so that we
                // don't need to change the data structure.
                let runtime = predict_run_time(area).ceil() as usize;
                cell.nr_vertices += runtime;
                if !cell.feature_ids.contains(&fid) {
                    cell.feature_ids.push(fid)
                }
                fid += 1;
            }
        }
        Ok(())
    }

    fn centroid(vertices: Vec<[f64; 2]>) -> [f64; 2] {
        let mut x_sum: f64 = 0.0;
        let mut y_sum: f64 = 0.0;
        for [x, y] in vertices.iter() {
            x_sum += *x;
            y_sum += *y;
        }
        [
            (x_sum / vertices.len() as f64),
            (y_sum / vertices.len() as f64),
        ]
    }

    pub fn export_grid(&self) -> std::io::Result<()> {
        self.grid.export_db()
    }
}

/// Predict the reconstruction time with geoflow for the given feature, from the area
/// of the footprint. Uses a simple linear model. The model parameters were determined
/// from a random sample of 3D BAG version v210908_fd2cee53.
pub fn predict_run_time(footprint_area: f64) -> f64 {
    let intercept: f64 = -347.033256;
    let coefficient: f64 = 6.222044;
    intercept + coefficient * footprint_area
}

pub fn parse_wkt_polygon(wkt: &str) -> Vec<[f64; 2]> {
    let mut res: Vec<[f64; 2]> = Vec::new();
    let wlower = wkt.to_lowercase();
    let coords_str: Vec<&str> = wlower
        .strip_prefix("polygon(")
        .unwrap()
        .strip_suffix(')')
        .unwrap()
        .split(',')
        .collect();
    for coord in coords_str {
        let c = coord.trim_matches(|c| c == '(' || c == ')');
        let cv: Vec<&str> = c.split(' ').collect();
        let x: f64 = cv[0].parse().unwrap();
        let y: f64 = cv[1].parse().unwrap();
        res.push([x, y]);
    }
    res
}

#[derive(Clone, Debug, Default)]
pub struct FeatureDb {
    pub primary_key: i32,
}
pub type FeatureSetDb = Vec<FeatureDb>;

mod tests_worlddb {
    use super::parse_wkt_polygon;
    use super::WorldDb;
    use crate::spatial_structs::QuadTree;
    use postgres::{Client, Error, NoTls};

    #[test]
    fn test_parse_wkt_polygon() {
        let wkt = "POLYGON((138019.09 456805.535,138018.832 456802.644,138018.376 456802.686,138018.071 456799.248,138019.684 456799.105,138019.397 456795.847,138018.282 456795.945,138017.819 456790.686,138016.563 456790.796,138016.481 456789.807,138016.244 456787.169,138018.609 456786.992,138018.602 456786.87,138028.334 456786.172,138029.43 456798.371,138029.32 456798.381,138030.339 456810.397,138030.449 456810.387,138031.126 456818.864,138031.016 456818.872,138031.937 456830.906,138032.047 456830.898,138033.001 456842.031,138032.1 456843.22,138018.425 456844.181,138018.364 456843.352,138015.454 456843.577,138015.514 456844.395,138000 456845.528,137976.701 456847.27,137976.769 456848.187,137974.655 456848.345,137974.587 456847.428,137905.446 456852.655,137904.177 456834.915,137903.721 456828.958,137902.747 456829.031,137902.571 456826.648,137903.539 456826.575,137901.901 456804.718,137901.851 456804.722,137901.158 456795.545,137906.121 456795.193,137906.117 456795.143,137910.327 456794.843,137910.332 456794.913,137913.624 456794.69,137920.358 456794.163,137920.375 456794.203,137932.105 456793.348,137932.193 456795.172,137934.521 456795.015,137960.807 456793.026,137960.824 456793.245,137973.585 456792.254,137973.592 456792.354,137979.735 456791.891,137979.71 456791.562,137988.448 456790.938,137988.471 456791.257,137994.625 456790.813,137994.616 456790.693,138000 456790.295,138007.631 456789.738,138007.594 456789.229,138010.608 456789.012,138010.663 456789.78,138012.578 456789.641,138012.901 456794.089,138007.949 456794.447,138008.418 456801.412,138001.94 456801.879,138002 456802.871,138000 456803.004,137995.976 456803.293,137995.945 456802.982,137992.083 456803.267,137992.035 456802.422,137989.355 456802.623,137989.521 456804.369,137985.081 456804.722,137984.885 456802.486,137980.579 456802.808,137980.633 456803.665,137976.911 456803.982,137976.881 456803.576,137974.137 456803.801,137974.204 456804.741,137970.742 456805.027,137970.675 456804.056,137969.85 456804.11,137969.91 456804.98,137968.754 456805.084,137968.687 456804.207,137967.37 456804.319,137967.701 456809.208,137961.919 456809.698,137961.642 456805.243,137961.6 456804.525,137957.376 456804.866,137957.339 456804.39,137957.003 456804.42,137957.069 456805.398,137954.916 456805.575,137954.879 456805.071,137948.274 456805.573,137948.236 456805.074,137941.842 456805.561,137942.079 456808.888,137935.55 456809.354,137935.21 456804.477,137932.892 456804.634,137933.031 456806.73,137929.849 456806.96,137929.861 456807.182,137921.109 456807.831,137921.715 456816.575,137926.073 456816.273,137926.054 456816.057,137936.14 456815.305,137936.57 456821.718,137943.658 456821.311,137943.766 456823.514,137941.22 456823.353,137935.014 456823.839,137934.875 456821.82,137926.968 456822.243,137921.563 456822.63,137921.733 456825.075,137916.721 456825.451,137916.885 456827.851,137915.64 456827.942,137915.84 456830.978,137914.239 456831.074,137914.452 456834.122,137913.226 456834.212,137913.529 456838.613,137915.143 456838.479,137915.27 456840.251,137913.651 456840.385,137914.107 456846.627,137918.013 456846.342,137917.587 456840.652,137930.004 456839.693,137929.846 456837.356,137930.85 456837.264,137931.033 456839.614,137933.072 456839.456,137932.889 456837.113,137936.026 456836.876,137935.935 456835.784,137936.767 456835.717,137937.029 456839.146,137939.125 456839.001,137938.986 456836.842,137942.111 456836.586,137942.25 456838.752,137954.844 456837.793,137954.674 456835.348,137955.525 456835.281,137955.714 456837.726,137957.912 456837.556,137957.742 456835.206,137961.743 456834.915,137961.913 456837.252,137969.205 456836.713,137968.84 456831.817,137973.229 456831.47,137973.198 456831.197,137974.43 456831.1,137974.813 456836.301,137979.869 456835.944,137979.705 456833.613,137981.82 456833.462,137981.984 456835.761,137992.146 456835.015,137992.068 456833.822,137995.224 456833.572,137995.309 456834.772,137998.593 456834.517,137998.618 456834.822,138000 456834.716,138000.083 456834.71,138000 456833.842,137999.825 456832.012,138000 456831.996,138002.046 456831.803,138002.177 456833.2,138007.731 456832.7,138007.642 456831.703,138011.028 456831.4,138011.611 456837.933,138021.969 456837.007,138021.727 456834.187,138022.842 456834.087,138022.569 456831.029,138020.248 456831.236,138019.977 456828.208,138021.172 456828.101,138020.004 456815.063,138020.166 456815.071,138020.024 456814.057,138021 456813.969,138020.752 456811.2,138018.202 456811.428,138017.884 456807.88,138017.688 456805.659,138019.09 456805.535))";
        let res = parse_wkt_polygon(wkt);
        println!("{:?}", res);
    }

    #[test]
    fn test_extent() -> Result<(), Box<dyn std::error::Error>> {
        let uri = "postgresql://db3dbag_user:db3dbag_1234@localhost:5560/baseregisters";
        let mut client = Client::connect(uri.as_ref(), NoTls)?;
        if let Ok(bbox) = WorldDb::extent_rw(&mut client, "top10nl.gebouw", "geometrie_vlak") {
            println!("{:?}", bbox);
        };
        client.close()?;
        Ok(())
    }

    #[test]
    fn test_epsg() -> Result<(), Box<dyn std::error::Error>> {
        let uri = "postgresql://db3dbag_user:db3dbag_1234@localhost:5560/baseregisters";
        let mut client = Client::connect(uri.as_ref(), NoTls)?;
        if let Ok(epsg) = WorldDb::epsg(&mut client, "top10nl.gebouw", "geometrie_vlak") {
            assert_eq!(epsg, 28992_u16);
        };
        client.close()?;
        Ok(())
    }

    #[test]
    fn test_new() -> Result<(), Box<dyn std::error::Error>> {
        let uri = "postgresql://db3dbag_user:db3dbag_1234@localhost:5560/baseregisters";
        let world = WorldDb::new(&uri, "top10nl.gebouw", "geometrie_vlak", "fid", 250)?;
        Ok(())
    }
}

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
pub struct World {
    pub cityobject_types: Option<Vec<CityObjectType>>,
    pub crs: Crs,
    pub features: FeatureSet,
    pub grid: crate::spatial_structs::SquareGrid,
    pub path_features_root: PathBuf,
    pub path_metadata: PathBuf,
    pub transform: Transform,
}

impl World {
    pub fn new<P: AsRef<Path>>(
        path_metadata: P,
        path_features_root: P,
        cellsize: u16,
        cityobject_types: Option<Vec<CityObjectType>>,
        arg_minz: Option<i32>,
        arg_maxz: Option<i32>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let path_features_root = path_features_root.as_ref().to_path_buf();
        let path_metadata = path_metadata.as_ref().to_path_buf();
        let cm = CityJSONMetadata::from_file(&path_metadata)?;
        let crs = cm.metadata.reference_system;
        let transform = cm.transform;

        // FIXME: if cityobject_types is None, then all cityobject are ignored, instead of included
        // Compute the extent of the features and the number of features.
        // We don't store the computed extent explicitly, because the grid contains that info.
        let (extent_qc, nr_features, cityobject_types_ignored) =
            Self::extent_qc(&path_features_root, cityobject_types.as_ref());
        info!(
            "Found {} features of type {:?}",
            nr_features, &cityobject_types
        );
        info!("Ignored feature types: {:?}", &cityobject_types_ignored);
        debug!("extent_qc: {:?}", &extent_qc);
        let extent_rw = extent_qc.to_bbox(&transform, arg_minz, arg_maxz);
        debug!(
            "Computed extent from features in real-world coordinates: {:?}",
            &extent_rw
        );

        // Allocate the grid, but at this point it is still empty
        let epsg = crs.to_epsg()?;
        let grid = SquareGrid::new(&extent_rw, cellsize, epsg, Some(10.0));
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

    /// Compute the extent (in quantized coordinates), the number of features and the
    /// CityObject types that are present in the data but not selected.
    fn extent_qc<P: AsRef<Path>>(
        path_features: P,
        cityobject_types: Option<&Vec<CityObjectType>>,
    ) -> (BboxQc, usize, Vec<CityObjectType>) {
        info!(
            "Computing extent from the features of type {:?}",
            cityobject_types
        );
        // Do a first loop over the features to calculate their extent and their number.
        // Need a mutable iterator, because .next() consumes the next value and advances the iterator.
        let mut features_enum_iter = WalkDir::new(&path_features)
            .into_iter()
            .filter_map(Self::jsonl_path);
        // Init the extent with from the first feature of the requested types
        let mut extent_qc = BboxQc([0, 0, 0, 0, 0, 0]);
        let mut found_feature_type = false;
        let mut nr_features = 0;
        let mut cotypes_ignored: Vec<CityObjectType> = Vec::new();
        debug!("Searching for the first feature of the requested type...");
        loop {
            if let Some(feature_path) = features_enum_iter.next() {
                if let Ok(cf) = CityJSONFeatureVertices::from_file(&feature_path) {
                    if let Some(eqc) = cf.bbox_of_types(cityobject_types) {
                        extent_qc = eqc;
                        found_feature_type = true;
                        nr_features += 1;
                        break;
                    } else {
                        for (_, co) in cf.cityobjects.iter() {
                            if !cotypes_ignored.contains(&co.cotype) {
                                cotypes_ignored.push(co.cotype);
                            }
                        }
                    }
                } else {
                    error!("Failed to parse {:?}", &feature_path)
                }
            }
        }
        if !found_feature_type {
            panic!(
                "Did not find any CityJSONFeature of type {:?}",
                &cityobject_types
            );
        }
        debug!("First feature found. Iterating over all features to compute the extent.");
        for feature_path in features_enum_iter {
            if let Ok(cf) = CityJSONFeatureVertices::from_file(&feature_path) {
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
                    nr_features += 1;
                } else {
                    for (_, co) in cf.cityobjects.iter() {
                        if !cotypes_ignored.contains(&co.cotype) {
                            cotypes_ignored.push(co.cotype);
                        }
                    }
                }
            } else {
                error!("Failed to parse {:?}", &feature_path);
            }
        }
        (extent_qc, nr_features, cotypes_ignored)
    }

    /// Return the file path if the 'DirEntry' is a .jsonl file (eg. .city.jsonl).
    pub fn jsonl_path(walkdir_res: Result<walkdir::DirEntry, walkdir::Error>) -> Option<PathBuf> {
        if let Ok(entry) = walkdir_res {
            if let Some(ext) = entry.path().extension() {
                if ext == "jsonl" {
                    Some(entry.path().to_path_buf())
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            // TODO: notify the user if some path cannot be accessed (eg. permission), https://docs.rs/walkdir/latest/walkdir/struct.Error.html
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

    // Export the grid of the World into the working directory.
    pub fn export_grid(&self) -> std::io::Result<()> {
        self.grid
            .export(Some(&self.features), Some(&self.transform))
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

#[derive(Deserialize, Debug)]
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
#[derive(Deserialize, Debug)]
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
    /// FIXME: u16 is not enough for EPSG codes, because we can have values in 90k+
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
#[derive(Debug, Default, Clone, Ord, PartialOrd, Eq, PartialEq)]
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

#[derive(Debug, Deserialize, clap::ValueEnum, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
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
