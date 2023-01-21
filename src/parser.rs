use std::collections::HashMap;
use std::fmt;
use std::fs::read_to_string;
use std::path::{Path, PathBuf};

use serde::Deserialize;
use serde_json::from_str;

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
    pub reference_system: CRS,
}

/// Coordinate Reference System as defined by the
/// [referenceSystem](https://www.cityjson.org/specs/1.1.3/#referencesystem-crs) CityJSON object.
#[derive(Deserialize, Debug)]
pub struct CRS(String);

impl CRS {
    /// Return the EPSG code from the CRS definition, if the CRS definition is indeed an EPSG.
    ///
    /// ## Examples
    /// ```
    /// let crs = CRS("https://www.opengis.net/def/crs/EPSG/0/7415");
    /// let epsg_code = crs.to_epsg().unwrap();
    /// assert_eq!(7415_u16, epsg_code);
    /// ```
    pub fn to_epsg(&self) -> Result<u16, Box<dyn std::error::Error>> {
        let parts: Vec<&str> = self.0.split("/").collect();
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
    fn centroid_quantized(&self) -> [i64; 2] {
        let mut x_sum: i64 = 0;
        let mut y_sum: i64 = 0;
        for [x, y, _z] in self.vertices.iter() {
            x_sum = x_sum + *x as i64;
            y_sum = y_sum + *y as i64;
        }
        // Yes, we divide an integer with an integer and we discard the decimals, but that's ok,
        // because the quantized coordinates (integers) already include the decimals of the
        // real-world coordinates. Thus, when the quantized centroid is scaled to the real-world
        // coordinate with a factor `< 0` (eg. 0.001), we will get accurate-enough coordinates
        // for the centroid.
        [
            (x_sum / self.vertices.len() as i64) as i64,
            (y_sum / self.vertices.len() as i64) as i64,
        ]
    }

    /// Feature centroid (2D) computed as the average coordinate.
    /// The centroid coordinates are real-world coordinates (thus they are transformed back to
    /// real-world coordinates from the quantized coordinates).
    fn centroid(&self, transform: &Transform) -> [f64; 2] {
        let [ctr_x, ctr_y] = self.centroid_quantized();
        [
            (ctr_x as f64 * transform.scale[0]) + transform.translate[0],
            (ctr_y as f64 * transform.scale[1]) + transform.translate[1],
        ]
    }

    /// Compute the 3D bounding box of the feature.
    /// Returns quantized coordinates.
    pub fn bbox(&self) -> [i64; 6] {
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
        [x_min, y_min, z_min, x_max, y_max, z_max]
    }

    /// Compute the 3D bounding box of only the provided CityObject types in the feature.
    /// Returns quantized coordinates.
    pub fn bbox_of_types(&self, cotypes: &[&CityObjectType]) -> Option<[i64; 6]> {
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0];
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0];
        let mut found_co_geometry = false;
        for (_, co) in self.cityobjects.iter() {
            if cotypes.contains(&&co.cotype) {
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
            Some([x_min, y_min, z_min, x_max, y_max, z_max])
        } else {
            None
        }
    }

    /// Compute the 2D quantized centroid and the 3D bounding box in one loop.
    ///
    /// Combines the [centroid_quantized] and [bbox] methods to compute the values in a single
    /// loop over the vertices.
    fn centroid_quantized_bbox(&self) -> [i64; 8] {
        let mut x_sum: i64 = 0;
        let mut y_sum: i64 = 0;
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0].clone();
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0].clone();
        for [x, y, z] in self.vertices.iter() {
            x_sum = x_sum + x;
            y_sum = y_sum + y;
            if *x < x_min {
                x_min = x.clone()
            } else if *x > x_max {
                x_max = x.clone()
            }
            if *y < y_min {
                y_min = y.clone()
            } else if *y > y_max {
                y_max = y.clone()
            }
            if *z < z_min {
                z_min = z.clone()
            } else if *z > z_max {
                z_max = z.clone()
            }
        }
        let x_ctr = x_sum / self.vertices.len() as i64;
        let y_ctr = y_sum / self.vertices.len() as i64;
        [x_ctr, y_ctr, x_min, y_min, z_min, x_max, y_max, z_max]
    }

    /// Extracts some information from the CityJSONFeature and returns a tuple with them.
    pub fn file_to_tuple<P: AsRef<Path>>(path: P) -> Result<Feature, Box<dyn std::error::Error>> {
        let cf: CityJSONFeatureVertices = Self::from_file(path.as_ref())?;
        Ok(cf.to_feature(path.as_ref()))
    }

    /// Sets the 'path_jsonl' to default.
    pub fn to_feature<P: AsRef<Path>>(&self, path: P) -> Feature {
        let ctr_bbox = self.centroid_quantized_bbox();
        Feature {
            centroid_quantized: [ctr_bbox[0], ctr_bbox[1]],
            nr_vertices: self.vertex_count(),
            path_jsonl: path.as_ref().to_path_buf(),
            bbox_quantized: [
                ctr_bbox[2],
                ctr_bbox[3],
                ctr_bbox[4],
                ctr_bbox[5],
                ctr_bbox[6],
                ctr_bbox[7],
            ],
        }
    }
}

/// Stores the information that is computed from a CityJSONFeature.
#[derive(Debug, Default, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct Feature {
    pub(crate) centroid_quantized: [i64; 2],
    pub(crate) nr_vertices: u16,
    pub path_jsonl: PathBuf,
    pub bbox_quantized: [i64; 6],
}

impl Feature {
    pub fn centroid(&self, cm: &CityJSONMetadata) -> [f64; 2] {
        let [ctr_x, ctr_y] = self.centroid_quantized;
        [
            (ctr_x as f64 * cm.transform.scale[0]) + cm.transform.translate[0],
            (ctr_y as f64 * cm.transform.scale[1]) + cm.transform.translate[1],
        ]
    }
}

#[derive(Debug, Deserialize, clap::ValueEnum, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
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
    MultiSurface {
        lod: String,
        boundaries: MultiSurface,
    },
    Solid {
        lod: String,
        boundaries: Solid,
    },
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

    fn test_output_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("output")
    }

    #[test]
    fn test_crs_to_epsg() {
        let crs = CRS("https://www.opengis.net/def/crs/EPSG/0/7415".to_string());
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
            println!("{:#?}", v.get(0));
        }
        Ok(())
    }

    #[test]
    fn test_centroid() -> serde_json::Result<()> {
        let pb: PathBuf = test_data_dir().join("3dbag_feature_x71.city.jsonl");
        let cf: CityJSONFeatureVertices = CityJSONFeatureVertices::from_file(&pb).unwrap();
        let ctr_quantized = cf.centroid_quantized();
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
