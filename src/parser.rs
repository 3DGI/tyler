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
    vertices: Vec<[i64; 3]>,
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
    fn centroid_quantized(&self) -> [i32; 2] {
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
            (x_sum / self.vertices.len() as i64) as i32,
            (y_sum / self.vertices.len() as i64) as i32,
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
        let [mut x_min, mut y_min, mut z_min] = self.vertices[0].clone();
        let [mut x_max, mut y_max, mut z_max] = self.vertices[0].clone();
        for [x, y, z] in self.vertices.iter() {
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
        [x_min, y_min, z_min, x_max, y_max, z_max]
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
        Ok(Feature {
            centroid_quantized: cf.centroid_quantized(),
            nr_vertices: cf.vertex_count(),
            path_jsonl: path.as_ref().to_path_buf(),
        })
    }
}

/// Stores the information that is computed from a CityJSONFeature.
pub struct Feature {
    centroid_quantized: [i32; 2],
    nr_vertices: u16,
    pub path_jsonl: PathBuf,
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
        let cityjsonfeature_str = r#"{
            "type": "CityJSONFeature",
            "id": "id-1",
            "CityObjects": {},
            "vertices": [[1,1,1], [2,2,2], [3,3,3]]
        }"#;
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
