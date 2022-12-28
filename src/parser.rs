use crate::quadtree;

use std::fs::{File, read_to_string};
use std::io::BufReader;
use std::path::{Path, PathBuf};

use serde::Deserialize;
use serde_json::{from_reader, from_str};

/// A partial [CityJSON object](https://www.cityjson.org/specs/1.1.3/#cityjson-object).
/// It is partial, because we only store the metadata that is necessary for parsing the
/// CityJSONFeatures.
#[derive(Deserialize, Debug)]
pub struct CityJSONMetadata {
    transform: Transform,
    metadata: Metadata,
}

#[derive(Deserialize, Debug)]
struct Transform {
    scale: [f64; 3],
    translate: [f64; 3],
}

#[derive(Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
struct Metadata {
    reference_system: String,
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
struct CityJSONFeatureVertices {
    vertices: Vec<[i32; 3]>,
}


impl CityJSONFeatureVertices {
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
    fn centroid_quantized(&self) -> (i32, i32) {
        let mut x_sum: i32 = 0;
        let mut y_sum: i32 = 0;
        for [x,y,z] in self.vertices.iter() {
            x_sum = x_sum + x;
            y_sum = y_sum + y;
        }
        // Yes, we divide an integer with an integer and we discard the decimals, but that's ok,
        // because the quantized coordinates (integers) already include the decimals of the
        // real-world coordinates. Thus, when the quantized centroid is scaled to the real-world
        // coordinate with a factor `< 0` (eg. 0.001), we will get accurate-enough coordinates
        // for the centroid.
        (x_sum / self.vertices.len() as i32, y_sum / self.vertices.len() as i32)
    }

    /// Feature centroid (2D) computed as the average coordinate.
    /// The centroid coordinates are real-world coordinates (thus they are transformed back to
    /// real-world coordinates from the quantized coordinates).
    fn centroid(&self, transform: &Transform) -> (f64, f64) {
        let ctr = self.centroid_quantized();
        ((ctr.0 as f64 * transform.scale[0]) + transform.translate[0],
         (ctr.1 as f64 * transform.scale[1]) + transform.translate[1])
    }
}

/// Stores the information that is computed from a CityJSONFeature.
/// (morton code of centroid, vertex count, path to the file)
pub type FeatureTuple = (u128, u16, PathBuf);


pub fn parse_metadata(pb: &PathBuf) -> Result<CityJSONMetadata, Box<dyn std::error::Error>> {
    let file = File::open(&pb)?;
    let reader = BufReader::new(&file);
    let cm: CityJSONMetadata = from_reader(reader)?;
    Ok(cm)
}

/// Extracts some information from the CityJSONFeature and returns a tuple with them.
/// The function is meant to be used when iterating over a directory of `.city.jsonl` files.
pub fn feature_to_tuple<P: AsRef<Path>>(path: P, cm: &CityJSONMetadata) -> Result<FeatureTuple, Box<dyn std::error::Error>> {
    let cf_str = read_to_string(path.as_ref())?;
    let cf: CityJSONFeatureVertices = from_str(&cf_str)?;
    let center = cf.centroid(&cm.transform);
    Ok((quadtree::morton_encode(&center.0, &center.1), cf.vertex_count(), path.as_ref().to_path_buf()))
}


#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::from_str;
    use std::fs::read_to_string;

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
        let cf_str = read_to_string(&pb).unwrap();
        let cf: CityJSONFeatureVertices = from_str(&cf_str)?;
        let ctr_quantized = cf.centroid_quantized();
        println!("quantized centroid: {:#?}", ctr_quantized);

        let pb: PathBuf = test_data_dir().join("3dbag_x00.city.json");
        let cm_str = read_to_string(&pb).unwrap();
        let cm: CityJSONMetadata = from_str(&cm_str)?;

        let ctr_real_world: (f64, f64) = ((ctr_quantized.0 as f64 * cm.transform.scale[0]) + cm.transform.translate[0],
                                          (ctr_quantized.1 as f64 * cm.transform.scale[1]) + cm.transform.translate[1]);
        println!("real-world centroid: {:#?}", ctr_real_world);

        Ok(())
    }
}