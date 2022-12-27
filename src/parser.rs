use serde::Deserialize;
use zerovec::{VarZeroVec, ZeroSlice};
use std::path::PathBuf;

/// A partial [CityJSON object](https://www.cityjson.org/specs/1.1.3/#cityjson-object).
/// It is partial, because we only store the metadata that is necessary for parsing the
/// CityJSONFeatures.
#[derive(Deserialize, Debug)]
struct CityJSONMetadata {
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
/// It borrows from the JSON string that is loaded into memory. In order to achieve this zero-copy
/// (or zero-allocation) deserialization of vectors with serde, we need the
/// [zerovec](https://crates.io/crates/zerovec) crate.
#[derive(Deserialize, Debug)]
struct CityJSONFeatureVertices<'a> {
    #[serde(borrow)]
    vertices: VarZeroVec<'a, ZeroSlice<i64>>,
}


impl CityJSONFeatureVertices<'_> {
    fn to_feature_tuple(&self) {
        todo!()
    }

    /// Return the number of vertices of the feature.
    /// We assume that the number of vertices in a feature does not exceed 65535 (thus `u16`).
    fn vertex_count(&self) -> u16 {
        self.vertices.len() as u16
    }

    /// Feature centroid (2D) computed as the average coordinate.
    fn centroid(&self, transform: &Transform) -> (f64, f64) {
        let mut x_sum: f64 = 0.0;
        let mut y_sum: f64 = 0.0;
        for v in self.vertices.iter() {
            let x = v.get(0).unwrap() as f64;
            let y = v.get(1).unwrap() as f64;
            x_sum = x_sum + ((x * transform.scale[0]) + transform.translate[0]);
            y_sum = y_sum + ((y * transform.scale[1]) + transform.translate[1]);
        }
        let x = x_sum / self.vertices.len() as f64;
        let y = y_sum / self.vertices.len() as f64;
        (x, y)
    }
}

type FeatureTuple = (u64, u16, PathBuf);

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
        let pb: PathBuf = test_data_dir().join("3dbag_x00.city.json");
        let cm_str = read_to_string(&pb).unwrap();
        let cm: CityJSONMetadata = from_str(&cm_str)?;

        let pb: PathBuf = test_data_dir().join("3dbag_feature_x71.city.jsonl");
        let cf_str = read_to_string(&pb).unwrap();
        let cf: CityJSONFeatureVertices = from_str(&cf_str)?;
        let ctr_real_world = cf.centroid(&cm.transform);
        println!("real-world centroid: {:#?}", ctr_real_world);

        Ok(())
    }
}