use serde::Deserialize;
use zerovec::{VarZeroVec, ZeroSlice};

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
    reference_system: String
}

/// Container for storing the CityJSONFeature vertices.
/// It borrows from the JSON string that is loaded into memory. In order to achieve this zero-copy
/// (or zero-allocation) deserialization of vectors with serde, we need the
/// [zerovec](https://crates.io/crates/zerovec) crate.
#[derive(Deserialize, Debug)]
struct CityJSONFeatureVertices<'a> {
    #[serde(borrow)]
    vertices: VarZeroVec<'a, ZeroSlice<i64>>
}


#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::from_str;

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
}