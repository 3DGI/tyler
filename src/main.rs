mod parser;
mod quadtree;

use std::path::{Path};
use std::fs::{self, DirEntry};
use std::io;

use clap::{crate_version, Arg, Command};
use walkdir::WalkDir;

static FORMATS: [&str; 1] = [
    "3dtiles",
];


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cmd = Command::new("tyler")
        .about("Create tiles from 3D city objects encoded as CityJSONFeatures.")
        .version(crate_version!())
        .arg(
            Arg::new("metadata")
                .value_name("METADATA")
                .short('m')
                .long("metadata")
                .required(true)
                .help("Main CityJSON file (.city.json), containing the coordinate reference system and transformation properties.")
        )
        .arg(
            Arg::new("features")
                .value_name("FEATURES")
                .short('f')
                .long("features")
                .required(true)
                .help("Directory of CityJSONFeatures (.city.json.).")
        )
        .arg(
            Arg::new("output")
                .value_name("OUTPUT")
                .short('o')
                .long("output")
                .required(true)
                .help("Directory for the output.")
        )
        .arg(
            Arg::new("format")
                .value_name("FORMAT")
                .long("format")
                .required(true)
                .help("Output format.")
                .value_parser(FORMATS)
        );
    let matches = cmd.get_matches();

    let path_metadata = Path::new(matches.get_one::<String>("metadata").expect("required"))
        .canonicalize()
        .expect("Could not find the METADATA file.");
    if !path_metadata.is_file() { panic!("METADATA must be an existing file") }
    let path_features = Path::new(matches.get_one::<String>("features").expect("required"))
        .canonicalize()
        .expect("Could not find the FEATURES directory.");
    if !path_features.is_dir() { panic!("FEATURES must be an existing directory") }
    let path_output = Path::new(matches.get_one::<String>("output").expect("required")).to_path_buf();
    if !path_output.is_dir() {
        fs::create_dir_all(&path_output)?;
        println!("Created output directory {:#?}", &path_output);
    }

    let cm: parser::CityJSONMetadata = parser::parse_metadata(&path_metadata)?;

    let feature_set_iter = WalkDir::new(&path_features).into_iter()
        .filter_map(|res| {
            if let Ok(entry) = res {
                if let Some(ext) = entry.path().extension() {
                    if ext == "jsonl" {
                        Some(entry.path().to_path_buf())
                    } else { None }
                } else { None }
            } else { None }
        })
        .map(|path| parser::feature_to_tuple(path, &cm))
        .filter_map(|res| res.ok());

    let _: Vec<parser::FeatureTuple> = feature_set_iter.collect();

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_visit_dir() {
        let p = Path::new("/home/balazs/Development/cjio_dbexport/tests/data/output");
        let mut cnt = 0;
        for entry in WalkDir::new(p) {
            if let Some(ext) = entry.unwrap().path().extension() {
                    if ext == "jsonl" {
                        cnt = cnt + 1;
                    }
                }
        }
        println!("{}", cnt);
    }
}