use std::path::{Path, PathBuf};

use clap::Parser;

#[derive(Parser)]
#[command(author, version, about)]
pub struct Cli {
    /// Main CityJSON file (.city.json), containing the coordinate reference system and
    /// transformation properties.
    #[arg(short, long, value_parser = existing_canonical_path)]
    pub metadata: PathBuf,
    /// Directory of CityJSONFeatures (.city.jsonl). The directory and all its
    /// subdirectories are searched recursively for feature files.
    #[arg(short, long, value_parser = existing_canonical_path)]
    pub features: PathBuf,
    /// Directory for the output.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Output format.
    #[arg(long, value_parser = FORMATS)]
    pub format: String,
    /// The CityObject type to use for the 3D Tiles
    /// (https://www.cityjson.org/specs/1.1.3/#the-different-city-objects).
    #[arg(long, value_enum)]
    pub object_type: Option<Vec<crate::parser::CityObjectType>>,
    /// Export the grid and the feature centroids in to .tsv files in the working
    /// directory. Used for debugging.
    #[arg(long)]
    pub grid_export: bool,
    /// Set the cell size for the grid that is used for constructing the quadtree.
    #[arg(long, default_value = "200")]
    pub grid_cellsize: Option<u16>,
    /// Limit the minimum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely small z
    /// coordinates.
    #[arg(long)]
    pub grid_minz: Option<i32>,
    /// Limit the maximum z coordinate for the bounding box that is computed from the
    /// features. Useful if the features contain errors with extremely large z
    /// coordinates.
    #[arg(long)]
    pub grid_maxz: Option<i32>,
    /// What to count for the quadtree leaf capacity.
    #[arg(long, value_enum, default_value = "vertices")]
    pub qtree_capacity_type: Option<crate::spatial_structs::QuadTreeCapacityType>,
    /// The capacity of a leaf of the quadtree. If a quadrant has less than or equal
    /// the capacity, its subtiles are merged.
    #[arg(long, default_value = "10000")]
    pub qtree_capacity: Option<usize>,
    /// Path to the gltfpack executable (https://meshoptimizer.org/gltf/).
    #[arg(long, value_parser = existing_path)]
    pub exe_gltfpack: Option<PathBuf>,
    /// Path to the geoflow executable for clipping and exporting the gltf files.
    #[arg(long, value_parser = existing_path)]
    pub exe_geof: Option<PathBuf>,
    /// Path to the python interpreter (>=3.8) to use for generating CityJSON tiles.
    /// The interpreter must have a recent cjio (https://github.com/cityjson/cjio)
    /// installed.
    #[arg(long, value_parser = existing_path)]
    pub exe_python: Option<PathBuf>,
}

// Would be more elegant to have an enum for this, but we cannot have enum variants
// that begin with a number.
static FORMATS: [&str; 2] = ["3dtiles", "cityjson"];

fn existing_canonical_path(s: &str) -> Result<PathBuf, String> {
    if let Ok(c) = Path::new(s).canonicalize() {
        if c.exists() {
            Ok(c)
        } else {
            Err(format!("path {:?} does not exist", &c))
        }
    } else {
        Err(format!("could not resolve the path {:?}", s))
    }
}

/// We don't want to canonicalize paths to executables, especially a python exe from a
/// virtualenv, because the symlink would get resolved and we would end up with a path
/// to the python interpreter that was used for creating the virtualenv, and not the
/// interpreter that links to the virtualenv.
fn existing_path(s: &str) -> Result<PathBuf, String> {
    let p = Path::new(s).to_path_buf();
    if p.exists() {
        Ok(p)
    } else {
        Err(format!("path {:?} does not exist", &p))
    }
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};

    fn required_args() -> Vec<&'static str> {
        vec![
            "tyler",
            "-m",
            "metadata.city.json",
            "-f",
            env!("CARGO_MANIFEST_DIR"),
            "-o",
            env!("CARGO_MANIFEST_DIR"),
            "--format",
            "3dtiles",
        ]
    }

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert()
    }

    /// Can we pass multiple CityObject types?
    #[test]
    fn verify_object_types() {
        let mut types: Vec<&'static str> =
            vec!["--object-type", "Building", "--object-type", "PlantCover"];
        let mut args = required_args();
        args.append(&mut types);
        let cli = Cli::try_parse_from(args).unwrap();
        let otypes = &cli.object_type.unwrap();
        assert!(otypes.contains(&crate::parser::CityObjectType::Building));
        assert!(otypes.contains(&crate::parser::CityObjectType::PlantCover));
    }
}
