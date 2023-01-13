mod formats;
mod parser;
mod proj;
mod spatial_structs;

use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use clap::{crate_version, Arg, ArgAction, Command};
use log::{debug, error, info};
use rayon::prelude::*;
use subprocess::{Exec, Redirection};
use walkdir::WalkDir;

static FORMATS: [&str; 2] = ["3dtiles", "cityjson"];

type FeatureSet = Vec<parser::Feature>;

/// 3D bounding box.
/// [min x, min y, min z, max x, max y, max z]
pub type Bbox = [f64; 6];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

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
        )
        .arg(
            Arg::new("export")
                .long("export-grid")
                .action(ArgAction::SetTrue)
                .help("Export the grid and the feature centroids in to .tsv files in the working directory.")
        )
        .arg(
            Arg::new("cellsize")
                .long("cellsize")
                .default_value("200")
                .value_parser(clap::value_parser!(u16))
                .help("Set the cell size for the square grid.")
        )
        .arg(
            Arg::new("qtree_limit")
                .long("quadtree-limit")
                .default_value("10000")
                .value_parser(clap::value_parser!(usize))
                .help("The threshold for merging cells in the quadtree. If a quadrant has <= the limit its subtiles are merged.")
        )
        .arg(
            Arg::new("qtree_criteria")
                .long("quadtree-criteria")
                .required(true)
                .value_parser(["objects", "vertices"])
                .help("The criteria to use when evaluating the quadtree-limit.")
        )
        .arg(
            Arg::new("python")
                .long("python-bin")
                .default_value("python3")
                .value_parser(clap::value_parser!(String))
                .help("Path to the python interpreter (>=3.8) to use for converting the CityJSONFeatures to the target format. The interpreter must have a recent [cjio](https://github.com/cityjson/cjio) installed.")
        );
    let matches = cmd.get_matches();

    let path_metadata = Path::new(matches.get_one::<String>("metadata").expect("required"))
        .canonicalize()
        .expect("Could not find the METADATA file.");
    if !path_metadata.is_file() {
        panic!("METADATA must be an existing file")
    }
    let path_features = Path::new(matches.get_one::<String>("features").expect("required"))
        .canonicalize()
        .expect("Could not find the FEATURES directory.");
    if !path_features.is_dir() {
        panic!("FEATURES must be an existing directory")
    }
    let path_output =
        Path::new(matches.get_one::<String>("output").expect("required")).to_path_buf();
    if !path_output.is_dir() {
        fs::create_dir_all(&path_output)?;
        info!("Created output directory {:#?}", &path_output);
    }
    let output_format = matches
        .get_one::<String>("format")
        .expect("format is required")
        .as_str();
    let do_export = matches.get_one::<bool>("export").copied().unwrap();
    let cellsize: u16 = *matches
        .get_one::<u16>("cellsize")
        .expect("could not parse the 'cellsize' argument to u16");
    let python_bin = matches
        .get_one::<String>("python")
        .expect("could not parse the python interpreter path from the arguments")
        .as_str();
    if !Path::new(&python_bin).is_file() {
        panic!("Python interpreter does not exist {}", python_bin);
    }
    let ql: usize = *matches
        .get_one::<usize>("qtree_limit")
        .expect("could not parse the quadtree-limit from the arguments");
    let qt_crit = matches
        .get_one::<String>("qtree_criteria")
        .expect("quadtree-criteria is required")
        .as_str();
    let mut quadtree_limit = spatial_structs::QuadTreeLimit::Objects(4000 as usize);
    match qt_crit {
        "objects" => quadtree_limit = spatial_structs::QuadTreeLimit::Objects(ql),
        "vertices" => quadtree_limit = spatial_structs::QuadTreeLimit::Vertices(ql),
        &_ => {}
    }

    let cm = parser::CityJSONMetadata::from_file(&path_metadata)?;

    // Return the file path if the 'DirEntry' is a .jsonl file (eg. .city.jsonl).
    let jsonl_path_closure = |res: Result<walkdir::DirEntry, walkdir::Error>| {
        if let Ok(entry) = res {
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
    };

    // Do a first loop over the features to calculate their extent and their number.
    // Need a mutable iterator, because .next() consumes the next value and advances the iterator.
    let mut features_enum_iter = WalkDir::new(&path_features)
        .into_iter()
        .filter_map(jsonl_path_closure)
        .enumerate();
    // Init the extent with from the first feature.
    let (mut nr_features, feature_path) = features_enum_iter
        .next()
        .expect(".jsonl file should be accessible");
    let cf = parser::CityJSONFeatureVertices::from_file(&feature_path)?;
    let mut extent_qc = cf.bbox();
    for (nf, fp) in features_enum_iter {
        let cf = parser::CityJSONFeatureVertices::from_file(&fp)?;
        let [x_min, y_min, z_min, x_max, y_max, z_max] = cf.bbox();
        if x_min < extent_qc[0] {
            extent_qc[0] = x_min
        } else if x_max > extent_qc[3] {
            extent_qc[3] = x_max
        }
        if y_min < extent_qc[1] {
            extent_qc[1] = y_min
        } else if y_max > extent_qc[4] {
            extent_qc[4] = y_max
        }
        if z_min < extent_qc[2] {
            extent_qc[2] = z_min
        } else if z_max > extent_qc[5] {
            extent_qc[5] = z_max
        }
        nr_features = nf;
    }
    info!("Found {} features", nr_features);
    debug!("extent_qc: {:?}", &extent_qc);
    // Get the real-world coordinates for the extent
    let extent_rw_min = extent_qc[0..3]
        .into_iter()
        .enumerate()
        .map(|(i, qc)| (*qc as f64 * cm.transform.scale[i]) + cm.transform.translate[i]);
    let extent_rw_max = extent_qc[3..6]
        .into_iter()
        .enumerate()
        .map(|(i, qc)| (*qc as f64 * cm.transform.scale[i]) + cm.transform.translate[i]);
    let extent_rw: [f64; 6] = extent_rw_min
        .chain(extent_rw_max)
        .collect::<Vec<f64>>()
        .try_into()
        .expect("should be able to create an [f64; 6] from the extent_rw vector");
    debug!("extent real-world: {:?}", &extent_rw);

    // Init the grid from the extent
    let epsg = cm.metadata.reference_system.to_epsg()?;
    let mut grid = spatial_structs::SquareGrid::new(&extent_rw, cellsize, epsg, Some(10.0));
    debug!("{}", grid);

    let feature_set_iter = WalkDir::new(&path_features)
        .into_iter()
        .filter_map(jsonl_path_closure)
        .map(|feature_path| parser::CityJSONFeatureVertices::file_to_tuple(&feature_path))
        .filter_map(|res| res.ok());
    let mut feature_set: FeatureSet = Vec::with_capacity(nr_features);
    for (fid, mut feature) in feature_set_iter.enumerate() {
        let centroid = feature.centroid(&cm);
        grid.insert(&centroid, fid);
        feature_set.push(feature);
    }

    // Debug
    if do_export {
        debug!("exporting grid to working directory");
        grid.export(&feature_set, &cm)?;
    }

    // Build quadtree
    info!("Building quadtree");
    let quadtree = spatial_structs::QuadTree::from_grid(&grid, &feature_set, quadtree_limit);

    // 3D Tiles
    info!("Generating 3D Tiles tileset");
    let tileset_path = path_output.join("tileset.json");
    let tileset =
        formats::cesium3dtiles::Tileset::from_quadtree(&quadtree, &grid, &cm, &feature_set);
    tileset.to_file(tileset_path)?;

    let path_output_tiles = path_output.join("tiles");
    if !path_output_tiles.is_dir() {
        fs::create_dir_all(&path_output_tiles)?;
        info!("Created output directory {:#?}", &path_output_tiles);
    }
    let path_features_input_dir = path_output.join("inputs");
    if !path_features_input_dir.is_dir() {
        fs::create_dir_all(&path_features_input_dir)?;
        info!("Created output directory {:#?}", &path_features_input_dir);
    }

    // Export by calling a python subprocess to merge the .jsonl files and convert them to the
    // target format
    let python_script = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("resources")
        .join("python")
        .join("convert_cityjsonfeatures.py");

    let output_extension = match output_format {
        "3dtiles" => "glb",
        "cityjson" => "city.json",
        _ => "unknown",
    };

    let mut cellids: Vec<spatial_structs::CellId> = Vec::with_capacity(grid.length * grid.length);
    cellids = grid.into_iter().map(|(cellid, _cell)| cellid).collect();

    let leaves: Vec<&spatial_structs::QuadTree> = quadtree.collect_leaves();
    info!("Exporting {} tiles", leaves.len());
    leaves.into_par_iter().for_each(|tile| {
        if tile.nr_items > 0 {
            let tileid = tile.id();
            let file_name = format!("{}", &tileid);
            let output_file = path_output_tiles
                .join(&file_name)
                .with_extension(output_extension);
            // We write the list of feature paths for a tile into a text file, instead of passing
            // super long paths-string to the subprocess, because with very long arguments we can
            // get an 'Argument list too long' error.
            let path_features_input_file = path_features_input_dir
                .join(&file_name)
                .with_extension("input");
            fs::create_dir_all(path_features_input_file.parent().unwrap()).unwrap_or_else(|_| {
                panic!(
                    "should be able to create the directory {:?}",
                    path_features_input_file.parent().unwrap()
                )
            });
            let mut feature_input = File::create(&path_features_input_file).unwrap_or_else(|_| {
                panic!(
                    "should be able to create a file {:?}",
                    &path_features_input_file
                )
            });
            for cellid in tile.cells.iter() {
                let cell = grid.cell(cellid);
                for fid in cell.iter() {
                    let fp = feature_set[*fid]
                        .path_jsonl
                        .clone()
                        .into_os_string()
                        .into_string()
                        .unwrap();
                    writeln!(feature_input, "{}", fp)
                        .expect("should be able to write feature path to the input file");
                }
            }

            let b = tile.bbox(&grid);
            let bbox = format!("{},{},{},{}", b[0], b[1], b[3], b[4]);

            debug!("converting {}", &tileid);
            let res_exit_status = Exec::cmd(python_bin)
                .arg(&python_script)
                .arg(output_format)
                .arg(output_file)
                .arg(&path_metadata)
                .arg(&path_features_input_file)
                .arg(bbox)
                .stdout(Redirection::Pipe)
                .stderr(Redirection::Merge)
                .capture();
            if let Ok(capturedata) = res_exit_status {
                let stdout = capturedata.stdout_str();
                if !capturedata.success() {
                    error!("{} subprocess stdout: {}", &tileid, stdout);
                    error!(
                        "{} subprocess stderr: {}",
                        &tileid,
                        capturedata.stderr_str()
                    );
                } else if !stdout.is_empty() && stdout != "\n" {
                    debug!("{} subproces stdout {}", &tileid, capturedata.stdout_str());
                }
            } else if let Err(popen_error) = res_exit_status {
                error!("{}", popen_error);
            }
        } else {
            debug!("tile {} is empty", &tile.id())
        }
    });
    info!("Done");
    Ok(())
}
