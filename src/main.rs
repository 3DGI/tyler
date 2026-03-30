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
mod cli;
mod formats;
mod parser;
mod proj;
mod spatial_structs;

use core::time::Duration;
use std::collections::HashSet;
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::formats::cesium3dtiles::{Tile, TileId};
use clap::Parser;
use log::{debug, info, log_enabled, warn, Level};
use rayon::prelude::*;
use subprocess::{Exec, Redirection};

#[derive(Debug, Default, Clone)]
struct SubprocessConfig {
    output_extension: String,
    exe: PathBuf,
    script: PathBuf,
    timeout: Option<Duration>,
    verbose: bool,
}

#[derive(Debug, Clone, clap::ValueEnum, Eq, PartialEq)]
#[clap(rename_all = "lower")]
pub enum Formats {
    _3DTiles,
    CityJSON,
}

impl ToString for Formats {
    fn to_string(&self) -> String {
        match self {
            Formats::_3DTiles => "3DTiles".to_string(),
            Formats::CityJSON => "CityJSON".to_string(),
        }
    }
}

#[derive(Default, Debug)]
struct DebugData {
    world: Option<PathBuf>,
    quadtree: Option<PathBuf>,
    tiles_results: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct PreparedInput {
    source: parser::InputSource,
    metadata_path: PathBuf,
    feature_base_document: Option<Vec<u8>>,
}

fn prepare_input(
    cli: &crate::cli::Cli,
    output_dir: &Path,
) -> Result<PreparedInput, Box<dyn std::error::Error>> {
    match cjindex::resolve_dataset(&cli.features, None) {
        Ok(resolved) => {
            if cli.metadata.is_some() {
                info!(
                    "Ignoring --metadata for cjindex dataset input; using the dataset metadata from {}",
                    resolved.dataset_root.display()
                );
            }
            let inspection = resolved.inspect()?;
            let mut city_index =
                cjindex::CityIndex::open(resolved.storage_layout(), &resolved.index_path)?;
            if !inspection.index.exists || inspection.index.fresh != Some(true) {
                info!(
                    "Rebuilding cjindex sidecar at {}",
                    resolved.index_path.display()
                );
                city_index.reindex()?;
            }
            let feature_base_document = derive_base_document(&city_index)?;
            let metadata_dir = output_dir.join("metadata");
            fs::create_dir_all(&metadata_dir)?;
            let metadata_path = metadata_dir.join("cjindex-metadata.city.json");
            fs::write(&metadata_path, &feature_base_document)?;
            Ok(PreparedInput {
                source: parser::InputSource::from_cjindex_resolved(&resolved),
                metadata_path,
                feature_base_document: Some(feature_base_document),
            })
        }
        Err(_error) => {
            let metadata_path = cli.metadata.clone().ok_or_else(|| {
                "--metadata is required when --features points at a legacy feature-file tree"
                    .to_string()
            })?;
            Ok(PreparedInput {
                source: parser::InputSource::LegacyFeatureFiles {
                    features_root: cli.features.clone(),
                },
                metadata_path,
                feature_base_document: None,
            })
        }
    }
}

fn derive_base_document(
    city_index: &cjindex::CityIndex,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let metadata = city_index.metadata()?;
    let Some(base_document) = metadata.first() else {
        return Err("cjindex dataset does not contain any source metadata".into());
    };
    if metadata
        .iter()
        .skip(1)
        .any(|candidate| candidate.as_ref() != base_document.as_ref())
    {
        return Err(
            "cjindex dataset contains multiple metadata documents; tyler requires one shared base document".into(),
        );
    }
    Ok(serde_json::to_vec(base_document.as_ref())?)
}

fn collect_tile_feature_ids(
    world: &parser::World,
    qtree_node: &spatial_structs::QuadTree,
) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut feature_ids = Vec::new();
    for cellid in qtree_node.cells() {
        let cell = world.grid.cell(cellid);
        for fid in &cell.feature_ids {
            if seen.insert(*fid) {
                feature_ids.push(*fid);
            }
        }
    }
    feature_ids
}

/// Write the list of feature paths for a tile into a text file, instead of passing
/// super long paths-string to the subprocess, because with very long arguments we can
/// get an 'Argument list too long' error.
fn write_inputs(
    world: &parser::World,
    path_features_input_dir: &Path,
    qtree_node: &spatial_structs::QuadTree,
    file_name: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let path_features_input_file = path_features_input_dir
        .join(file_name)
        .with_extension("input");
    let path_tile_ndjson = path_features_input_dir
        .join(file_name)
        .with_extension("city.jsonl");
    fs::create_dir_all(path_features_input_file.parent().unwrap()).unwrap_or_else(|_| {
        panic!(
            "should be able to create the directory {:?}",
            path_features_input_file.parent().unwrap()
        )
    });
    let ndjson_file = File::create(&path_tile_ndjson)
        .unwrap_or_else(|_| panic!("should be able to create a file {:?}", &path_tile_ndjson));
    let mut feature_output = BufWriter::new(ndjson_file);
    let feature_ids = collect_tile_feature_ids(world, qtree_node);

    match &world.input_source {
        parser::InputSource::LegacyFeatureFiles { features_root } => {
            for fid in feature_ids {
                let parser::FeatureReference::LegacyPath(relative_path) =
                    &world.features[fid].reference
                else {
                    return Err("legacy input unexpectedly referenced a cjindex feature".into());
                };
                let feature_path = features_root.join(relative_path);
                let bytes = fs::read(&feature_path)?;
                feature_output.write_all(&bytes)?;
                if !bytes.ends_with(b"\n") {
                    feature_output.write_all(b"\n")?;
                }
            }
        }
        parser::InputSource::CjIndexDataset { .. } => {
            let city_index = world.input_source.open_index()?;
            for fid in feature_ids {
                let parser::FeatureReference::CjIndexId(feature_id) =
                    &world.features[fid].reference
                else {
                    return Err(
                        "cjindex input unexpectedly referenced a legacy feature path".into(),
                    );
                };
                let model = city_index.get(feature_id)?.ok_or_else(|| {
                    format!("feature {feature_id} could not be resolved from cjindex")
                })?;
                cjlib::json::to_feature_writer(&mut feature_output, &model)?;
                feature_output.write_all(b"\n")?;
            }
        }
    }

    let _fi_file = File::create(&path_features_input_file).unwrap_or_else(|_| {
        panic!(
            "should be able to create a file {:?}",
            &path_features_input_file
        )
    });
    let mut feature_input = BufWriter::new(_fi_file);
    writeln!(feature_input, "{}", path_tile_ndjson.display())?;
    Ok(path_features_input_file)
}

fn run_subprocess(
    subprocess_config: &SubprocessConfig,
    tile: Tile,
    output_file: PathBuf,
    cmd: Exec,
) -> Option<Tile> {
    let cmd_string = cmd.to_cmdline_lossy();
    debug!("{cmd_string}");
    let redirection_stdout = Redirection::Pipe; // Redirection::Pipe | subprocess::NullFile
    let redirection_stderr = Redirection::Pipe; // Redirection::Merge
    let exec = cmd.stdout(redirection_stdout).stderr(redirection_stderr);
    let popen_res = exec.popen();
    match popen_res {
        Ok(mut popen) => {
            let (mut stdout_opt, mut stderr_opt): (Option<String>, Option<String>) = (None, None);
            let mut _exit_status = subprocess::ExitStatus::Undetermined;
            if let Some(timeout) = subprocess_config.timeout {
                let mut communicator = popen.communicate_start(None);
                if let Some(status) = popen.wait_timeout(timeout).unwrap() {
                    if let Ok(s) = communicator.read_string() {
                        (stdout_opt, stderr_opt) = s;
                    };
                    _exit_status = status;
                } else {
                    warn!(
                        "Tile {} timed out, conversion subprocess command:\n{}",
                        &tile.id, cmd_string
                    );
                    popen.kill().unwrap();
                    popen.wait().unwrap();
                    _exit_status = popen.exit_status().unwrap();
                }
            } else {
                (stdout_opt, stderr_opt) = popen.communicate(None).unwrap();
                _exit_status = popen.wait().unwrap();
            }

            // The stderr is Redirection::Merge-d into the stdout
            if !output_file.exists() {
                if subprocess_config.verbose {
                    warn!(
                        "Tile {} conversion failed, conversion subprocess command:\n{}\nsubprocess stdout:\n{}\nsubprocess stderr:\n{}",
                        tile.id, cmd_string, stdout_opt.unwrap_or_default(), stderr_opt.unwrap_or_default(),
                    );
                } else {
                    warn!(
                        "Tile {} conversion failed, conversion subprocess command:\n{}",
                        tile.id, cmd_string
                    );
                }
                return Some(tile);
            }
        }
        Err(popen_error) => {
            warn!("{}", popen_error);
            return Some(tile);
        }
    }
    None
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // --- Begin argument parsing
    let cli = crate::cli::Cli::parse();
    debug!("{:?}", &cli);
    info!("tyler version: {}", clap::crate_version!());
    if !cli.output.is_dir() {
        fs::create_dir_all(&cli.output)?;
        info!("Created output directory {:#?}", &cli.output);
    }
    // Since we have a default value, we can safely unwrap.
    let grid_cellsize = cli.grid_cellsize.unwrap();
    let geometric_error_above_leaf = cli.geometric_error_above_leaf.unwrap();
    let format = Formats::_3DTiles; // override --format
    let subprocess_config = match format {
        Formats::_3DTiles => {
            #[allow(unused)]
            let mut exe = PathBuf::new();
            if let Some(ref exe_g) = cli.exe_geof {
                assert!(exe_g.exists() && exe_g.is_file(), "geoflow executable must be an existing file for generating 3D Tiles, exe_geof: {:?}", &exe_g);
                exe = exe_g.clone();
            } else {
                debug!(
                    "exe_geof is not set for generating 3D Tiles, defaulting to 'geof' in the filesystem PATH"
                );
                exe = PathBuf::from("geof");
            }
            if !cli.cesium3dtiles_tileset_only {
                let res = Exec::cmd(&exe)
                    .arg("--version")
                    .arg("--verbose")
                    .stdout(Redirection::Pipe)
                    .stderr(Redirection::Merge)
                    .capture();
                let res_plugins = Exec::cmd(&exe)
                    .arg("--list-plugins")
                    .arg("--verbose")
                    .stdout(Redirection::Pipe)
                    .stderr(Redirection::Merge)
                    .capture();
                if let Ok(capture_data) = res {
                    let plugins_stdout_str = res_plugins.unwrap().stdout_str();
                    info!(
                        "geof version:\n{}{}",
                        capture_data.stdout_str(),
                        plugins_stdout_str
                    );
                } else if let Err(popen_error) = res {
                    panic!("Could not execute geof ({:?}):\n{}", &exe, popen_error)
                }
            }
            let geof_flowchart_path = match env::var("TYLER_RESOURCES_DIR") {
                Ok(val) => PathBuf::from(val).join("geof").join("createGLB.json"),
                Err(_) => PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("resources")
                    .join("geof")
                    .join("createGLB.json"),
            };
            let timeout = cli.timeout.map(|t| Duration::new(t, 0));
            SubprocessConfig {
                output_extension: "glb".to_string(),
                exe,
                script: geof_flowchart_path,
                timeout,
                verbose: cli.verbose_geof,
            }
        }
        Formats::CityJSON => {
            // TODO: refactor parallel loop
            panic!("cityjson output is not supported");
            // if let Some(exe) = cli.exe_python {
            //     SubprocessConfig {
            //         output_extension: "city.json".to_string(),
            //         exe,
            //         script: PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            //             .join("resources")
            //             .join("python")
            //             .join("convert_cityjsonfeatures.py"),
            //     }
            // } else {
            //     panic!("exe_python must be set for generating CityJSON tiles")
            // }
        }
    };
    debug!("{:?}", &subprocess_config);
    // Since we have a default value, it is safe to unwrap
    // let qtree_capacity = 0; // override cli.qtree_capacity
    let qtree_criteria = spatial_structs::QuadTreeCriteria::Vertices; // override --qtree-criteria
    let quadtree_capacity = match qtree_criteria {
        spatial_structs::QuadTreeCriteria::Objects => {
            spatial_structs::QuadTreeCapacity::Objects(cli.qtree_capacity.unwrap())
        }
        spatial_structs::QuadTreeCriteria::Vertices => {
            spatial_structs::QuadTreeCapacity::Vertices(cli.qtree_capacity.unwrap())
        }
    };
    let metadata_class: String = match format {
        Formats::_3DTiles => {
            if cli.cesium3dtiles_tileset_only {
                String::new()
            } else if cli.cesium3dtiles_metadata_class.is_none() {
                panic!("metadata_class must be set for writing 3D Tiles")
            } else {
                cli.cesium3dtiles_metadata_class.clone().unwrap()
            }
        }
        Formats::CityJSON => "".to_string(),
    };
    if cli.cesium3dtiles_content_bv_from_tile && !cli.cesium3dtiles_content_add_bv {
        warn!("cesium3dtiles_content_bv_from_tile is true, but cesium3dtiles_content_add_bv is false. The tile content bounding volumes are not going to be added, unless you set --3dtiles-content-add-bv");
    }
    let proj_data = match env::var("PROJ_DATA") {
        Ok(val) => {
            debug!("PROJ_DATA: {}", &val);
            Some(val)
        }
        Err(_val) => {
            warn!("PROJ_DATA environment variable is not set");
            None
        }
    };
    let debug_data = match cli.debug_load_data {
        None => DebugData::default(),
        Some(ref dir_path) => {
            if dir_path.is_dir() {
                let world_path = dir_path.join("world.bincode");
                let quadtree_path = dir_path.join("quadtree.bincode");
                let _tileset_path = dir_path.join("tileset.bincode");
                let tiles_results_path = dir_path.join("tiles_results.bincode");
                DebugData {
                    world: world_path.exists().then_some(world_path),
                    quadtree: quadtree_path.exists().then_some(quadtree_path),
                    tiles_results: tiles_results_path.exists().then_some(tiles_results_path),
                }
            } else {
                warn!(
                    "debug_load_data {dir_path:?} is not a directory, cannot load .bincode files"
                );
                DebugData::default()
            }
        }
    };
    debug!("{:?}", debug_data);
    let debug_data_output_path = cli.output.join("debug");
    if (cli.grid_export || log_enabled!(Level::Debug)) && !debug_data_output_path.exists() {
        fs::create_dir(&debug_data_output_path)?;
    }
    // --- end of argument parsing

    // Populate the World with features
    // Primitive types that implement Copy are efficiently copied into the function and
    // and it is cleaner to avoid the indirection. However, heap-allocated container
    // types are best passed by reference, because it is "expensive" to Clone them
    // (they don't implement Copy). When we move a value, we explicitly transfer
    // ownership of the value (eg cli.object_type).
    let prepared_input = if debug_data.world.is_none() {
        Some(prepare_input(&cli, &cli.output)?)
    } else {
        None
    };
    let cityobject_types = cli.object_type.clone();

    let world: parser::World = match debug_data.world {
        None => {
            let prepared_input = prepared_input
                .as_ref()
                .expect("prepared input must exist when world is built from source");
            let mut world = match &prepared_input.feature_base_document {
                Some(feature_base_document) => parser::World::from_cjindex(
                    prepared_input.source.clone(),
                    prepared_input.metadata_path.clone(),
                    feature_base_document.clone(),
                    grid_cellsize,
                    cityobject_types,
                    cli.grid_minz,
                    cli.grid_maxz,
                )?,
                None => parser::World::new(
                    &prepared_input.metadata_path,
                    &cli.features,
                    grid_cellsize,
                    cityobject_types,
                    cli.grid_minz,
                    cli.grid_maxz,
                )?,
            };
            world.index_with_grid()?; // todo input: in general, build a line index
            world
        }
        Some(world_path) => {
            info!("Loading world from bincode {world_path:?}");
            let world_file = File::open(world_path)?;
            bincode::deserialize_from(world_file)?
        }
    };

    info!(
        "Computed grid statistics: {}",
        world.grid.compute_statistics()
    );

    if cli.grid_export {
        info!("Exporting the grid to TSV to {:?}", &debug_data_output_path);
        world.export_grid(cli.grid_export_features, Some(&debug_data_output_path))?;
    }
    if log_enabled!(Level::Debug) {
        debug!(
            "Exporting the world instance to bincode to {:?}",
            &debug_data_output_path
        );
        world.export_bincode(Some("world"), Some(&debug_data_output_path))?;
    }

    // Build quadtree
    let quadtree: spatial_structs::QuadTree = match debug_data.quadtree {
        None => {
            info!("Building quadtree");
            spatial_structs::QuadTree::from_world(&world, quadtree_capacity)
        }
        Some(quadtree_path) => {
            info!("Loading quadtree from bincode {quadtree_path:?}");
            let quadtree_file = File::open(quadtree_path)?;
            bincode::deserialize_from(quadtree_file)?
        }
    };

    if cli.grid_export {
        info!(
            "Exporting the quadtree to TSV to {:?}",
            &debug_data_output_path
        );
        quadtree.export(&world, Some(&debug_data_output_path))?;
    }
    if log_enabled!(Level::Debug) {
        debug!(
            "Exporting the quadtree instance to bincode to {:?}",
            &debug_data_output_path
        );
        quadtree.export_bincode(Some("quadtree"), Some(&debug_data_output_path))?;
    }

    // 3D Tiles

    let tileset_path = cli.output.join("tileset.json");
    let subtrees_path = cli.output.join("subtrees");
    let tileset_path_unpruned = cli.output.join("tileset_unpruned.json");
    let subtrees_path_unpruned = cli.output.join("subtrees_unpruned");
    info!("Generating 3D Tiles tileset");
    let mut tileset = formats::cesium3dtiles::Tileset::from_quadtree(
        &quadtree,
        &world,
        geometric_error_above_leaf,
        grid_cellsize,
        cli.grid_minz,
        cli.grid_maxz,
        cli.cesium3dtiles_content_bv_from_tile,
        cli.cesium3dtiles_content_add_bv,
    );

    if cli.grid_export {
        info!(
            "Exporting the explicit tileset to TSV files to {:?}",
            &debug_data_output_path
        );
        tileset.export(Some(&debug_data_output_path))?;
    }

    let (tiles, _subtrees) = match cli.cesium3dtiles_implicit {
        true => {
            let mut tileset_implicit = tileset.clone();
            // FIXME: here we have a Vec<(Tile, TileId)> in 'tiles' instead of Vec<&Tile>, because of the
            //  mess with the implicit/explicit tile id-s.
            info!("Converting to implicit tiling");
            // Tileset.make_implicit() outputs the tiles that have content. If only the leaves have
            //  content, then only the leaves are outputted.
            let components: Vec<_> = subtrees_path_unpruned
                .components()
                .map(|comp| comp.as_os_str())
                .collect();
            let subtrees_dir_option = components.last().cloned().unwrap().to_str();
            let tiles_subtrees = tileset_implicit.make_implicit(
                &world.grid,
                &quadtree,
                cli.grid_export,
                subtrees_dir_option,
                Some(&debug_data_output_path),
            );

            if cli.cesium3dtiles_tileset_only || log_enabled!(Level::Debug) {
                info!("Writing unpruned 3D Tiles tileset");
                tileset_implicit.to_file(&tileset_path_unpruned)?;

                info!("Writing unpruned subtrees for implicit tiling");
                fs::create_dir_all(&subtrees_path_unpruned)?;
                for (subtree_id, subtree_bytes) in &tiles_subtrees.1 {
                    fs::create_dir_all(
                        subtrees_path_unpruned
                            .join(format!("{}/{}", subtree_id.level, subtree_id.x)),
                    )
                    .unwrap();
                    let out_path = subtrees_path_unpruned
                        .join(&subtree_id.to_string())
                        .with_extension("subtree");
                    let mut subtree_file = File::create(&out_path)
                        .unwrap_or_else(|_| panic!("could not create {:?} for writing", &out_path));
                    if let Err(_e) = subtree_file.write_all(subtree_bytes) {
                        warn!("Failed to write subtree {} content", subtree_id);
                    }
                }
            }

            tiles_subtrees
        }
        false => {
            let just_tiles = tileset.collect_leaves();
            // FIXME: here we need Vec<(Tile, TileId)> instead of Vec<&Tile>, for the same reason
            //  as above
            let tiles: Vec<(Tile, TileId)> = just_tiles
                .into_iter()
                .map(|tile_ref| (tile_ref.clone(), tile_ref.id.clone()))
                .collect();

            info!("Writing unpruned 3D Tiles tileset");
            tileset.to_file(&tileset_path_unpruned)?;

            (tiles, vec![])
        }
    };

    // Export by calling a subprocess to merge the .jsonl files and convert them to the
    // target format
    let cotypes_str: Vec<String> = match &world.cityobject_types {
        None => Vec::new(),
        Some(cotypes) => cotypes.iter().map(|co| co.to_string()).collect(),
    };
    let cotypes_arg = cotypes_str.join(",");

    let attribute_spec: String = match &cli.object_attribute {
        None => "".to_string(),
        Some(attributes) => attributes.join(","),
    };

    let path_output_tiles = cli.output.join("t");
    let path_features_input_dir = cli.output.join("inputs");
    // TODO: need to refactor this parallel loop somehow that it does not only read the
    //  3d tiles tiles, but also works with cityjson output
    if !cli.cesium3dtiles_tileset_only {
        fs::create_dir_all(&path_output_tiles)?;
        info!("Created output directory {:#?}", &path_output_tiles);
        fs::create_dir_all(&path_features_input_dir)?;
        info!("Created output directory {:#?}", &path_features_input_dir);

        let tiles_len = tiles.len();
        let tiles_failed_iter = tiles.into_par_iter().map(|(tile, tileid)| {
            #[allow(unused)]
            let mut tile_failed: Option<Tile> = None;
            let tileid_grid = &tile.id;
            let qtree_nodeid: spatial_structs::QuadTreeNodeId = tileid_grid.into();
            let qtree_node = quadtree
                .node(&qtree_nodeid)
                .unwrap_or_else(|| panic!("did not find tile {} in quadtree", tileid_grid));
            if qtree_node.nr_items == 0 {
                // The Tileset.prune() method removes the empty tiles from the tileset,
                //  so skipping the tile conversion without failure is ok if it's empty.
                debug!("Tile is empty ({}), skipping conversion", tileid_grid);
                return tile_failed;
            }
            let tileid_string = tileid.to_string();
            let file_name = tileid_string;
            let output_file = path_output_tiles
                .join(&file_name)
                .with_extension(&subprocess_config.output_extension);
            let path_features_input_file = match write_inputs(
                &world,
                &path_features_input_dir,
                qtree_node,
                file_name.as_str(),
            ) {
                Ok(path) => path,
                Err(error) => {
                    warn!(
                        "Failed to write NDJSON input for tile {}: {}",
                        tileid_grid, error
                    );
                    return Some(tile);
                }
            };

            // We use the quadtree node bbox here instead of the Tileset.Tile bounding
            // volume, because the Tile is in EPSG:4979 and we need the input data CRS
            let b = qtree_node.bbox(&world.grid);
            // We need to string-format all the arguments with an = separator, because that's what
            // geof can accept.
            // TODO: maybe replace the subprocess carte with std::process to remove the dependency
            let mut cmd = Exec::cmd(&subprocess_config.exe)
                .arg(&subprocess_config.script)
                .arg(format!(
                    "--output_format={}",
                    &format.to_string().to_lowercase()
                ))
                .arg(format!("--output_file={}", &output_file.to_str().unwrap()))
                .arg(format!(
                    "--path_metadata={}",
                    &world.path_metadata.to_str().unwrap()
                ))
                .arg(format!(
                    "--path_features_input_file={}",
                    &path_features_input_file.to_str().unwrap()
                ))
                .arg(format!("--min_x={}", b[0]))
                .arg(format!("--min_y={}", b[1]))
                .arg(format!("--min_z={}", b[2]))
                .arg(format!("--max_x={}", b[3]))
                .arg(format!("--max_y={}", b[4]))
                .arg(format!("--max_z={}", b[5]))
                .arg(format!("--cotypes={}", &cotypes_arg))
                .arg(format!("--metadata_class={}", &metadata_class))
                .arg(format!("--attribute_spec={}", &attribute_spec))
                .arg(format!("--geometric_error={}", &tile.geometric_error))
                .arg(format!("--bag3dBuildingsMode={}", cli.bag3d_buildings_mode))
                .arg(format!(
                    "--bag3dAttributesPerPart={}",
                    cli.bag3d_attributes_per_part
                ));

            if cli.verbose_geof {
                cmd = cmd.arg("--verbose".to_string())
            }

            if format == Formats::_3DTiles {
                // geof specific args
                // colors
                if cli.color_building.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBuilding={}",
                        cli.color_building.as_ref().unwrap()
                    ));
                }
                if cli.color_building_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBuildingPart={}",
                        cli.color_building_part.as_ref().unwrap()
                    ));
                }
                if cli.color_building_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBuildingInstallation={}",
                        cli.color_building_installation.as_ref().unwrap()
                    ));
                }
                if cli.color_tin_relief.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorTINRelief={}",
                        cli.color_tin_relief.as_ref().unwrap()
                    ));
                }
                if cli.color_road.is_some() {
                    cmd = cmd.arg(format!("--colorRoad={}", cli.color_road.as_ref().unwrap()));
                }
                if cli.color_railway.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorRailway={}",
                        cli.color_railway.as_ref().unwrap()
                    ));
                }
                if cli.color_transport_square.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorTransportSquare={}",
                        cli.color_transport_square.as_ref().unwrap()
                    ));
                }
                if cli.color_water_body.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorWaterBody={}",
                        cli.color_water_body.as_ref().unwrap()
                    ));
                }
                if cli.color_plant_cover.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorPlantCover={}",
                        cli.color_plant_cover.as_ref().unwrap()
                    ));
                }
                if cli.color_solitary_vegetation_object.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorSolitaryVegetationObject={}",
                        cli.color_solitary_vegetation_object.as_ref().unwrap()
                    ));
                }
                if cli.color_land_use.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorLandUse={}",
                        cli.color_land_use.as_ref().unwrap()
                    ));
                }
                if cli.color_city_furniture.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorCityFurniture={}",
                        cli.color_city_furniture.as_ref().unwrap()
                    ));
                }
                if cli.color_bridge.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBridge={}",
                        cli.color_bridge.as_ref().unwrap()
                    ));
                }
                if cli.color_bridge_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBridgePart={}",
                        cli.color_bridge_part.as_ref().unwrap()
                    ));
                }
                if cli.color_bridge_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBridgeInstallation={}",
                        cli.color_bridge_installation.as_ref().unwrap()
                    ));
                }
                if cli.color_bridge_construction_element.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorBridgeConstructionElement={}",
                        cli.color_bridge_construction_element.as_ref().unwrap()
                    ));
                }
                if cli.color_tunnel.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorTunnel={}",
                        cli.color_tunnel.as_ref().unwrap()
                    ));
                }
                if cli.color_tunnel_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorTunnelPart={}",
                        cli.color_tunnel_part.as_ref().unwrap()
                    ));
                }
                if cli.color_tunnel_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorTunnelInstallation={}",
                        cli.color_tunnel_installation.as_ref().unwrap()
                    ));
                }
                if cli.color_generic_city_object.is_some() {
                    cmd = cmd.arg(format!(
                        "--colorGenericCityObject={}",
                        cli.color_generic_city_object.as_ref().unwrap()
                    ));
                }

                // lod filter
                if cli.lod_building.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBuilding={}",
                        cli.lod_building.as_ref().unwrap()
                    ));
                }
                if cli.lod_building_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBuildingPart={}",
                        cli.lod_building_part.as_ref().unwrap()
                    ));
                }
                if cli.lod_building_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBuildingInstallation={}",
                        cli.lod_building_installation.as_ref().unwrap()
                    ));
                }
                if cli.lod_tin_relief.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodTINRelief={}",
                        cli.lod_tin_relief.as_ref().unwrap()
                    ));
                }
                if cli.lod_road.is_some() {
                    cmd = cmd.arg(format!("--lodRoad={}", cli.lod_road.as_ref().unwrap()));
                }
                if cli.lod_railway.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodRailway={}",
                        cli.lod_railway.as_ref().unwrap()
                    ));
                }
                if cli.lod_transport_square.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodTransportSquare={}",
                        cli.lod_transport_square.as_ref().unwrap()
                    ));
                }
                if cli.lod_water_body.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodWaterBody={}",
                        cli.lod_water_body.as_ref().unwrap()
                    ));
                }
                if cli.lod_plant_cover.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodPlantCover={}",
                        cli.lod_plant_cover.as_ref().unwrap()
                    ));
                }
                if cli.lod_solitary_vegetation_object.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodSolitaryVegetationObject={}",
                        cli.lod_solitary_vegetation_object.as_ref().unwrap()
                    ));
                }
                if cli.lod_land_use.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodLandUse={}",
                        cli.lod_land_use.as_ref().unwrap()
                    ));
                }
                if cli.lod_city_furniture.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodCityFurniture={}",
                        cli.lod_city_furniture.as_ref().unwrap()
                    ));
                }
                if cli.lod_bridge.is_some() {
                    cmd = cmd.arg(format!("--lodBridge={}", cli.lod_bridge.as_ref().unwrap()));
                }
                if cli.lod_bridge_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBridgePart={}",
                        cli.lod_bridge_part.as_ref().unwrap()
                    ));
                }
                if cli.lod_bridge_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBridgeInstallation={}",
                        cli.lod_bridge_installation.as_ref().unwrap()
                    ));
                }
                if cli.lod_bridge_construction_element.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodBridgeConstructionElement={}",
                        cli.lod_bridge_construction_element.as_ref().unwrap()
                    ));
                }
                if cli.lod_tunnel.is_some() {
                    cmd = cmd.arg(format!("--lodTunnel={}", cli.lod_tunnel.as_ref().unwrap()));
                }
                if cli.lod_tunnel_part.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodTunnelPart={}",
                        cli.lod_tunnel_part.as_ref().unwrap()
                    ));
                }
                if cli.lod_tunnel_installation.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodTunnelInstallation={}",
                        cli.lod_tunnel_installation.as_ref().unwrap()
                    ));
                }
                if cli.lod_generic_city_object.is_some() {
                    cmd = cmd.arg(format!(
                        "--lodGenericCityObject={}",
                        cli.lod_generic_city_object.as_ref().unwrap()
                    ));
                }

                if let Some(ref cotypes) = world.cityobject_types {
                    if cotypes.contains(&parser::CityObjectType::Building)
                        || cotypes.contains(&parser::CityObjectType::BuildingPart)
                    {
                        cmd = cmd.arg("--simplify_error=0.0").arg("--skip_clip=true");
                    } else if cli.simplification_max_error.is_some() {
                        cmd = cmd.arg(format!(
                            "--simplify_error={}",
                            cli.simplification_max_error.as_ref().unwrap()
                        ));
                    }
                }

                cmd = cmd.arg(format!("--smooth_normals={}", cli.smooth_normals));
            }

            if let Some(pd) = &proj_data {
                cmd = cmd.env("PROJ_DATA", pd);
            }

            tile_failed = run_subprocess(&subprocess_config, tile, output_file, cmd);
            tile_failed
        });

        let mut tiles_results: Vec<Option<Tile>> = Vec::with_capacity(tiles_len + 2);
        if let Some(tiles_results_path) = debug_data.tiles_results {
            info!("Loading tiles_results from {tiles_results_path:?}");
            let tiles_results_file = File::open(tiles_results_path)?;
            tiles_results = bincode::deserialize_from(tiles_results_file)?
        } else {
            info!("Converting and optimizing {tiles_len} tiles");
            tiles_failed_iter.collect_into_vec(&mut tiles_results);
            if log_enabled!(Level::Debug) {
                debug!(
                    "Exporting the tiles_results instance to bincode to {:?}",
                    &debug_data_output_path
                );
                let outpath = debug_data_output_path.join("tiles_results.bincode");
                let tiles_results_file = File::create(outpath)?;
                bincode::serialize_into(tiles_results_file, &tiles_results)?;
            }
        }
        let tiles_failed: Vec<Tile> = tiles_results.into_iter().flatten().collect();
        info!("Done");

        if !log_enabled!(Level::Debug) {
            fs::remove_dir_all(path_features_input_dir)?;
        }

        info!("Pruning tileset of {} failed tiles", tiles_failed.len());
        for (i, failed) in tiles_failed.iter().enumerate() {
            debug!("{}, removing failed from the tileset: {}", i, failed.id);
        }
        // Remove tiles that failed the gltf conversion
        tileset.prune(&tiles_failed, &quadtree);
        if cli.cesium3dtiles_implicit {
            // FIXME: here we re-create the implicit tileset from the pruned tileset,
            //  because it is simpler than flipping the bits of the unavailable tiles,
            //  because of the mixed up explicit/implicit tile IDs. But ideally, we
            //  flip the bits, so we won't need to duplicate the tileset here.
            let components: Vec<_> = subtrees_path
                .components()
                .map(|comp| comp.as_os_str())
                .collect();
            let subtrees_dir_option = components.last().cloned().unwrap().to_str();
            let (_, subtrees) = tileset.make_implicit(
                &world.grid,
                &quadtree,
                cli.grid_export,
                subtrees_dir_option,
                Some(&debug_data_output_path),
            );
            info!("Writing subtrees for implicit tiling");
            fs::create_dir_all(&subtrees_path)?;
            for (subtree_id, subtree_bytes) in subtrees {
                fs::create_dir_all(
                    subtrees_path.join(format!("{}/{}", subtree_id.level, subtree_id.x)),
                )
                .unwrap();
                let out_path = subtrees_path
                    .join(&subtree_id.to_string())
                    .with_extension("subtree");
                let mut subtree_file = File::create(&out_path)
                    .unwrap_or_else(|_| panic!("could not create {:?} for writing", &out_path));
                if let Err(_e) = subtree_file.write_all(&subtree_bytes) {
                    warn!("Failed to write subtree {} content", subtree_id);
                }
            }
        } else {
            let available_levels = tileset.available_levels();
            // A five level deep tree is still managable in size.
            if available_levels > 5 {
                // Try to find the split where each child tileset starts to have more tiles in their
                // tree, than the ancestor tree. This way, the main tileset is smaller in size than
                // the child tilesets, so it loads faster. This method is not very accurate, because
                // it doesn't account for the actual number of tiles on each level, it only
                // calculates with the theoretical maximum.
                let mut split_at_level = 0;
                for level in (0..available_levels).rev() {
                    let subtree_depth: u32 = (available_levels - level) as u32;
                    let nr_tiles_subtree = (4_usize.pow(subtree_depth) - 1) / 3;
                    let ancestor_tree_depth: u32 =
                        (available_levels - (available_levels - level)) as u32;
                    let nr_tiles_ancestor = (4_usize.pow(ancestor_tree_depth) - 1) / 3;
                    if nr_tiles_ancestor < nr_tiles_subtree {
                        split_at_level = level;
                        break;
                    }
                }
                info!(
                    "Splitting the explicit tileset into external tilesets at level {}",
                    split_at_level
                );
                let external_tilesets = tileset.split(split_at_level);
                for (filename, child_tileset) in &external_tilesets {
                    let tileset_path = cli.output.join(filename);
                    child_tileset.to_file(&tileset_path)?;
                }
            }
        }
        info!("Writing 3D Tiles tileset");
        tileset.to_file(&tileset_path)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_test_dir(prefix: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("tyler-{prefix}-{unique}"));
        fs::create_dir_all(&path).expect("create test dir");
        path
    }

    fn resource_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("data")
            .join(name)
    }

    fn build_quadtree(world: &parser::World) -> spatial_structs::QuadTree {
        spatial_structs::QuadTree::from_world(world, spatial_structs::QuadTreeCapacity::Objects(1))
    }

    fn exported_ndjson_path(input_file: &Path) -> PathBuf {
        fs::read_to_string(input_file)
            .expect("read input file")
            .trim()
            .into()
    }

    #[test]
    fn write_inputs_exports_legacy_features_as_ndjson() {
        let dataset_dir = unique_test_dir("legacy");
        let features_dir = dataset_dir.join("features");
        fs::create_dir_all(&features_dir).expect("create features dir");
        let metadata_path = dataset_dir.join("metadata.city.json");
        let feature_path = features_dir.join("sample.city.jsonl");
        fs::copy(resource_path("3dbag_x00.city.json"), &metadata_path).expect("copy metadata");
        fs::copy(resource_path("3dbag_feature_x71.city.jsonl"), &feature_path)
            .expect("copy feature");

        let mut world = parser::World::new(
            &metadata_path,
            &features_dir,
            200,
            Some(vec![parser::CityObjectType::Building]),
            None,
            None,
        )
        .expect("build legacy world");
        world.index_with_grid().expect("index legacy world");
        let quadtree = build_quadtree(&world);
        let inputs_dir = dataset_dir.join("inputs");
        let input_file =
            write_inputs(&world, &inputs_dir, &quadtree, "tile").expect("write inputs");
        let ndjson_path = exported_ndjson_path(&input_file);
        let ndjson = fs::read_to_string(ndjson_path).expect("read exported ndjson");

        assert!(ndjson.contains("\"type\":\"CityJSONFeature\""));
        assert_eq!(ndjson.lines().count(), 1);
    }

    #[test]
    fn write_inputs_exports_cjindex_ndjson_as_ndjson() {
        let dataset_dir = unique_test_dir("cjindex-ndjson");
        let metadata =
            fs::read_to_string(resource_path("3dbag_x00.city.json")).expect("read metadata");
        let feature = fs::read_to_string(resource_path("3dbag_feature_x71.city.jsonl"))
            .expect("read feature");
        let ndjson_source = dataset_dir.join("source.city.jsonl");
        fs::write(&ndjson_source, format!("{metadata}\n{feature}\n")).expect("write ndjson source");

        let resolved =
            cjindex::resolve_dataset(&dataset_dir, None).expect("resolve ndjson dataset");
        let mut city_index =
            cjindex::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex ndjson dataset");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            Some(vec![parser::CityObjectType::Building]),
            None,
            None,
        )
        .expect("build cjindex ndjson world");
        world.index_with_grid().expect("index cjindex ndjson world");
        let quadtree = build_quadtree(&world);
        let inputs_dir = dataset_dir.join("inputs");
        let input_file =
            write_inputs(&world, &inputs_dir, &quadtree, "tile").expect("write inputs");
        let ndjson_path = exported_ndjson_path(&input_file);
        let ndjson = fs::read_to_string(ndjson_path).expect("read exported ndjson");

        assert!(ndjson.contains("\"type\":\"CityJSONFeature\""));
        assert_eq!(ndjson.lines().count(), 1);
    }

    #[test]
    fn write_inputs_exports_cjindex_cityjson_as_ndjson() {
        let dataset_dir = unique_test_dir("cjindex-cityjson");
        let metadata: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_x00.city.json")).expect("read metadata"),
        )
        .expect("parse metadata");
        let feature: Value = serde_json::from_slice(
            &fs::read(resource_path("3dbag_feature_x71.city.jsonl")).expect("read feature"),
        )
        .expect("parse feature");
        let mut cityjson = metadata;
        cityjson["CityObjects"] = feature["CityObjects"].clone();
        cityjson["vertices"] = feature["vertices"].clone();
        let cityjson_path = dataset_dir.join("source.city.json");
        fs::write(
            &cityjson_path,
            serde_json::to_vec(&cityjson).expect("serialize cityjson"),
        )
        .expect("write cityjson source");

        let resolved =
            cjindex::resolve_dataset(&dataset_dir, None).expect("resolve cityjson dataset");
        let mut city_index =
            cjindex::CityIndex::open(resolved.storage_layout(), &resolved.index_path)
                .expect("open index");
        city_index.reindex().expect("reindex cityjson dataset");
        let feature_base_document = derive_base_document(&city_index).expect("derive base doc");
        let metadata_path = dataset_dir.join("metadata.city.json");
        fs::write(&metadata_path, &feature_base_document).expect("write metadata");

        let mut world = parser::World::from_cjindex(
            parser::InputSource::from_cjindex_resolved(&resolved),
            metadata_path,
            feature_base_document,
            200,
            Some(vec![parser::CityObjectType::Building]),
            None,
            None,
        )
        .expect("build cjindex cityjson world");
        world
            .index_with_grid()
            .expect("index cjindex cityjson world");
        let quadtree = build_quadtree(&world);
        let inputs_dir = dataset_dir.join("inputs");
        let input_file =
            write_inputs(&world, &inputs_dir, &quadtree, "tile").expect("write inputs");
        let ndjson_path = exported_ndjson_path(&input_file);
        let ndjson = fs::read_to_string(ndjson_path).expect("read exported ndjson");

        assert!(ndjson.contains("\"type\":\"CityJSONFeature\""));
        assert_eq!(ndjson.lines().count(), 1);
    }
}
