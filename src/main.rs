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
use std::env;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

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
    tiles_failed: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // --- Begin argument parsing
    let cli = crate::cli::Cli::parse();
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
            let mut exe = PathBuf::new();
            if let Some(exe_g) = cli.exe_geof {
                assert!(exe_g.exists() && exe_g.is_file(), "geoflow executable must be an existing file for generating 3D Tiles, exe_geof: {:?}", &exe_g);
                exe = exe_g;
            } else {
                debug!(
                    "exe_geof is not set for generating 3D Tiles, defaulting to 'geof' in the filesystem PATH"
                );
                exe = PathBuf::from("geof");
            }
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
            let geof_flowchart_path = match env::var("TYLER_RESOURCES_DIR") {
                Ok(val) => PathBuf::from(val).join("geof").join("createGLB.json"),
                Err(_) => PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("resources")
                    .join("geof")
                    .join("createGLB.json"),
            };
            let timeout = cli.timeout.map(|_t| Duration::new(5, 0));
            SubprocessConfig {
                output_extension: "glb".to_string(),
                exe,
                script: geof_flowchart_path,
                timeout,
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
            if cli.cesium3dtiles_metadata_class.is_none() {
                panic!("metadata_class must be set for writing 3D Tiles")
            } else {
                cli.cesium3dtiles_metadata_class.unwrap()
            }
        }
        Formats::CityJSON => "".to_string(),
    };
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
        Some(dir_path) => {
            if dir_path.is_dir() {
                let world_path = dir_path.join("world.bincode");
                let quadtree_path = dir_path.join("quadtree.bincode");
                let _tileset_path = dir_path.join("tileset.bincode");
                let tiles_failed_path = dir_path.join("tiles_failed.bincode");
                DebugData {
                    world: world_path.exists().then_some(world_path),
                    quadtree: quadtree_path.exists().then_some(quadtree_path),
                    tiles_failed: tiles_failed_path.exists().then_some(tiles_failed_path),
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
    // --- end of argument parsing

    // Populate the World with features
    // Primitive types that implement Copy are efficiently copied into the function and
    // and it is cleaner to avoid the indirection. However, heap-allocated container
    // types are best passed by reference, because it is "expensive" to Clone them
    // (they don't implement Copy). When we move a value, we explicitly transfer
    // ownership of the value (eg cli.object_type).
    let world: parser::World = match debug_data.world {
        None => {
            let mut world = parser::World::new(
                &cli.metadata,
                &cli.features,
                grid_cellsize,
                cli.object_type,
                cli.grid_minz,
                cli.grid_maxz,
            )?;
            world.index_with_grid();
            world
        }
        Some(world_path) => {
            info!("Loading world from {world_path:?}");
            let world_file = File::open(world_path)?;
            bincode::deserialize_from(world_file)?
        }
    };

    if cli.grid_export {
        info!("Exporting the grid to the working directory");
        world.export_grid(cli.grid_export_features)?;
    }
    if log_enabled!(Level::Debug) {
        debug!("Exporting the world instance to the working directory");
        world.export_bincode(Some("world"))?;
    }

    // Build quadtree
    let quadtree: spatial_structs::QuadTree = match debug_data.quadtree {
        None => {
            info!("Building quadtree");
            spatial_structs::QuadTree::from_world(&world, quadtree_capacity)
        }
        Some(quadtree_path) => {
            info!("Loading quadtree from {quadtree_path:?}");
            let quadtree_file = File::open(quadtree_path)?;
            bincode::deserialize_from(quadtree_file)?
        }
    };

    if cli.grid_export {
        info!("Exporting the quadtree to the working directory");
        quadtree.export(&world.grid)?;
    }
    if log_enabled!(Level::Debug) {
        debug!("Exporting the quadtree instance to the working directory");
        quadtree.export_bincode(Some("quadtree"))?;
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
    );

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
            );

            if log_enabled!(Level::Debug) {
                debug!("Writing unpruned 3D Tiles tileset");
                tileset_implicit.to_file(&tileset_path_unpruned)?;

                debug!("Writing unpruned subtrees for implicit tiling");
                fs::create_dir_all(&subtrees_path_unpruned)?;
                for (subtree_id, subtree_bytes) in &tiles_subtrees.1 {
                    fs::create_dir_all(
                        subtrees_path.join(format!("{}/{}", subtree_id.level, subtree_id.x)),
                    )
                    .unwrap();
                    let out_path = subtrees_path
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

    let path_output_tiles = cli.output.join("tiles");
    let path_features_input_dir = cli.output.join("inputs");
    // TODO: need to refactor this parallel loop somehow that it does not only read the
    //  3d tiles tiles, but also works with cityjson output
    if !cli.cesium3dtiles_tileset_only {
        fs::create_dir_all(&path_output_tiles)?;
        info!("Created output directory {:#?}", &path_output_tiles);
        fs::create_dir_all(&path_features_input_dir)?;
        info!("Created output directory {:#?}", &path_features_input_dir);

        let tiles_len = tiles.len();
        let tiles_failed_iter = tiles
            .into_par_iter()
            .map(|(tile, tileid)| {
                let mut tile_failed: Option<Tile> = None;
                let tileid_grid = &tile.id;
                let qtree_nodeid: spatial_structs::QuadTreeNodeId = tileid_grid.into();
                let qtree_node = quadtree
                    .node(&qtree_nodeid)
                    .unwrap_or_else(|| panic!("did not find tile {} in quadtree", tileid_grid));
                let tileid_string = tileid.to_string();
                let file_name = tileid_string;
                let output_file = path_output_tiles
                    .join(&file_name)
                    .with_extension(&subprocess_config.output_extension);
                // We write the list of feature paths for a tile into a text file, instead of passing
                // super long paths-string to the subprocess, because with very long arguments we can
                // get an 'Argument list too long' error.
                let path_features_input_file = path_features_input_dir
                    .join(&file_name)
                    .with_extension("input");
                fs::create_dir_all(path_features_input_file.parent().unwrap()).unwrap_or_else(
                    |_| {
                        panic!(
                            "should be able to create the directory {:?}",
                            path_features_input_file.parent().unwrap()
                        )
                    },
                );
                let mut feature_input =
                    File::create(&path_features_input_file).unwrap_or_else(|_| {
                        panic!(
                            "should be able to create a file {:?}",
                            &path_features_input_file
                        )
                    });
                for cellid in qtree_node.cells() {
                    let cell = world.grid.cell(cellid);
                    for fid in cell.feature_ids.iter() {
                        let fp = world.features[*fid]
                            .path_jsonl
                            .clone()
                            .into_os_string()
                            .into_string()
                            .unwrap();
                        writeln!(feature_input, "{}", fp)
                            .expect("should be able to write feature path to the input file");
                    }
                }

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
                    .arg(format!("--geometric_error={}", &tile.geometric_error));

                if format == Formats::_3DTiles {
                    // geof specific args
                    // verbose. we always run geof in --verbose, because otherwise it won't print 
                    // the logs when we get an error from geof
                    cmd = cmd.arg("--verbose");
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
                            cmd = cmd.arg("--simplify_ratio=1.0").arg("--skip_clip=true");
                        } else if cli.reduce_vertices.is_some() {
                            cmd = cmd.arg(format!(
                                "--simplify_ratio={}",
                                cli.reduce_vertices.as_ref().unwrap()
                            ));
                        }
                    }
                }

                if let Some(pd) = &proj_data {
                    cmd = cmd.env("PROJ_DATA", pd);
                }

                let cmd_string = cmd.to_cmdline_lossy();
                let exec = cmd.stdout(Redirection::Pipe).stderr(Redirection::Merge);
                let popen_res = exec.popen();
                match popen_res {
                    Ok(mut popen) => {
                        let (mut stdout_opt, mut stderr_opt): (Option<String>, Option<String>) =
                            (None, None);
                        let mut exit_status = subprocess::ExitStatus::Undetermined;
                        if let Some(timeout) = subprocess_config.timeout {
                            let mut communicator = popen.communicate_start(None);
                            if let Some(status) = popen.wait_timeout(timeout).unwrap() {
                                if let Ok(s) = communicator.read_string() {
                                    (stdout_opt, stderr_opt) = s;
                                };
                                // (stdout_opt, stderr_opt) = popen.communicate(None).unwrap();
                                exit_status = status;
                            } else {
                                warn!(
                                    "tile {} timed out, conversion subprocess command:\n{}",
                                    &tile.id, cmd_string
                                );
                                popen.kill().unwrap();
                                popen.wait().unwrap();
                                exit_status = popen.exit_status().unwrap();
                            }
                        } else {
                            (stdout_opt, stderr_opt) = popen.communicate(None).unwrap();
                            exit_status = popen.wait().unwrap();
                        }

                        let stdout = stdout_opt.unwrap_or_default();
                        // The stderr is Redirection::Merge-d into the stdout
                        if !exit_status.success() {
                            warn!("{} conversion subprocess failed\ncommand: {}\nwith stdout and stderr:\n{}", &tileid, &cmd_string, &stdout);
                        }
                        if !output_file.exists() {
                            warn!(
                                "{} output {:?} was not written by the subprocess, conversion subprocess command:\n{}",
                                &tileid, &output_file, &cmd_string
                            );
                            tile_failed = Some(tile);
                        }
                    }
                    Err(popen_error) => {
                        warn!("{}", popen_error);
                        tile_failed = Some(tile);
                    }
                }
                tile_failed
            })
            .filter_map(|failed_tile| failed_tile);

        let tiles_failed: Vec<Tile> = match debug_data.tiles_failed {
            None => {
                info!("Exporting and optimizing {tiles_len} tiles");
                tiles_failed_iter.collect()
            }
            Some(tiles_failed_path) => {
                info!("Loading tiles_failed from {tiles_failed_path:?}");
                let tiles_failed_file = File::open(tiles_failed_path)?;
                bincode::deserialize_from(tiles_failed_file)?
            }
        };
        info!("Done");

        if !log_enabled!(Level::Debug) {
            fs::remove_dir_all(path_features_input_dir)?;
        } else {
            let tiles_failed_file = File::create("tiles_failed.bincode")?;
            bincode::serialize_into(tiles_failed_file, &tiles_failed)?;
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
            let (_, subtrees) =
                tileset.make_implicit(&world.grid, &quadtree, cli.grid_export, subtrees_dir_option);
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
        }
        info!("Writing 3D Tiles tileset");
        tileset.to_file(&tileset_path)?;
    }

    Ok(())
}
