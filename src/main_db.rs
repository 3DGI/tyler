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
mod cli_db;
mod formats;
mod parser;
mod proj;
mod spatial_structs;

use std::env;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::formats::cesium3dtiles::{Tile, TileId};
use clap::Parser;
use log::{debug, error, info, log_enabled, warn, Level};
use rayon::prelude::*;
use subprocess::{Exec, Redirection};

#[derive(Debug, Default, Clone)]
struct SubprocessConfig {
    output_extension: String,
    exe: PathBuf,
    script: PathBuf,
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // --- Begin argument parsing
    let cli = crate::cli_db::Cli::parse();
    // Since we have a default value, we can safely unwrap.
    let grid_cellsize = cli.grid_cellsize.unwrap();
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
    // --- end of argument parsing

    let mut world = parser::WorldDb::new(
        cli.uri.as_ref(),
        cli.table.as_ref(),
        cli.geometry_column.as_ref(),
        cli.primary_key.as_ref(),
        grid_cellsize,
    )?;
    world.index_with_grid()?;

    // Build quadtree
    info!("Building quadtree");
    let quadtree = spatial_structs::QuadTree::from_worlddb(&world, quadtree_capacity);

    Ok(())
}
