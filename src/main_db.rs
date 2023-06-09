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
use postgres::{Client, Error, NoTls};
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
    // FIXME: need to sanitize the user input database relations
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
        &cli.uri,
        &cli.table,
        &cli.geometry_column,
        &cli.primary_key,
        grid_cellsize,
    )?;
    world.index_with_grid()?;

    // Build quadtree
    info!("Building quadtree");
    let quadtree = spatial_structs::QuadTree::from_worlddb(&world, quadtree_capacity);

    world.export_grid()?;
    quadtree.export(&world.grid)?;

    info!("Writing quadtree leaves to the database");
    let mut client = Client::connect(&cli.uri, NoTls)?;
    let table_tiles = format!("{}.tiles", cli.output_schema);
    let table_index = format!("{}.index", cli.output_schema);
    if cli.drop_existing {
        let q = format!("DROP TABLE IF EXISTS {table_index} CASCADE",);
        client.batch_execute(&q)?;
        let q = format!("DROP TABLE IF EXISTS {table_tiles} CASCADE",);
        client.batch_execute(&q)?;
    }
    let q = format!(
        "CREATE TABLE {table_tiles} (tile_id text, cnt bigint, boundary geometry(Polygon, {}))",
        world.grid.epsg
    );
    let q =
        format!("COMMENT ON TABLE {table_tiles} IS 'Tile boundaries, generated with tyler-db.'",);
    client.batch_execute(&q)?;
    let q = format!(
        "CREATE TABLE {table_index} ({} int, tile_id text)",
        &cli.primary_key
    );
    client.batch_execute(&q)?;
    let q = format!(
        "COMMENT ON TABLE {table_index} IS 'Feature-Tile index, generated with tyler-db.'",
    );
    client.batch_execute(&q)?;
    let q_tiles = format!(
        "INSERT INTO {table_tiles} VALUES ($1, $2, st_geomfromtext($3, {}))",
        world.grid.epsg
    );
    let statement_tiles = client.prepare(&q_tiles)?;
    let statement_copy = client.prepare(format!("COPY {table_index} FROM STDIN").as_str())?;

    let leaves = quadtree.collect_leaves();
    for leaf in leaves {
        let wkt = leaf.to_wkt(&world.grid);
        let tile_id = leaf.id.to_string();
        let nr_items_i64 = leaf.nr_items as i64;
        let mut nr_features_inserted = 0;
        let mut data = String::new();
        for cellid in leaf.cells() {
            let cell = world.grid.cell(cellid);
            for feature_id in &cell.feature_ids {
                let pk = world.features[*feature_id].primary_key;
                data.push_str(format!("{}\t{}\n", &pk, &tile_id).as_str());
                nr_features_inserted += 1;
            }
        }
        if nr_features_inserted > 0 {
            let mut writer = client.copy_in(&statement_copy)?;
            writer.write_all(data.as_bytes())?;
            writer.finish()?;
            client.execute(&statement_tiles, &[&tile_id, &nr_items_i64, &wkt])?;
        }
    }
    info!("Done");
    Ok(())
}
