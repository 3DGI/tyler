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
use std::path::{Path, PathBuf};

use clap::Parser;

#[derive(Parser)]
#[command(author, version, about)]
pub struct Cli {
    /// PostgreSQL connection URI (postgresql://[userspec@][hostspec][/dbname][?paramspec]).
    #[arg(long)]
    pub uri: String,
    /// Database table name.
    #[arg(long)]
    pub table: String,
    /// Database table geometry column.
    #[arg(long)]
    pub geometry_column: String,
    /// Database table primary key column.
    #[arg(long)]
    pub primary_key: String,
    /// The schema where to write the 'tiles' and 'index' tables.
    #[arg(long)]
    pub output_schema: String,
    /// Drop the 'tiles' and 'index' tables if they exist.
    #[arg(long)]
    pub drop_existing: bool,
    /// Set the 2D cell size for the grid that is used for constructing the quadtree. In input units (eg. meters).
    #[arg(long, default_value = "250")]
    pub grid_cellsize: Option<u16>,
    /// The maximum number of vertices in a leaf of the quadtree.
    #[arg(long, default_value = "18000")]
    pub qtree_capacity: Option<usize>,
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert()
    }
}
