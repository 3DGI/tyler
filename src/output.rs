use std::path::PathBuf;

use crate::TileCoord;

/// Naming shared by published files and the references to those files.
#[derive(Clone, Copy)]
pub struct OutputPaths<'a> {
    pub tiles_dir: &'a str,
    pub suffix: Option<&'a str>,
}

impl OutputPaths<'_> {
    pub fn filename(&self, stem: &str, extension: &str) -> String {
        cityjson_convert::output_filename(stem, extension, self.suffix)
    }

    pub fn tile_path(&self, tile: &TileCoord, extension: &str) -> String {
        format!(
            "{}/{}/{}",
            tile.level,
            tile.x,
            self.filename(&format!("{}-{}-{}", tile.level, tile.x, tile.y), extension)
        )
    }

    pub fn tile_uri(&self, tile: &TileCoord, extension: &str) -> String {
        format!("{}/{}", self.tiles_dir, self.tile_path(tile, extension))
    }

    pub fn tsv_uri(&self, tile: &TileCoord) -> String {
        format!(
            "{}/{}/{}/{}-{}-{}/{}",
            self.tiles_dir,
            tile.level,
            tile.x,
            tile.level,
            tile.x,
            tile.y,
            self.filename("cityobjects", "tsv")
        )
    }

    pub fn implicit_content_uri(&self) -> String {
        format!(
            "{}/{{level}}/{{x}}/{}",
            self.tiles_dir,
            self.filename("{level}-{x}-{y}", "glb")
        )
    }

    pub fn subtree_uri(&self, directory: &str) -> String {
        format!(
            "{directory}/{{level}}/{{x}}/{}",
            self.filename("{y}", "subtree")
        )
    }

    pub fn subtree_path(&self, tile: &crate::formats::cesium3dtiles::TileId) -> PathBuf {
        PathBuf::from(tile.level.to_string())
            .join(tile.x.to_string())
            .join(self.filename(&tile.y.to_string(), "subtree"))
    }
}
