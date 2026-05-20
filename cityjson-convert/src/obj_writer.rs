use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use cityjson_lib::CityModel;
use log::debug;

use crate::triangle_mesh::{build_triangle_mesh, TriangleMeshOptions};
use crate::ObjExportOptions;

pub fn write_city_model_obj<P: AsRef<Path>>(
    model: &CityModel,
    output_path: P,
    options: &ObjExportOptions,
) -> Result<()> {
    let mesh = build_triangle_mesh(
        model,
        &TriangleMeshOptions {
            clip_bbox: options.clip_bbox,
            clip_geographic_region: None,
        },
    )?;

    let file = File::create(output_path.as_ref())?;
    let mut writer = BufWriter::new(file);
    let mut next_index = 1_u64;
    let mut vertex_count = 0_u64;
    let mut face_count = 0_u64;

    for object in mesh.objects {
        if object.triangles.is_empty() {
            continue;
        }
        writeln!(writer, "o {}", sanitize_object_name(&object.object_id))?;
        for triangle in object.triangles {
            for position in triangle.source_positions {
                writeln!(writer, "v {} {} {}", position[0], position[1], position[2])?;
            }
            writeln!(
                writer,
                "f {} {} {}",
                next_index,
                next_index + 1,
                next_index + 2
            )?;
            next_index += 3;
            vertex_count += 3;
            face_count += 1;
        }
    }

    writer.flush()?;
    debug!(
        "OBJ Summary: {} vertices and {} faces written to {}",
        vertex_count,
        face_count,
        output_path.as_ref().display()
    );
    Ok(())
}

fn sanitize_object_name(value: &str) -> String {
    let sanitized = value
        .chars()
        .map(|ch| {
            if ch.is_whitespace() || ch.is_control() {
                '_'
            } else {
                ch
            }
        })
        .collect::<String>();
    if sanitized.is_empty() {
        "_".to_string()
    } else {
        sanitized
    }
}

#[cfg(test)]
mod tests {
    use super::sanitize_object_name;

    #[test]
    fn sanitize_replaces_whitespace_and_control_characters() {
        assert_eq!(sanitize_object_name("a b\tc\n"), "a_b_c_");
    }
}
