use crate::geometry::Triangle;
use anyhow::Result;
use std::io::Write;
use byteorder::{WriteBytesExt, LittleEndian};

pub fn write_stl<W: Write>(writer: &mut W, triangles: &[Triangle]) -> Result<()> {
    // 80 bytes header
    let header = [0u8; 80];
    writer.write_all(&header)?;

    // Number of triangles (u32)
    let count = triangles.len() as u32;
    writer.write_u32::<LittleEndian>(count)?;

    for tri in triangles {
        // Normal
        writer.write_f32::<LittleEndian>(tri.normal.x)?;
        writer.write_f32::<LittleEndian>(tri.normal.y)?;
        writer.write_f32::<LittleEndian>(tri.normal.z)?;

        // V1
        writer.write_f32::<LittleEndian>(tri.v1.x)?;
        writer.write_f32::<LittleEndian>(tri.v1.y)?;
        writer.write_f32::<LittleEndian>(tri.v1.z)?;

        // V2
        writer.write_f32::<LittleEndian>(tri.v2.x)?;
        writer.write_f32::<LittleEndian>(tri.v2.y)?;
        writer.write_f32::<LittleEndian>(tri.v2.z)?;

        // V3
        writer.write_f32::<LittleEndian>(tri.v3.x)?;
        writer.write_f32::<LittleEndian>(tri.v3.y)?;
        writer.write_f32::<LittleEndian>(tri.v3.z)?;

        // Attribute byte count (u16)
        writer.write_u16::<LittleEndian>(0)?;
    }

    Ok(())
}
