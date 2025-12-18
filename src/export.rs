use crate::geometry::Triangle;
use anyhow::Result;
use std::io::{Write, BufWriter};
use byteorder::{WriteBytesExt, LittleEndian};

pub fn write_stl<W: Write>(writer: &mut W, triangles: &[Triangle]) -> Result<()> {
    // Wrap in BufWriter for better I/O performance (64KB buffer)
    let mut buf = BufWriter::with_capacity(64 * 1024, writer);
    
    // 80 bytes header
    let header = [0u8; 80];
    buf.write_all(&header)?;

    // Number of triangles (u32)
    let count = triangles.len() as u32;
    buf.write_u32::<LittleEndian>(count)?;

    for tri in triangles {
        // Normal
        buf.write_f32::<LittleEndian>(tri.normal.x)?;
        buf.write_f32::<LittleEndian>(tri.normal.y)?;
        buf.write_f32::<LittleEndian>(tri.normal.z)?;

        // V1
        buf.write_f32::<LittleEndian>(tri.v1.x)?;
        buf.write_f32::<LittleEndian>(tri.v1.y)?;
        buf.write_f32::<LittleEndian>(tri.v1.z)?;

        // V2
        buf.write_f32::<LittleEndian>(tri.v2.x)?;
        buf.write_f32::<LittleEndian>(tri.v2.y)?;
        buf.write_f32::<LittleEndian>(tri.v2.z)?;

        // V3
        buf.write_f32::<LittleEndian>(tri.v3.x)?;
        buf.write_f32::<LittleEndian>(tri.v3.y)?;
        buf.write_f32::<LittleEndian>(tri.v3.z)?;

        // Attribute byte count (u16)
        buf.write_u16::<LittleEndian>(0)?;
    }

    buf.flush()?;
    Ok(())
}
