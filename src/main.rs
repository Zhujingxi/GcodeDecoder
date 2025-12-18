use anyhow::{Context, Result};
use std::fs;
use std::path::PathBuf;
use clap::Parser;

use gcodedecoder::{parser, processor, geometry, export};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Input G-code file
    input: PathBuf,

    /// Output STL file
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    let content = fs::read_to_string(&cli.input)
        .with_context(|| format!("Failed to read input file: {:?}", cli.input))?;

    // Logic will be similar to lib.rs but without WASM types
    // We might want to move the high-level orchestration to a shared function in lib.rs 
    // to avoid strict duplication, but for now duplicate to keep lib.rs clean for WASM.
    
    let segments = parser::parse_gcode(&content)?;
    
    let mut processor = processor::GCodeProcessor::new();
    let processed_segments = processor.process(&segments);
    
    let mesh = geometry::generate_mesh(&processed_segments);
    
    let output_path = cli.output.unwrap_or_else(|| {
        let mut p = cli.input.clone();
        p.set_extension("stl");
        p
    });

    let file = fs::File::create(&output_path)?;
    let mut writer = std::io::BufWriter::new(file);
    export::write_stl(&mut writer, &mesh)?;
    
    println!("Exported to {:?}", output_path);

    Ok(())
}
