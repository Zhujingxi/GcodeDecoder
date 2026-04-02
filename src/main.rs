use anyhow::{Context, Result};
use clap::Parser;
use std::fs;
use std::path::PathBuf;

use gcodedecoder::{config, export, geometry, parser, processor};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Input G-code file
    input: PathBuf,

    /// Output STL file
    #[arg(short, long)]
    output: Option<PathBuf>,

    // === Configuration Options ===
    /// Nozzle diameter in mm
    #[arg(long, default_value = "0.4")]
    nozzle_diameter: f32,

    /// Layer height in mm
    #[arg(long, default_value = "0.2")]
    layer_height: f32,

    /// Filament diameter in mm
    #[arg(long, default_value = "1.75")]
    filament_diameter: f32,

    /// Number of sides for tube mesh (4=low, 8=medium, 16=high quality)
    #[arg(long, default_value = "8")]
    mesh_sides: usize,
}

fn main() -> Result<()> {
    // Configure Rayon to use all available CPU cores for maximum parallelism
    let num_threads = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4);

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to configure Rayon thread pool");

    let cli = Cli::parse();

    // Build configuration from CLI arguments
    let config = config::ConfigBuilder::new()
        .nozzle_diameter(cli.nozzle_diameter)
        .layer_height(cli.layer_height)
        .filament_diameter(cli.filament_diameter)
        .mesh_sides(cli.mesh_sides)
        .build();

    println!("Using {} threads for parallel processing", num_threads);
    println!("Configuration:");
    println!(
        "  Nozzle: {}mm | Layer: {}mm | Filament: {}mm | Mesh Sides: {}",
        config.nozzle_diameter, config.layer_height, config.filament_diameter, config.mesh_sides
    );

    let content = fs::read_to_string(&cli.input)
        .with_context(|| format!("Failed to read input file: {:?}", cli.input))?;

    let segments = parser::parse_gcode(&content)?;

    let mut processor = processor::GCodeProcessor::with_config(&config);
    let processed_segments = processor.process(&segments);

    let mesh = geometry::generate_mesh(&processed_segments, &config);

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
