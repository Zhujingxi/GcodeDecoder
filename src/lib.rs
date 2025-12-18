pub mod config;
pub mod parser;
pub mod processor;
pub mod geometry;
pub mod export;
pub mod optimizer;

use wasm_bindgen::prelude::*;
use config::{Config, ConfigBuilder};
use processor::GCodeProcessor;

#[wasm_bindgen]
pub fn set_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Convert G-code to STL with default configuration (backward compatible)
#[wasm_bindgen]
pub fn convert_gcode_to_stl(gcode: &str) -> Result<Vec<u8>, JsValue> {
    convert_gcode_to_stl_with_config_internal(gcode, &Config::default())
}

/// Convert G-code to STL with custom configuration
/// 
/// Parameters:
/// - gcode: G-code content as string
/// - nozzle_diameter: Nozzle diameter in mm (e.g., 0.4)
/// - layer_height: Layer height in mm (e.g., 0.2)
/// - filament_diameter: Filament diameter in mm (e.g., 1.75)
/// - mesh_sides: Number of sides for tube mesh (e.g., 4, 8, 16)
#[wasm_bindgen]
pub fn convert_gcode_to_stl_with_config(
    gcode: &str,
    nozzle_diameter: f32,
    layer_height: f32,
    filament_diameter: f32,
    mesh_sides: usize,
) -> Result<Vec<u8>, JsValue> {
    let config = ConfigBuilder::new()
        .nozzle_diameter(nozzle_diameter)
        .layer_height(layer_height)
        .filament_diameter(filament_diameter)
        .mesh_sides(mesh_sides)
        .build();
    
    convert_gcode_to_stl_with_config_internal(gcode, &config)
}

fn convert_gcode_to_stl_with_config_internal(gcode: &str, config: &Config) -> Result<Vec<u8>, JsValue> {
    let segments = parser::parse_gcode(gcode)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
        
    let mut processor = GCodeProcessor::with_config(config);
    let processed_segments = processor.process(&segments);
    
    let mesh = geometry::generate_mesh(&processed_segments, config);
    
    let mut buffer = Vec::new();
    export::write_stl(&mut buffer, &mesh)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
        
    Ok(buffer)
}


