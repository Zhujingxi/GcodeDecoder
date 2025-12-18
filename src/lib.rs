pub mod parser;
pub mod processor;
pub mod geometry;
pub mod export;

use wasm_bindgen::prelude::*;
use processor::GCodeProcessor;



#[wasm_bindgen]
pub fn set_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn convert_gcode_to_stl(gcode: &str) -> Result<Vec<u8>, JsValue> {
    let segments = parser::parse_gcode(gcode)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
        
    let mut processor = GCodeProcessor::new();
    let processed_segments = processor.process(&segments);
    
    let mesh = geometry::generate_mesh(&processed_segments);
    
    let mut buffer = Vec::new();
    export::write_stl(&mut buffer, &mesh)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
        
    Ok(buffer)
}
