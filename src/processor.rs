use crate::config::Config;
use crate::parser::GCodeLine;
use glam::Vec3;

#[derive(Debug, Clone)]
pub struct ExtrusionPath {
    pub nodes: Vec<PathNode>,
}

#[derive(Debug, Clone, Copy)]
pub struct PathNode {
    pub pos: Vec3,
    pub width: f32,
    pub height: f32,
}

pub struct GCodeProcessor {
    current_pos: Vec3,
    current_e: f32,
    
    // Modes
    relative_positioning: bool, // G91
    relative_extrusion: bool, // M83
    
    // Config parameters (from Config)
    default_width: f32,
    default_height: f32,
    filament_diameter: f32,
    min_layer_height: f32,
    max_layer_height: f32,
    
    // State for path building
    current_path: Option<ExtrusionPath>,
    first_extrusion_handled: bool,
}

impl GCodeProcessor {
    /// Create a new processor with the given configuration
    pub fn with_config(config: &Config) -> Self {
        Self {
            current_pos: Vec3::ZERO,
            current_e: 0.0,
            relative_positioning: false,
            relative_extrusion: false, // Default is usually absolute M82
            default_width: config.nozzle_diameter,
            default_height: config.layer_height,
            filament_diameter: config.filament_diameter,
            min_layer_height: config.min_layer_height(),
            max_layer_height: config.max_layer_height(),
            current_path: None,
            first_extrusion_handled: false,
        }
    }

    /// Create a new processor with default configuration
    pub fn new() -> Self {
        Self::with_config(&Config::default())
    }

    pub fn process(&mut self, lines: &[GCodeLine]) -> Vec<ExtrusionPath> {
        // Pre-allocate: estimate ~1 path per 20 G-code lines
        let mut paths = Vec::with_capacity(lines.len() / 20 + 1);

        for line in lines {
            match line.command.as_str() {
                "G0" | "G1" => {
                    self.process_move(line, &mut paths);
                }
                "G90" => self.relative_positioning = false,
                "G91" => self.relative_positioning = true,
                "M82" => self.relative_extrusion = false,
                "M83" => self.relative_extrusion = true,
                "G92" => {
                    // Reset positions
                    for (param, val) in &line.params {
                        match param {
                            'X' => self.current_pos.x = *val,
                            'Y' => self.current_pos.y = *val,
                            'Z' => self.current_pos.z = *val,
                            'E' => self.current_e = *val,
                            _ => {}
                        }
                    }
                    // Break current path on reset
                    self.finish_path(&mut paths);
                }
                _ => {}
            }
        }
        self.finish_path(&mut paths);

        paths
    }

    fn finish_path(&mut self, paths: &mut Vec<ExtrusionPath>) {
        if let Some(path) = self.current_path.take() {
            if path.nodes.len() > 1 {
                paths.push(path);
            }
        }
    }

    fn process_move(&mut self, line: &GCodeLine, paths: &mut Vec<ExtrusionPath>) {
        let mut target_pos = self.current_pos;
        let mut target_e = self.current_e;
        let mut has_move = false;
        let mut has_extrude = false;

        for (param, val) in &line.params {
            match param {
                'X' => {
                    if self.relative_positioning { target_pos.x += val } else { target_pos.x = *val }
                    has_move = true;
                }
                'Y' => {
                    if self.relative_positioning { target_pos.y += val } else { target_pos.y = *val }
                    has_move = true;
                }
                'Z' => {
                    if self.relative_positioning { target_pos.z += val } else { target_pos.z = *val }
                    has_move = true;
                }
                'E' => {
                    if self.relative_extrusion { target_e += val } else { target_e = *val }
                    has_extrude = true;
                }
                _ => {}
            }
        }

        let is_extruding = has_extrude && target_e > self.current_e;
        
        // Adaptive Layer Height Logic
        // If this is the FIRST extrusion, set the default height to the current Z level.
        // This handles cases where the first layer is 0.1, 0.24, 0.3 but we default to 0.2.
        if is_extruding && !self.first_extrusion_handled {
            // Check if current Z is reasonable for a layer height (derived from nozzle diameter)
            if self.current_pos.z > self.min_layer_height && self.current_pos.z < self.max_layer_height {
                self.default_height = self.current_pos.z;
                // println!("Adjusted base layer height to: {}", self.default_height); // Debug
            }
            self.first_extrusion_handled = true;
        }

        // Calculate extrusion logic
        // Use default width unless calculated otherwise
        let mut segment_width = self.default_width;
        let segment_height = self.default_height;

        if is_extruding && has_move {
            let dist = target_pos.distance(self.current_pos);
            let extrude_amount = target_e - self.current_e;
            
            // Try to use momentary Z height if it looks like a layer height?
            // Actually, simply using the detected default_height is safer than trusting Z-pos 
            // because of Z-hops. 
            // Ideally we tracked Z-per-layer, but for now fixed first layer detection handles the "Bottom Layer" user issue.

            if dist > 0.0001 && extrude_amount > 0.0 {
                 // Volumetric Calculation
                 // Filament Volume = PI * r^2 * length
                 let r = self.filament_diameter / 2.0;
                 let vol_in = std::f32::consts::PI * r * r * extrude_amount;
                 
                 // Line Volume = Width * Height * Distance 
                 segment_width = vol_in / (segment_height * dist);
            }
        }

        if has_move {
            if is_extruding {
                // If we don't have a current path, start one with the CURRENT position (start point)
                if self.current_path.is_none() {
                    self.current_path = Some(ExtrusionPath {
                        nodes: vec![PathNode {
                            pos: self.current_pos,
                            width: segment_width, // Start node width? 
                            // Using the segment width for the start node is a reasonable approximation for the first segment
                            height: segment_height,
                        }]
                    });
                }
                
                // Add the target point
                if let Some(ref mut path) = self.current_path {
                     path.nodes.push(PathNode {
                        pos: target_pos,
                        width: segment_width,
                        height: segment_height,
                    });
                }
            } else {
                // Travel move: finish current path and do NOT start a new one (until we extrude again)
                self.finish_path(paths);
            }
            
            self.current_pos = target_pos;
        }
        
        if has_extrude {
            self.current_e = target_e;
        }
    }
}
