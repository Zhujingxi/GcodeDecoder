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
    
    // Config
    default_width: f32,
    default_height: f32,
    
    // State for path building
    current_path: Option<ExtrusionPath>,
}

impl GCodeProcessor {
    pub fn new() -> Self {
        Self {
            current_pos: Vec3::ZERO,
            current_e: 0.0,
            relative_positioning: false,
            relative_extrusion: false, // Default is usually absolute M82
            default_width: 0.4, // Default nozzle diam
            default_height: 0.2, // Default layer height
            current_path: None,
        }
    }

    pub fn process(&mut self, lines: &[GCodeLine]) -> Vec<ExtrusionPath> {
        let mut paths = Vec::new();

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
        
        if has_move {
            if is_extruding {
                // If we don't have a current path, start one with the CURRENT position (start point)
                if self.current_path.is_none() {
                    self.current_path = Some(ExtrusionPath {
                        nodes: vec![PathNode {
                            pos: self.current_pos,
                            width: self.default_width,
                            height: self.default_height,
                        }]
                    });
                }
                
                // Add the target point
                if let Some(ref mut path) = self.current_path {
                     path.nodes.push(PathNode {
                        pos: target_pos,
                        width: self.default_width,
                        height: self.default_height,
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
