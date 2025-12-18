use std::time::Instant;
use glam::Vec3;
use gcodedecoder::config::Config;
use gcodedecoder::processor::{ExtrusionPath, PathNode};
use gcodedecoder::optimizer::optimize_paths;

fn main() {
    let mut config = Config::default();
    config.voxel_size_modifier = 0.25; 

    let mut paths = Vec::new();
    let gap = 0.3;
    let grid_size = 80; 
    
    println!("Generating {}x{}x{} grid of paths...", grid_size, grid_size, grid_size);
    for z in 0..grid_size {
        for y in 0..grid_size {
            for x in 0..grid_size {
                paths.push(ExtrusionPath {
                    nodes: vec![
                        PathNode {
                            pos: Vec3::new(x as f32 * gap, y as f32 * gap, z as f32 * gap),
                            width: config.nozzle_diameter,
                            height: config.layer_height,
                        },
                         PathNode {
                            pos: Vec3::new(x as f32 * gap, y as f32 * gap, (z as f32 * gap) + 0.2),
                            width: config.nozzle_diameter,
                            height: config.layer_height,
                        }
                    ]
                });
            }
        }
    }
    
    println!("Generated {} paths", paths.len());
    
    let start = Instant::now();
    let _optimized = optimize_paths(paths, &config);
    println!("Optimization took: {:.2?}", start.elapsed());
}
