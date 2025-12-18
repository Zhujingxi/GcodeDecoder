use crate::config::Config;
use crate::processor::ExtrusionPath;
use glam::Vec3;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

pub struct Triangle {
    pub v1: Vec3,
    pub v2: Vec3,
    pub v3: Vec3,
    pub normal: Vec3,
}

pub fn generate_mesh(paths: &[ExtrusionPath], config: &Config) -> Vec<Triangle> {
    let num_sides = config.mesh_sides; // From quality preset
    
    // Precompute unit circle (Cos, Sin) table
    let unit_circle: Vec<(f32, f32)> = (0..num_sides)
        .map(|j| {
            let angle = (j as f32 / num_sides as f32) * std::f32::consts::TAU;
            (angle.cos(), angle.sin())
        })
        .collect();

    // Helper closure to process one path
    let process_path = |path: &ExtrusionPath| -> Vec<Triangle> {
        if path.nodes.len() < 2 {
            return Vec::new();
        }
        
        let node_count = path.nodes.len();
        // Pre-allocate: (nodes-1) segments * num_sides quads * 2 tris + 2 caps * num_sides tris
        let estimated_tris = (node_count - 1) * num_sides * 2 + num_sides * 2;
        let mut local_triangles = Vec::with_capacity(estimated_tris);

        // Buffers for rings to avoid re-allocating them constantly
        // We only need the previous ring and the current ring to connect them
        let mut prev_ring: Vec<Vec3> = Vec::with_capacity(num_sides);
        let mut curr_ring: Vec<Vec3> = Vec::with_capacity(num_sides);
        
        // We'll also need to store the start ring for the start cap
        let mut start_ring: Vec<Vec3> = Vec::with_capacity(num_sides);

        for i in 0..node_count {
            let current = path.nodes[i].pos;
            let width = path.nodes[i].width;
            let height = path.nodes[i].height;
            
            // Calculate direction vectors to determine orientation
            let dir: Vec3;
            
            if i == 0 {
                dir = (path.nodes[i+1].pos - current).normalize();
            } else if i == node_count - 1 {
                dir = (current - path.nodes[i-1].pos).normalize();
            } else {
                let dir_prev = (current - path.nodes[i-1].pos).normalize();
                let dir_next = (path.nodes[i+1].pos - current).normalize();
                
                let sum = dir_prev + dir_next;
                if sum.length_squared() < config.direction_epsilon {
                    dir = dir_prev;
                } else {
                    dir = sum.normalize();
                }
            }
            
            // Define Frame
            let up = Vec3::Z;
            let safe_up = if dir.abs_diff_eq(Vec3::Z, config.up_vector_epsilon) || dir.abs_diff_eq(Vec3::NEG_Z, config.up_vector_epsilon) {
                Vec3::X
            } else {
                up
            };
            
            let miter_normal = dir;
            let right = miter_normal.cross(safe_up).normalize();
            let real_up = right.cross(miter_normal).normalize();
            
            // Generate current ring
            curr_ring.clear();
            for &(c, s) in &unit_circle {
                // Ellipse points on the 2D plane defined by right/real_up
                let local_x = c * width * 0.5;
                let local_y = s * height * 0.5;
                
                let point = current + right * local_x + real_up * local_y;
                curr_ring.push(point);
            }
            
            // If this is the first ring, save it for the cap
            if i == 0 {
                start_ring.clear();
                start_ring.extend_from_slice(&curr_ring);
            }
            
            // If we have a previous ring, stitch it to the current one
            if i > 0 {
                // Stitch rings
                for j in 0..num_sides {
                    let next_j = (j + 1) % num_sides;
                    
                    let p1 = prev_ring[j];
                    let p2 = curr_ring[j];
                    let p3 = curr_ring[next_j];
                    let p4 = prev_ring[next_j];
                    
                    add_quad(&mut local_triangles, p1, p2, p3, p4);
                }
            }
            
            // Swap buffers: curr becomes prev for next iteration
            // We can't simple swap because we need curr to be potentially reused or cleared.
            // But actually we can just swap the vectors themselves
            std::mem::swap(&mut prev_ring, &mut curr_ring);
        }
        
        // After loop, prev_ring actually holds the *last* ring generated (because of the swap at end of loop)
        let end_ring = &prev_ring;
        
        // Caps
        // Start Cap (i=0) - flip order to face out
        let center_start = path.nodes[0].pos;
        for j in 0..num_sides {
            let next_j = (j + 1) % num_sides;
            local_triangles.push(Triangle {
                v1: center_start,
                v2: start_ring[next_j],
                v3: start_ring[j],
                normal: (start_ring[j] - center_start).normalize(),
            });
        }
        
        // End Cap - normal order
        let center_end = path.nodes.last().unwrap().pos;
        for j in 0..num_sides {
            let next_j = (j + 1) % num_sides;
             local_triangles.push(Triangle {
                v1: center_end,
                v2: end_ring[j],
                v3: end_ring[next_j],
                normal: (end_ring[j] - center_end).normalize(),
            });
        }

        local_triangles
    };

    #[cfg(not(target_arch = "wasm32"))]
    {
        paths.par_iter()
             .flat_map(process_path)
             .collect()
    }

    #[cfg(target_arch = "wasm32")]
    {
        paths.iter()
             .flat_map(process_path)
             .collect()
    }
}

#[inline(always)]
fn add_quad(triangles: &mut Vec<Triangle>, p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) {
    // p1, p2, p3, p4 in CCW order
    // Tri 1: p1, p2, p3
    // Tri 2: p1, p3, p4
    
    // Compute face normal
    let normal = (p2 - p1).cross(p3 - p1).normalize_or_zero();
    if normal == Vec3::ZERO { return; }

    triangles.push(Triangle { v1: p1, v2: p2, v3: p3, normal });
    triangles.push(Triangle { v1: p1, v2: p3, v3: p4, normal });
}
