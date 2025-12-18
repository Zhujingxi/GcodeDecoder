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

pub fn generate_mesh(paths: &[ExtrusionPath]) -> Vec<Triangle> {
    let num_sides = 8; // Resolution of the tube (matches user edit)
    
    // Precompute unit circle (Cos, Sin) table
    let unit_circle: Vec<(f32, f32)> = (0..num_sides)
        .map(|j| {
            let angle = (j as f32 / num_sides as f32) * std::f32::consts::TAU;
            (angle.cos(), angle.sin())
        })
        .collect();

    // Helper closure to process one path
    let process_path = |path: &ExtrusionPath| -> Vec<Triangle> {
        let mut local_triangles = Vec::new();
        if path.nodes.len() < 2 {
            return local_triangles;
        }

        // We generate a ring of vertices for each node in the path.
        let mut rings: Vec<Vec<Vec3>> = Vec::with_capacity(path.nodes.len());
        
        for i in 0..path.nodes.len() {
            let current = path.nodes[i].pos;
            let width = path.nodes[i].width;
            let height = path.nodes[i].height;
            
            // Calculate direction vectors to determine orientation
            let dir: Vec3;
            
            if i == 0 {
                dir = (path.nodes[i+1].pos - current).normalize();
            } else if i == path.nodes.len() - 1 {
                dir = (current - path.nodes[i-1].pos).normalize();
            } else {
                let dir_prev = (current - path.nodes[i-1].pos).normalize();
                let dir_next = (path.nodes[i+1].pos - current).normalize();
                
                let sum = dir_prev + dir_next;
                if sum.length_squared() < 0.001 {
                    dir = dir_prev;
                } else {
                    dir = sum.normalize();
                }
            }
            
            // Define Frame
            let up = Vec3::Z;
            let safe_up = if dir.abs_diff_eq(Vec3::Z, 0.01) || dir.abs_diff_eq(Vec3::NEG_Z, 0.01) {
                Vec3::X
            } else {
                up
            };
            
            let miter_normal = dir;
            let right = miter_normal.cross(safe_up).normalize();
            let real_up = right.cross(miter_normal).normalize();
            
            let mut ring = Vec::with_capacity(num_sides);
            
            for &(c, s) in &unit_circle {
                // Ellipse points on the 2D plane defined by right/real_up
                let local_x = c * width * 0.5;
                let local_y = s * height * 0.5;
                
                let point = current + right * local_x + real_up * local_y;
                ring.push(point);
            }
            rings.push(ring);
        }
        
        // Stitch rings
        for i in 0..rings.len() - 1 {
            let ring_a = &rings[i];
            let ring_b = &rings[i+1];
            
            for j in 0..num_sides {
                let next_j = (j + 1) % num_sides;
                
                let p1 = ring_a[j];
                let p2 = ring_b[j];
                let p3 = ring_b[next_j];
                let p4 = ring_a[next_j];
                
                add_quad(&mut local_triangles, p1, p2, p3, p4);
            }
        }
        
        // Caps
        // Start Cap (i=0) - flip order to face out
        let start_ring = &rings[0];
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
        let end_ring = rings.last().unwrap();
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
