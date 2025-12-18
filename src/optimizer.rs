use crate::config::Config;
use crate::processor::ExtrusionPath;
use glam::{Vec3, UVec3};
use std::sync::atomic::{AtomicU8, Ordering};
use rayon::prelude::*;
use std::time::Instant;

// Voxel State Constants
const VOXEL_EMPTY: u8 = 0;
const VOXEL_OCCUPIED: u8 = 1;
const VOXEL_EXTERNAL: u8 = 2;

struct VoxelGrid {
    voxels: Vec<AtomicU8>,
    dimensions: UVec3,
    min_bound: Vec3,
    voxel_size: f32,
    coverage_factor: f32,
    visibility_radius_mult: f32,
}

impl VoxelGrid {
    fn new(min_bound: Vec3, max_bound: Vec3, voxel_size: f32, config: &Config) -> Self {
        let start = Instant::now();
        
        // Add padding to ensure external flood fill can reach around
        let padding = voxel_size * config.padding_factor();
        let padded_min = min_bound - Vec3::splat(padding);
        let size = (max_bound - min_bound) + Vec3::splat(padding * 2.0);
        
        let dims = (size / voxel_size).ceil().as_uvec3();
        // Clamp dimensions to avoid extreme memory usage
        let dims = dims.min(UVec3::splat(config.max_grid_dim()));
        
        // Use usize for total calculation to avoid u32 overflow
        let total_voxels = dims.x as usize * dims.y as usize * dims.z as usize;

        // Parallel Initialization
        let voxels: Vec<AtomicU8> = (0..total_voxels)
            .into_par_iter()
            .map(|_| AtomicU8::new(VOXEL_EMPTY))
            .collect();

        println!("  - Grid Init: {:.2?} ({} voxels)", start.elapsed(), total_voxels);

        Self {
            voxels,
            dimensions: dims,
            min_bound: padded_min,
            voxel_size,
            coverage_factor: config.coverage_factor(),
            visibility_radius_mult: config.visibility_radius_mult(),
        }
    }

    fn world_to_grid(&self, pos: Vec3) -> Option<UVec3> {
        if pos.cmplt(self.min_bound).any() {
            return None;
        }
        let offset = pos - self.min_bound;
        let idx = (offset / self.voxel_size).floor().as_uvec3();
        
        if idx.cmplt(self.dimensions).all() {
            Some(idx)
        } else {
            None
        }
    }

    fn get_index(&self, coord: UVec3) -> usize {
        (coord.z as usize * self.dimensions.y as usize * self.dimensions.x as usize) + 
        (coord.y as usize * self.dimensions.x as usize) + 
        (coord.x as usize)
    }
    
    fn coord_from_index(&self, idx: usize) -> UVec3 {
        let stride_y = self.dimensions.x as usize;
        let stride_z = stride_y * self.dimensions.y as usize;
        
        let z = idx / stride_z;
        let rem = idx % stride_z;
        let y = rem / stride_y;
        let x = rem % stride_y;
        
        UVec3::new(x as u32, y as u32, z as u32)
    }

    fn get_voxel(&self, idx: usize) -> u8 {
        self.voxels[idx].load(Ordering::Relaxed)
    }
    
    // Returns true if the state was updated from EMPTY to TARGET
    fn try_set_state(&self, idx: usize, target: u8) -> bool {
        self.voxels[idx].compare_exchange(
            VOXEL_EMPTY, 
            target, 
            Ordering::Relaxed, 
            Ordering::Relaxed
        ).is_ok()
    }
    
    fn set_occupied(&self, idx: usize) {
        self.voxels[idx].store(VOXEL_OCCUPIED, Ordering::Relaxed);
    }

    /// Mark occupied voxels - uses volumetric rasterization
    /// Parallelized using &self and Rayon
    fn mark_occupied(&self, paths: &[ExtrusionPath]) {
        let start = Instant::now();
        let dim_x = self.dimensions.x as usize;
        let dim_y = self.dimensions.y as usize;
        let stride_z = dim_x * dim_y;
        
        paths.par_iter().for_each(|path| {
            if path.nodes.len() < 2 { return; }

            for i in 0..path.nodes.len() - 1 {
                let p1 = path.nodes[i].pos;
                let p2 = path.nodes[i+1].pos;
                
                // Use the max width of the segment for radius
                let width = path.nodes[i].width.max(path.nodes[i+1].width);
                let height = path.nodes[i].height.max(path.nodes[i+1].height);
                
                // Radius for rasterization: width/2 + conservative voxel coverage
                let radius = width * 0.5 + self.voxel_size * self.coverage_factor; 
                let vertical_radius = height * 0.5 + self.voxel_size * self.coverage_factor;
                let radius_sq = radius * radius;
                
                // Bounding box of the segment expanded by radius
                let seg_min = p1.min(p2) - Vec3::new(radius, radius, vertical_radius);
                let seg_max = p1.max(p2) + Vec3::new(radius, radius, vertical_radius);
                
                // Convert to grid bounds
                let start_idx = match self.world_to_grid(seg_min) {
                    Some(idx) => idx,
                    None => UVec3::ZERO, 
                };
                
                let end_idx = match self.world_to_grid(seg_max) {
                    Some(idx) => idx,
                    None => self.dimensions - UVec3::ONE, 
                };

                // Precompute XY projection of the segment
                let p1_xy = Vec3::new(p1.x, p1.y, 0.0);
                let p2_xy = Vec3::new(p2.x, p2.y, 0.0);
                let segment_vec = p2_xy - p1_xy;
                let l2 = segment_vec.length_squared();
                let inv_l2 = if l2 == 0.0 { 0.0 } else { 1.0 / l2 };

                // Iterate over the bounding box
                for z in start_idx.z..=end_idx.z {
                    let voxel_z = self.min_bound.z + z as f32 * self.voxel_size + self.voxel_size * 0.5;
                    let base_idx_z = (z as usize) * stride_z;
                    
                    for y in start_idx.y..=end_idx.y {
                        let voxel_y = self.min_bound.y + y as f32 * self.voxel_size + self.voxel_size * 0.5;
                        let base_idx_y = base_idx_z + (y as usize) * dim_x;
                        
                        for x in start_idx.x..=end_idx.x {
                            let idx = base_idx_y + (x as usize);
                            
                            // Optimization: Check before writing to reduce cache coherence traffic
                            if self.get_voxel(idx) == VOXEL_OCCUPIED {
                                continue;
                            }

                            // Calculate voxel center X
                            let voxel_x = self.min_bound.x + x as f32 * self.voxel_size + self.voxel_size * 0.5;
                            let center_xy = Vec3::new(voxel_x, voxel_y, 0.0);
                            
                            // Distance check
                            let t = (center_xy - p1_xy).dot(segment_vec) * inv_l2;
                            let t_clamped = t.clamp(0.0, 1.0);
                            let closest_xy = p1_xy + segment_vec * t_clamped;
                            
                            let dist_xy_sq = center_xy.distance_squared(closest_xy);
                            
                            if dist_xy_sq <= radius_sq {
                                // Check Z height
                                let closest_z = p1.z + (p2.z - p1.z) * t_clamped;
                                let dist_z = (voxel_z - closest_z).abs();
                                
                                if dist_z <= vertical_radius {
                                    self.set_occupied(idx);
                                }
                            }
                        }
                    }
                }
            }
        });
        println!("  - Mark Occupied: {:.2?}", start.elapsed());
    }

    fn flood_fill_external(&self) {
        let start = Instant::now();
        
        // Pre-allocate frontier with estimated capacity for all 6 faces
        let estimated_frontier_size = ((self.dimensions.x * self.dimensions.y * 2 + 
                                        self.dimensions.y * self.dimensions.z * 2 + 
                                        self.dimensions.x * self.dimensions.z * 2) as usize).min(1_000_000);
        let mut frontier = Vec::with_capacity(estimated_frontier_size);
        
        // Parallel seeding: Process each face's rows in parallel
        // X-faces (YZ planes at x=0 and x=max)
        let x_seeds: Vec<usize> = (0..self.dimensions.z).into_par_iter().flat_map(|z| {
            let mut local = Vec::new();
            for y in 0..self.dimensions.y {
                let idx0 = self.get_index(UVec3::new(0, y, z));
                let idx1 = self.get_index(UVec3::new(self.dimensions.x - 1, y, z));
                if self.try_set_state(idx0, VOXEL_EXTERNAL) { local.push(idx0); }
                if self.try_set_state(idx1, VOXEL_EXTERNAL) { local.push(idx1); }
            }
            local
        }).collect();
        
        // Y-faces (XZ planes at y=0 and y=max)
        let y_seeds: Vec<usize> = (0..self.dimensions.z).into_par_iter().flat_map(|z| {
            let mut local = Vec::new();
            for x in 0..self.dimensions.x {
                let idx0 = self.get_index(UVec3::new(x, 0, z));
                let idx1 = self.get_index(UVec3::new(x, self.dimensions.y - 1, z));
                if self.try_set_state(idx0, VOXEL_EXTERNAL) { local.push(idx0); }
                if self.try_set_state(idx1, VOXEL_EXTERNAL) { local.push(idx1); }
            }
            local
        }).collect();
        
        // Z-faces (XY planes at z=0 and z=max)
        let z_seeds: Vec<usize> = (0..self.dimensions.y).into_par_iter().flat_map(|y| {
            let mut local = Vec::new();
            for x in 0..self.dimensions.x {
                let idx0 = self.get_index(UVec3::new(x, y, 0));
                let idx1 = self.get_index(UVec3::new(x, y, self.dimensions.z - 1));
                if self.try_set_state(idx0, VOXEL_EXTERNAL) { local.push(idx0); }
                if self.try_set_state(idx1, VOXEL_EXTERNAL) { local.push(idx1); }
            }
            local
        }).collect();
        
        // Merge all seeds
        frontier.extend(x_seeds);
        frontier.extend(y_seeds);
        frontier.extend(z_seeds);
        
        println!("    > Initial seeds: {}", frontier.len());

        // Helper to get neighbor indices
        
        let stride_x = 1usize;
        let stride_y = self.dimensions.x as usize;
        let stride_z = (self.dimensions.x * self.dimensions.y) as usize;
        
        let dim_x = self.dimensions.x as i32;
        let dim_y = self.dimensions.y as i32;
        let dim_z = self.dimensions.z as i32;



        // Parallel Frontier BFS
        while !frontier.is_empty() {
            // Map current frontier to next frontier in parallel
            let next_frontier: Vec<usize> = frontier.par_iter().map(|&idx| {
                let coord = self.coord_from_index(idx); // cost of calc
                let mut local_next = Vec::new();
                
                let cur_x = coord.x as i32;
                let cur_y = coord.y as i32;
                let cur_z = coord.z as i32;

                // Check 6 neighbors
                // Manually unrolled for explicit checks
                
                // X-
                if cur_x > 0 {
                    let n_idx = idx - stride_x;
                    if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }
                // X+
                if cur_x < dim_x - 1 {
                    let n_idx = idx + stride_x;
                    if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }
                
                // Y-
                if cur_y > 0 {
                    let n_idx = idx - stride_y;
                     if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }
                // Y+
                if cur_y < dim_y - 1 {
                    let n_idx = idx + stride_y;
                     if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }
                
                // Z-
                if cur_z > 0 {
                    let n_idx = idx - stride_z;
                     if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }
                // Z+
                if cur_z < dim_z - 1 {
                    let n_idx = idx + stride_z;
                     if self.try_set_state(n_idx, VOXEL_EXTERNAL) { local_next.push(n_idx); }
                }

                local_next
            })
            .flatten()
            .collect();
            
            frontier = next_frontier;
        }

        println!("  - Flood Fill: {:.2?}", start.elapsed());
    }

    /// Check visibility with 2x extrusion width radius expansion
    /// A path is visible if any voxel within 2x width/height distance is External
    fn is_visible_with_radius(&self, pos: Vec3, width: f32, height: f32) -> bool {
        if let Some(coord) = self.world_to_grid(pos) {
            // Use 1x extrusion dimensions for visibility check radius
            let check_radius = (width.max(height) * self.visibility_radius_mult / self.voxel_size).ceil() as i32;
            let check_radius = check_radius.max(1); // At least 1 voxel
            
            let dim_x = self.dimensions.x as i32;
            let dim_y = self.dimensions.y as i32;
            let dim_z = self.dimensions.z as i32;

            let start_x = (coord.x as i32 - check_radius).max(0);
            let end_x = (coord.x as i32 + check_radius).min(dim_x - 1);
            
            let start_y = (coord.y as i32 - check_radius).max(0);
            let end_y = (coord.y as i32 + check_radius).min(dim_y - 1);
            
            let start_z = (coord.z as i32 - check_radius).max(0);
            let end_z = (coord.z as i32 + check_radius).min(dim_z - 1);

            let stride_y = self.dimensions.x as usize;
            let stride_z = (self.dimensions.x * self.dimensions.y) as usize;

            // SIMD-friendly row-based iteration with early exit
            for z in start_z..=end_z {
                let base_z = (z as usize) * stride_z;
                for y in start_y..=end_y {
                    let base_y = base_z + (y as usize) * stride_y;
                    let row_start = base_y + (start_x as usize);
                    let row_end = base_y + (end_x as usize) + 1;
                    
                    // Check entire row at once - enables vectorization
                    if self.voxels[row_start..row_end].iter().any(|v| v.load(Ordering::Relaxed) == VOXEL_EXTERNAL) {
                        return true;
                    }
                }
            }
        }
        false
    }
}

pub fn optimize_paths(paths: Vec<ExtrusionPath>, config: &Config) -> Vec<ExtrusionPath> {
    if paths.is_empty() {
        return paths;
    }

    // 1. Calculate Bounds (Parallel Reduction)
    let (min_bound, max_bound) = paths.par_iter()
        .map(|path| {
            path.nodes.iter().fold(
                (Vec3::splat(f32::MAX), Vec3::splat(f32::MIN)),
                |(min, max), node| (min.min(node.pos), max.max(node.pos))
            )
        })
        .reduce(
            || (Vec3::splat(f32::MAX), Vec3::splat(f32::MIN)),
            |(min1, max1), (min2, max2)| (min1.min(min2), max1.max(max2))
        );

    // 2. Init Voxel Grid with FINER resolution
    // Using half the nozzle diameter for better accuracy
    let voxel_size = config.voxel_size() * config.voxel_size_modifier();
    let grid = VoxelGrid::new(min_bound, max_bound, voxel_size, config);

    // 3. Mark Occupied (Parallel)
    grid.mark_occupied(&paths);

    // 4. Flood Fill External (Parallel BFS)
    grid.flood_fill_external();

    // 5. Filter Paths (Parallel)
    let start_filter = Instant::now();
    // Keep path if ANY node is within 2x extrusion width of External air
    let result = paths.into_par_iter().filter(|path| {
        // Check all nodes with 2x radius visibility
        for node in &path.nodes {
            if grid.is_visible_with_radius(node.pos, node.width, node.height) {
                return true; 
            }
        }
        
        // Also check segment midpoints for long segments
        if path.nodes.len() > 1 {
            for i in 0..path.nodes.len() - 1 {
                let midpoint = (path.nodes[i].pos + path.nodes[i+1].pos) * 0.5;
                let avg_width = (path.nodes[i].width + path.nodes[i+1].width) * 0.5;
                let avg_height = (path.nodes[i].height + path.nodes[i+1].height) * 0.5;
                if grid.is_visible_with_radius(midpoint, avg_width, avg_height) {
                    return true;
                }
            }
        }
        
        false
    }).collect();
    
    println!("  - Filter Paths: {:.2?}", start_filter.elapsed());
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::processor::PathNode;

    #[test]
    fn test_internal_culling() {
        // Create a dense block of paths that forms a solid interior
        let config = Config::default();
        let mut paths = Vec::new();
        let gap = 0.3;
        let grid_size = 15;
        
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
        
        let initial_count = paths.len();
        let optimized = optimize_paths(paths, &config);
        let final_count = optimized.len();
        
        println!("Optimized: {} -> {}", initial_count, final_count);
        
        assert!(final_count <= initial_count); 
        assert!(final_count > 0);
    }
    
    #[test]
    fn test_surface_paths_preserved() {
        // Create a hollow shell - all paths should be kept
        let config = Config::default();
        let mut paths = Vec::new();
        let gap = 0.4;
        let size = 5;
        
        // Only create paths on the surface faces
        for x in 0..size {
            for y in 0..size {
                // Bottom face (z=0)
                paths.push(ExtrusionPath {
                    nodes: vec![
                        PathNode { pos: Vec3::new(x as f32 * gap, y as f32 * gap, 0.0), width: 0.4, height: 0.2 },
                        PathNode { pos: Vec3::new(x as f32 * gap + gap, y as f32 * gap, 0.0), width: 0.4, height: 0.2 },
                    ]
                });
                // Top face
                paths.push(ExtrusionPath {
                    nodes: vec![
                        PathNode { pos: Vec3::new(x as f32 * gap, y as f32 * gap, size as f32 * gap), width: 0.4, height: 0.2 },
                        PathNode { pos: Vec3::new(x as f32 * gap + gap, y as f32 * gap, size as f32 * gap), width: 0.4, height: 0.2 },
                    ]
                });
            }
        }
        
        let initial_count = paths.len();
        let optimized = optimize_paths(paths, &config);
        let final_count = optimized.len();
        
        println!("Surface test: {} -> {}", initial_count, final_count);
        
        // All surface paths should be preserved
        assert_eq!(final_count, initial_count);
    }
}
