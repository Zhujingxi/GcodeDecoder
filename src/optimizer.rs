use crate::config::Config;
use crate::processor::ExtrusionPath;
use glam::{Vec3, UVec3};
use std::sync::atomic::{AtomicU64, Ordering};
use rayon::prelude::*;
use std::time::Instant;

// Voxel State Constants
const VOXEL_EMPTY: u8 = 0;
const VOXEL_OCCUPIED: u8 = 1;
const VOXEL_EXTERNAL: u8 = 2;
const VOXEL_VISIBLE_REGION: u8 = 3;

// Bit Packing Constants
const BITS_PER_VOXEL: usize = 2;
const VOXELS_PER_WORD: usize = 64 / BITS_PER_VOXEL;
const BIT_MASK: u64 = 0b11;


struct VoxelGrid {
    voxels: Vec<AtomicU64>,
    dimensions: UVec3,
    min_bound: Vec3,
    voxel_size: f32,
    coverage_factor: f32,
}

impl VoxelGrid {
    fn new(min_bound: Vec3, max_bound: Vec3, voxel_size: f32, config: &Config) -> Self {
        let start = Instant::now();
        
        // Add padding to ensure external flood fill can reach around
        let padding = voxel_size * config.padding_factor;
        let padded_min = min_bound - Vec3::splat(padding);
        let size = (max_bound - min_bound) + Vec3::splat(padding * 2.0);
        
        let dims = (size / voxel_size).ceil().as_uvec3();
        // Clamp dimensions to avoid extreme memory usage
        let dims = dims.min(UVec3::splat(config.max_grid_dim));
        
        // Use usize for total calculation to avoid u32 overflow
        let total_voxels = dims.x as usize * dims.y as usize * dims.z as usize;
        let num_words = (total_voxels + VOXELS_PER_WORD - 1) / VOXELS_PER_WORD;

        // Parallel Initialization
        let voxels: Vec<AtomicU64> = (0..num_words)
            .into_par_iter()
            .map(|_| AtomicU64::new(0))
            .collect();

        println!("  - Grid Init: {:.2?} ({} voxels)", start.elapsed(), total_voxels);

        Self {
            voxels,
            dimensions: dims,
            min_bound: padded_min,
            voxel_size,
            coverage_factor: config.coverage_factor,
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

    #[inline(always)]
    fn get_index_unchecked(&self, x: u32, y: u32, z: u32) -> usize {
        (z as usize * self.dimensions.y as usize * self.dimensions.x as usize) + 
        (y as usize * self.dimensions.x as usize) + 
        (x as usize)
    }

    #[inline(always)]
    fn get_index_from_coord(&self, coord: UVec3) -> usize {
        self.get_index_unchecked(coord.x, coord.y, coord.z)
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
        let word_idx = idx / VOXELS_PER_WORD;
        let bit_offset = (idx % VOXELS_PER_WORD) * BITS_PER_VOXEL;
        let word = self.voxels[word_idx].load(Ordering::Relaxed);
        ((word >> bit_offset) & BIT_MASK) as u8
    }
    
    // Returns true if the state was updated from EMPTY to TARGET
    fn try_set_state(&self, idx: usize, target: u8) -> bool {
        let word_idx = idx / VOXELS_PER_WORD;
        let bit_offset = (idx % VOXELS_PER_WORD) * BITS_PER_VOXEL;
        let mask = BIT_MASK << bit_offset;
        let target_bits = (target as u64) << bit_offset;
        
        let atom = &self.voxels[word_idx];
        let mut old = atom.load(Ordering::Relaxed);
        
        loop {
            let current_val = (old >> bit_offset) & BIT_MASK;
            if current_val != VOXEL_EMPTY as u64 {
                return false; // Already occupied or external
            }
            
            let new = (old & !mask) | target_bits;
            match atom.compare_exchange_weak(old, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return true,
                Err(x) => old = x,
            }
        }
    }
    
    #[inline]
    fn set_occupied(&self, idx: usize) {
        let word_idx = idx / VOXELS_PER_WORD;
        let bit_offset = (idx % VOXELS_PER_WORD) * BITS_PER_VOXEL;
        let bits = (VOXEL_OCCUPIED as u64) << bit_offset;
        // fetch_or is atomic. We blindly set the bit. 
        // Since VOXEL_OCCUPIED is 1, and EMPTY is 0, this transitions 0->1.
        // It might conflict with EXTERNAL (2) which is 10 binary, resulting in 11 (3 -> Visible Region).
        // BUT mark_occupied happens BEFORE flood fill, so only state is 0 or 1.
        self.voxels[word_idx].fetch_or(bits, Ordering::Relaxed);
    }
    
    // Safe setter that updates 2 (External) to 3 (Visible) or 0 to 3.
    // Used in Dilation pass.
    fn set_visible_region(&self, idx: usize) {
         let word_idx = idx / VOXELS_PER_WORD;
        let bit_offset = (idx % VOXELS_PER_WORD) * BITS_PER_VOXEL;
        let mask = BIT_MASK << bit_offset;
        let target_bits = (VOXEL_VISIBLE_REGION as u64) << bit_offset;
        
        let atom = &self.voxels[word_idx];
        let mut old = atom.load(Ordering::Relaxed);
        
        loop {
            // If already visible region, nothing to do
            let current_val = (old >> bit_offset) & BIT_MASK;
             if current_val == VOXEL_VISIBLE_REGION as u64 {
                return;
            }
            
            // If it's occupied (1), do NOT overwrite it?
            // "Visible Region" means "Air near path".
            // If it's occupied, it's the path itself.
            // We want to mark EMPTY or EXTERNAL as VISIBLE_REGION?
            // Actually, we usually want to know if a path node is *inside* a Visible Region.
            // If the voxel at path node is Occupied, that's fine. 
            // We need to check if it's *near* External.
            // So we are marking the Air.
            // If we overwrite Occupied, we lose the path info?
            // Actually, `filter_paths` only checks if the path's location is in a "Visible" zone.
            // If we mark the path voxels themselves as Visible Region, that's fine too.
            // But let's assume we proceed:
            
            let new = (old & !mask) | target_bits;
            match atom.compare_exchange_weak(old, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return,
                Err(x) => old = x,
            }
        }
    }

    /// Mark occupied voxels - uses volumetric rasterization
    /// Parallelized using &self and Rayon
    fn mark_occupied(&self, paths: &[ExtrusionPath]) {
        let start = Instant::now();
        let dim_x = self.dimensions.x as usize;
        let dim_y = self.dimensions.y as usize;
        let stride_z = dim_x * dim_y;
        
        // Use chunks to reduce the overhead of parallel scheduling for small tasks
        // Sorting paths beforehand (spatial locality) makes this extremely cache efficient
        paths.par_chunks(64).for_each(|chunk| {
            for path in chunk {
                if path.nodes.len() < 2 { continue; }

                for i in 0..path.nodes.len() - 1 {
                    let p1 = path.nodes[i].pos;
                    let p2 = path.nodes[i+1].pos;
                    
                    let width = path.nodes[i].width.max(path.nodes[i+1].width);
                    let height = path.nodes[i].height.max(path.nodes[i+1].height);
                    
                    // Radius logic
                    let radius = width * 0.5 + self.voxel_size * self.coverage_factor; 
                    let vertical_radius = height * 0.5 + self.voxel_size * self.coverage_factor;
                    
                    // Precompute segment math
                    // let p1_xy = p1; // Unused
                    
                    let min_p = p1.min(p2) - Vec3::new(radius, radius, vertical_radius);
                    let max_p = p1.max(p2) + Vec3::new(radius, radius, vertical_radius);

                    // Grid bounds
                    let start_idx = match self.world_to_grid(min_p) {
                         Some(idx) => idx, None => UVec3::ZERO
                    };
                    let end_idx = match self.world_to_grid(max_p) {
                        Some(idx) => idx, None => self.dimensions - UVec3::ONE
                    };

                    // Optimizing the loop:
                    // For each Z and Y, we want to find the range of X that satisfies distance check.
                    
                    // Segment vector
                    let seg = p2 - p1;
                    let seg_len_sq = seg.length_squared();
                    let inv_seg_len_sq = if seg_len_sq > 0.0 { 1.0 / seg_len_sq } else { 0.0 };

                    for z in start_idx.z..=end_idx.z {
                        let voxel_z = self.min_bound.z + z as f32 * self.voxel_size + self.voxel_size * 0.5;
                        let base_z = (z as usize) * stride_z;
                        
                        // Z distance check (simple AABB prune first, but precise later)
                        // Actually, distance to segment is 3D.
                        // Let's stick to the separation of XY and Z for "Prism with rounded ends" shape?
                        // The original code treated it as specialized extrusion shape: 
                        // "Capsule in XY" vs "Rectangle in Z"? 
                        // Original code: dist_xy_sq <= radius_sq AND dist_z <= vertical_radius.
                        // This approximates the flattened extruded line.
                        
                        for y in start_idx.y..=end_idx.y {
                            let voxel_y = self.min_bound.y + y as f32 * self.voxel_size + self.voxel_size * 0.5;
                            let base_y = base_z + (y as usize) * dim_x;
                            
                            // To Optimize X:
                            // We need `dist_xy_sq( (x, y), segment_xy ) <= radius_sq`.
                            // This defines a range of X values.
                            // Solving for X is hard because `t` (closest point param) depends on X.
                            // But we can check bounds.
                            
                            // Fallback to "Check every X in box" but optimized math
                            for x in start_idx.x..=end_idx.x {
                                let idx = base_y + (x as usize);
                                
                                // Peek first
                                if self.get_voxel(idx) == VOXEL_OCCUPIED { continue; }
                                
                                let voxel_x = self.min_bound.x + x as f32 * self.voxel_size + self.voxel_size * 0.5;
                                let voxel_pos = Vec3::new(voxel_x, voxel_y, voxel_z);
                                
                                // Distance to finite line segment
                                let t = (voxel_pos - p1).dot(seg) * inv_seg_len_sq;
                                let t_clamped = t.clamp(0.0, 1.0);
                                let closest = p1 + seg * t_clamped;
                                
                                // We are using separate radius for Z and XY in original code.
                                // Let's replicate logic:
                                let closest_xy = Vec3::new(closest.x, closest.y, 0.0);
                                let voxel_xy = Vec3::new(voxel_x, voxel_y, 0.0);
                                let dist_xy_sq = voxel_xy.distance_squared(closest_xy);
                                
                                if dist_xy_sq <= radius * radius {
                                    let dist_z = (voxel_z - closest.z).abs();
                                    if dist_z <= vertical_radius {
                                        self.set_occupied(idx);
                                    }
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
        
        // Estimate frontier size
        let estimated_size = (self.dimensions.x * self.dimensions.y * 2) as usize;
        let mut frontier = Vec::with_capacity(estimated_size);
        
        // 1. Seed External from Faces
        // Using explicit loop for seeds is fast enough
        let dim = self.dimensions;
        
        // Function to process seed
        let check_and_push = |idx: usize, list: &mut Vec<usize>| {
            if self.try_set_state(idx, VOXEL_EXTERNAL) {
                list.push(idx);
            }
        };

        // Gather seeds in parallel stripes
        let seeds: Vec<usize> = (0..dim.z).into_par_iter().map(|z| {
            let mut local_seeds = Vec::new();
            for y in 0..dim.y {
                check_and_push(self.get_index_unchecked(0, y, z), &mut local_seeds);
                check_and_push(self.get_index_unchecked(dim.x-1, y, z), &mut local_seeds);
            }
            for x in 1..dim.x-1 {
                check_and_push(self.get_index_unchecked(x, 0, z), &mut local_seeds);
                check_and_push(self.get_index_unchecked(x, dim.y-1, z), &mut local_seeds);
            }
            local_seeds
        }).flatten().collect();
        
        frontier.extend(seeds);
        // Add Top/Bottom faces
        for y in 0..dim.y {
            for x in 0..dim.x {
                check_and_push(self.get_index_unchecked(x, y, 0), &mut frontier);
                check_and_push(self.get_index_unchecked(x, y, dim.z-1), &mut frontier);
            }
        }

        println!("    > Initial seeds: {}", frontier.len());
        
        // BFS
        let stride_x = 1usize;
        let stride_y = self.dimensions.x as usize;
        let stride_z = stride_y * self.dimensions.y as usize;
        
        while !frontier.is_empty() {
            let next_frontier: Vec<usize> = frontier.par_iter().map(|&idx| {
                 let mut local = Vec::new();
                 // Optimize: Don't decode coord if possible. Check bounds using remainder? 
                 // Cost of coord_from_index is division, which is slow.
                 // We can simply check neighbors if safe.
                 
                 let coord = self.coord_from_index(idx); // Still safer for boundary checks
                 
                 if coord.x > 0 { 
                     let n = idx - stride_x; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 if coord.x < dim.x - 1 { 
                     let n = idx + stride_x; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 
                 if coord.y > 0 { 
                     let n = idx - stride_y; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 if coord.y < dim.y - 1 { 
                     let n = idx + stride_y; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 
                 if coord.z > 0 { 
                     let n = idx - stride_z; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 if coord.z < dim.z - 1 { 
                     let n = idx + stride_z; 
                     if self.try_set_state(n, VOXEL_EXTERNAL) { local.push(n); } 
                 }
                 
                 local
            }).flatten().collect();
            
            frontier = next_frontier;
        }

        println!("  - Flood Fill: {:.2?}", start.elapsed());
    }
    
    // Dilate the External Region to create "Visible Region"
    // This effectively marks any voxel within Radius of Air as Visible.
    fn dilate_visible_region(&self, radius_voxels: i32) {
        if radius_voxels <= 0 { return; }
        let start = Instant::now();
        
        let dim = self.dimensions;
        let stride_y = dim.x as usize;
        let stride_z = stride_y * dim.y as usize;

        // Naive dilation is expensive (iterate all voxels, check neighbors).
        // Better: Iterate only EXTERNAL voxels and mark their neighbors.
        // Even better: Use the BFS approach again! 
        // Seed with all External voxels.
        // Run BFS for K steps.
        
        // But we already have the External voxels.
        // And we don't store "Frontier of External".
        // Re-scanning grid to find External frontier is fast.
        
        // 1. Find all External voxels that border Non-External (Occupied or Empty)
        // Actually, just find ALL External voxels (state=2) and use as seeds.
        // Since we want to mark Empty/Occupied as Visible.
        
        // Let's assume we just want to expand "External".
        // Any voxel reachable from External within K steps => Visible Region.
        
        // Collect all External voxels? That's huge.
        // But we only need the *boundary* of External region.
        // We can scan array.
        
        let frontier: Vec<usize> = (0..self.voxels.len() * VOXELS_PER_WORD).into_par_iter()
            .filter_map(|idx| {
                if self.get_voxel(idx) == VOXEL_EXTERNAL {
                    Some(idx)
                } else {
                    None
                }
            }).collect();
            
        // Optimization: Use a bitset for "Visited" or just rely on state change?
        // We want to change value to 3 (Visible). 
        // External (2) is already "Visible". We can leave it as 2.
        // We only change 0 or 1 to 3.
        
        // If we start with all External, we can BFS.
        // Depth-limited BFS.
        let mut next_frontier = frontier;
        
        for _step in 0..radius_voxels {
            if next_frontier.is_empty() { break; }
            
            next_frontier = next_frontier.par_iter().map(|&idx| {
                 let coord = self.coord_from_index(idx);
                 let mut local = Vec::new();
                 
                 let neighbors = [
                    if coord.x > 0 { Some(idx - 1) } else { None },
                    if coord.x < dim.x - 1 { Some(idx + 1) } else { None },
                    if coord.y > 0 { Some(idx - stride_y) } else { None },
                    if coord.y < dim.y - 1 { Some(idx + stride_y) } else { None },
                    if coord.z > 0 { Some(idx - stride_z) } else { None },
                    if coord.z < dim.z - 1 { Some(idx + stride_z) } else { None },
                 ];
                 
                 for n_opt in neighbors {
                     if let Some(n_idx) = n_opt {
                         let val = self.get_voxel(n_idx);
                         if val != VOXEL_EXTERNAL && val != VOXEL_VISIBLE_REGION {
                             // Mark this as visible
                             // We use a CAS loop similar to try_set_state but for Visible
                             // We don't care about return value, just try.
                             self.set_visible_region(n_idx);
                             local.push(n_idx);
                         }
                     }
                 }
                 local
            }).flatten().collect();
        }
        
        println!("  - Dilation: {:.2?}", start.elapsed());
    }
}

pub fn optimize_paths(mut paths: Vec<ExtrusionPath>, config: &Config) -> Vec<ExtrusionPath> {
    if paths.is_empty() {
        return paths;
    }

    // 0. Spatial Sorting to improve cache locality
    // Z-curve or just coordinate sort. Simple X/Y/Z sort is very effective.
    let start_sort = Instant::now();
    paths.par_sort_unstable_by(|a, b| {
        let p1 = if a.nodes.is_empty() { Vec3::ZERO } else { a.nodes[0].pos };
        let p2 = if b.nodes.is_empty() { Vec3::ZERO } else { b.nodes[0].pos };
        
        // Sort by Z (layer), then Y, then X
        p1.z.partial_cmp(&p2.z).unwrap_or(std::cmp::Ordering::Equal)
            .then(p1.y.partial_cmp(&p2.y).unwrap_or(std::cmp::Ordering::Equal))
            .then(p1.x.partial_cmp(&p2.x).unwrap_or(std::cmp::Ordering::Equal))
    });
    println!("  - Sort Paths: {:.2?}", start_sort.elapsed());

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
    let voxel_size = config.voxel_size();
    let grid = VoxelGrid::new(min_bound, max_bound, voxel_size, config);

    // 3. Mark Occupied (Parallel)
    grid.mark_occupied(&paths);

    // 4. Flood Fill External (Parallel BFS)
    grid.flood_fill_external();
    
    // 4b. Dilate External Region
    // Calculate dilation radius in voxels
    // We want to cover "visibility_radius_mult * width"
    // Approximate average width or use conservative max
    let avg_width = config.nozzle_diameter; // Approximation
    let check_radius = (avg_width * config.visibility_radius_mult / voxel_size).ceil() as i32;
    grid.dilate_visible_region(check_radius);

    // 5. Filter Paths (Parallel)
    let start_filter = Instant::now();
    // Keep path if ANY node is in VOXEL_VISIBLE_REGION or VOXEL_EXTERNAL
    let result = paths.into_par_iter().filter(|path| {
        // Fast check: Check sampled points
        for node in &path.nodes {
            if let Some(idx) = grid.world_to_grid(node.pos).map(|uv| grid.get_index_from_coord(uv)) {
                let v = grid.get_voxel(idx);
                if v == VOXEL_EXTERNAL || v == VOXEL_VISIBLE_REGION {
                    return true;
                }
            }
        }
        
        // Check midpoints for finer granularity
        if path.nodes.len() > 1 {
            for i in 0..path.nodes.len() - 1 {
                let midpoint = (path.nodes[i].pos + path.nodes[i+1].pos) * 0.5;
                if let Some(idx) = grid.world_to_grid(midpoint).map(|uv| grid.get_index_from_coord(uv)) {
                    let v = grid.get_voxel(idx);
                    if v == VOXEL_EXTERNAL || v == VOXEL_VISIBLE_REGION {
                        return true;
                    }
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
        let config = Config::default();
        let mut paths = Vec::new();
        let gap = 0.4;
        let size = 5;
        
        for x in 0..size {
            for y in 0..size {
                paths.push(ExtrusionPath {
                    nodes: vec![
                        PathNode { pos: Vec3::new(x as f32 * gap, y as f32 * gap, 0.0), width: 0.4, height: 0.2 },
                        PathNode { pos: Vec3::new(x as f32 * gap + gap, y as f32 * gap, 0.0), width: 0.4, height: 0.2 },
                    ]
                });
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
        
        assert_eq!(final_count, initial_count);
    }
}
