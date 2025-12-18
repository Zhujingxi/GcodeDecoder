/// Configuration for G-code processing and mesh generation
#[derive(Debug, Clone)]
pub struct Config {
    // === Basic Printer Parameters (user provides) ===
    /// Nozzle diameter in mm (e.g., 0.4)
    pub nozzle_diameter: f32,
    /// Default layer height in mm (e.g., 0.2)
    pub layer_height: f32,
    /// Filament diameter in mm (e.g., 1.75)
    pub filament_diameter: f32,

    // === Quality Settings ===
    /// Number of sides for tube mesh geometry
    pub mesh_sides: usize,

    // === Derived Values (calculated from basic params) ===
    /// Minimum detectable layer height
    min_layer_height: f32,
    /// Maximum detectable layer height
    max_layer_height: f32,
    /// Voxel size for optimizer
    voxel_size: f32,
    /// General epsilon for comparisons
    epsilon: f32,

    // === Advanced Optimizer Settings (with sensible defaults) ===
    /// Factor to pad the voxel grid bounding box (multiplier of voxel size). Default 3.0
    padding_factor: f32,
    /// Maximum dimension for voxel grid (to prevent OOM). Default 2000
    max_grid_dim: u32,
    /// Factor for conservative voxel coverage (multiplier of voxel size). Default 0.71 (~sqrt(2)/2)
    coverage_factor: f32,
    /// Multiplier for visibility check radius. Default 1.0
    visibility_radius_mult: f32,
    /// Modifier for voxel size relative to nozzle diameter. Default 0.5 (finer resolution)
    voxel_size_modifier: f32,

    // === Advanced Geometry Settings ===
    /// Epsilon for direction comparisons. Default 0.001
    direction_epsilon: f32,
    /// Epsilon for up-vector comparisons. Default 0.01
    up_vector_epsilon: f32,
}

impl Default for Config {
    fn default() -> Self {
        Self::new(0.4, 0.2, 1.75, 8)
    }
}

impl Config {
    /// Create a new config with basic parameters
    /// Derived values are calculated automatically
    pub fn new(
        nozzle_diameter: f32,
        layer_height: f32,
        filament_diameter: f32,
        mesh_sides: usize,
    ) -> Self {
        // Calculate derived values from basic parameters
        let min_layer_height = nozzle_diameter * 0.1;  // 10% of nozzle
        let max_layer_height = nozzle_diameter * 0.75; // 75% of nozzle
        
        // Defaults for improved optimizer
        let voxel_size_modifier = 0.5;
        let voxel_size = nozzle_diameter * voxel_size_modifier;              
        let epsilon = nozzle_diameter * 0.025;         // 2.5% of nozzle for comparisons

        Self {
            nozzle_diameter,
            layer_height,
            filament_diameter,
            mesh_sides,
            min_layer_height,
            max_layer_height,
            voxel_size,
            epsilon,
            // Defaults for advanced settings
            padding_factor: 3.0,
            max_grid_dim: 2000,
            coverage_factor: 0.71,
            visibility_radius_mult: 1.0,
            voxel_size_modifier,
            direction_epsilon: 0.001,
            up_vector_epsilon: 0.01,
        }
    }

    // === Getters for derived values ===
    
    /// Minimum detectable layer height (derived)
    pub fn min_layer_height(&self) -> f32 {
        self.min_layer_height
    }

    /// Maximum detectable layer height (derived)
    pub fn max_layer_height(&self) -> f32 {
        self.max_layer_height
    }

    /// Voxel size for optimizer (derived)
    pub fn voxel_size(&self) -> f32 {
        self.voxel_size
    }

    /// Epsilon for direction comparisons (derived)
    pub fn epsilon(&self) -> f32 {
        self.epsilon
    }

    /// Number of sides for tube geometry
    pub fn num_sides(&self) -> usize {
        self.mesh_sides
    }
    
    // === Getters for advanced settings ===
    pub fn padding_factor(&self) -> f32 { self.padding_factor }
    pub fn max_grid_dim(&self) -> u32 { self.max_grid_dim }
    pub fn coverage_factor(&self) -> f32 { self.coverage_factor }
    pub fn visibility_radius_mult(&self) -> f32 { self.visibility_radius_mult }
    pub fn voxel_size_modifier(&self) -> f32 { self.voxel_size_modifier }
    pub fn direction_epsilon(&self) -> f32 { self.direction_epsilon }
    pub fn up_vector_epsilon(&self) -> f32 { self.up_vector_epsilon }
}

/// Builder for Config with convenient defaults
#[derive(Debug, Clone)]
pub struct ConfigBuilder {
    nozzle_diameter: f32,
    layer_height: f32,
    filament_diameter: f32,
    mesh_sides: usize,
    // Optional overrides for advanced settings
    padding_factor: Option<f32>,
    max_grid_dim: Option<u32>,
    coverage_factor: Option<f32>,
    visibility_radius_mult: Option<f32>,
    voxel_size_modifier: Option<f32>,
}

impl Default for ConfigBuilder {
    fn default() -> Self {
        Self {
            nozzle_diameter: 0.4,
            layer_height: 0.2,
            filament_diameter: 1.75,
            mesh_sides: 32,
            padding_factor: None,
            max_grid_dim: None,
            coverage_factor: None,
            visibility_radius_mult: Some(1.3),
            voxel_size_modifier: None,
        }
    }
}

impl ConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn nozzle_diameter(mut self, value: f32) -> Self {
        self.nozzle_diameter = value;
        self
    }

    pub fn layer_height(mut self, value: f32) -> Self {
        self.layer_height = value;
        self
    }

    pub fn filament_diameter(mut self, value: f32) -> Self {
        self.filament_diameter = value;
        self
    }

    pub fn mesh_sides(mut self, sides: usize) -> Self {
        self.mesh_sides = sides;
        self
    }
    
    pub fn padding_factor(mut self, value: f32) -> Self {
        self.padding_factor = Some(value);
        self
    }

    pub fn build(self) -> Config {
        let mut config = Config::new(
            self.nozzle_diameter,
            self.layer_height,
            self.filament_diameter,
            self.mesh_sides,
        );
        
        if let Some(v) = self.padding_factor { config.padding_factor = v; }
        if let Some(v) = self.max_grid_dim { config.max_grid_dim = v; }
        if let Some(v) = self.coverage_factor { config.coverage_factor = v; }
        if let Some(v) = self.visibility_radius_mult { config.visibility_radius_mult = v; }
        
        // If voxel_size_modifier is changed, we need to recalculate voxel_size
        if let Some(v) = self.voxel_size_modifier { 
            config.voxel_size_modifier = v; 
            config.voxel_size = config.nozzle_diameter * v;
        }
        
        config
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = Config::default();
        assert_eq!(config.nozzle_diameter, 0.4);
        assert_eq!(config.layer_height, 0.2);
        assert_eq!(config.filament_diameter, 1.75);
        assert_eq!(config.num_sides(), 8);
        assert_eq!(config.padding_factor(), 3.0);
    }

    #[test]
    fn test_derived_values() {
        let config = Config::new(0.8, 0.3, 1.75, 16);
        
        // Check derived values scale with nozzle diameter
        assert!((config.min_layer_height() - 0.08).abs() < 0.001);  // 0.8 * 0.1
        assert!((config.max_layer_height() - 0.6).abs() < 0.001);   // 0.8 * 0.75
        assert!((config.voxel_size() - 0.4).abs() < 0.001); // 0.8 * 0.5 default modifier
        assert_eq!(config.num_sides(), 16);
    }

    #[test]
    fn test_builder() {
        let config = ConfigBuilder::new()
            .nozzle_diameter(0.6)
            .layer_height(0.25)
            .mesh_sides(4)
            .padding_factor(5.0)
            .build();
        
        assert_eq!(config.nozzle_diameter, 0.6);
        assert_eq!(config.layer_height, 0.25);
        assert_eq!(config.num_sides(), 4);
        assert_eq!(config.padding_factor(), 5.0);
    }
}
