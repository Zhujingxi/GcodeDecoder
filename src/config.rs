/// Configuration for G-code processing and mesh generation.
///
/// # Usage
///
/// ```rust
/// use gcodedecoder::config::Config;
/// // Quick setup with defaults
/// let config = Config::default();
///
/// // Custom basic parameters
/// let config = Config::new(0.4, 0.2, 1.75, 16);
///
/// // Full control via builder
/// let config = Config::builder()
///     .nozzle_diameter(0.6)
///     .layer_height(0.3)
///     .mesh_sides(32)
///     .build();
/// ```
#[derive(Debug, Clone)]
pub struct Config {
    // ─────────────────────────────────────────────────────────────
    // Printer Parameters
    // ─────────────────────────────────────────────────────────────
    /// Nozzle diameter in mm (default: 0.4)
    pub nozzle_diameter: f32,
    /// Default layer height in mm (default: 0.2)
    pub layer_height: f32,
    /// Filament diameter in mm (default: 1.75)
    pub filament_diameter: f32,

    // ─────────────────────────────────────────────────────────────
    // Mesh Quality
    // ─────────────────────────────────────────────────────────────
    /// Number of sides for tube mesh geometry (default: 16)
    pub mesh_sides: usize,

    // ─────────────────────────────────────────────────────────────
    // Optimizer Settings
    // ─────────────────────────────────────────────────────────────
    /// Voxel size modifier relative to nozzle diameter (default: 0.5)
    pub voxel_size_modifier: f32,
    /// Bounding box padding factor (multiplier of voxel size, default: 3.0)
    pub padding_factor: f32,
    /// Maximum grid dimension to prevent OOM (default: 2000)
    pub max_grid_dim: u32,
    /// Conservative voxel coverage factor (default: 0.71 ≈ √2/2)
    pub coverage_factor: f32,
    /// Visibility check radius multiplier (default: 1.3)
    pub visibility_radius_mult: f32,

    // ─────────────────────────────────────────────────────────────
    // Geometry Precision
    // ─────────────────────────────────────────────────────────────
    /// Epsilon for direction vector comparisons (default: 0.001)
    pub direction_epsilon: f32,
    /// Epsilon for up-vector comparisons (default: 0.01)
    pub up_vector_epsilon: f32,
}

// ═══════════════════════════════════════════════════════════════════
// Default values as constants for easy reference
// ═══════════════════════════════════════════════════════════════════

impl Config {
    pub const DEFAULT_NOZZLE_DIAMETER: f32 = 0.4;
    pub const DEFAULT_LAYER_HEIGHT: f32 = 0.2;
    pub const DEFAULT_FILAMENT_DIAMETER: f32 = 1.75;
    pub const DEFAULT_MESH_SIDES: usize = 8;
    pub const DEFAULT_VOXEL_SIZE_MODIFIER: f32 = 0.5;
    pub const DEFAULT_PADDING_FACTOR: f32 = 3.0;
    pub const DEFAULT_MAX_GRID_DIM: u32 = 2000;
    pub const DEFAULT_COVERAGE_FACTOR: f32 = 0.71;
    pub const DEFAULT_VISIBILITY_RADIUS_MULT: f32 = 1.3;
    pub const DEFAULT_DIRECTION_EPSILON: f32 = 0.001;
    pub const DEFAULT_UP_VECTOR_EPSILON: f32 = 0.01;
}

impl Default for Config {
    fn default() -> Self {
        Self {
            nozzle_diameter: Self::DEFAULT_NOZZLE_DIAMETER,
            layer_height: Self::DEFAULT_LAYER_HEIGHT,
            filament_diameter: Self::DEFAULT_FILAMENT_DIAMETER,
            mesh_sides: Self::DEFAULT_MESH_SIDES,
            voxel_size_modifier: Self::DEFAULT_VOXEL_SIZE_MODIFIER,
            padding_factor: Self::DEFAULT_PADDING_FACTOR,
            max_grid_dim: Self::DEFAULT_MAX_GRID_DIM,
            coverage_factor: Self::DEFAULT_COVERAGE_FACTOR,
            visibility_radius_mult: Self::DEFAULT_VISIBILITY_RADIUS_MULT,
            direction_epsilon: Self::DEFAULT_DIRECTION_EPSILON,
            up_vector_epsilon: Self::DEFAULT_UP_VECTOR_EPSILON,
        }
    }
}

impl Config {
    /// Create config with basic printer parameters (uses defaults for everything else)
    pub fn new(
        nozzle_diameter: f32,
        layer_height: f32,
        filament_diameter: f32,
        mesh_sides: usize,
    ) -> Self {
        Self {
            nozzle_diameter,
            layer_height,
            filament_diameter,
            mesh_sides,
            ..Default::default()
        }
    }

    /// Start building a config with the builder pattern
    pub fn builder() -> ConfigBuilder {
        ConfigBuilder::default()
    }

    // ─────────────────────────────────────────────────────────────
    // Computed Values (derived from other fields)
    // ─────────────────────────────────────────────────────────────

    /// Computed voxel size for optimizer (nozzle_diameter × voxel_size_modifier)
    #[inline]
    pub fn voxel_size(&self) -> f32 {
        self.nozzle_diameter * self.voxel_size_modifier
    }

    /// Minimum detectable layer height (10% of nozzle diameter)
    #[inline]
    pub fn min_layer_height(&self) -> f32 {
        self.nozzle_diameter * 0.1
    }

    /// Maximum detectable layer height (75% of nozzle diameter)
    #[inline]
    pub fn max_layer_height(&self) -> f32 {
        self.nozzle_diameter * 0.75
    }

    /// General epsilon for float comparisons (2.5% of nozzle diameter)
    #[inline]
    pub fn epsilon(&self) -> f32 {
        self.nozzle_diameter * 0.025
    }
}

// ═══════════════════════════════════════════════════════════════════
// Builder Pattern
// ═══════════════════════════════════════════════════════════════════

/// Fluent builder for `Config`
#[derive(Debug, Clone)]
pub struct ConfigBuilder(Config);

impl Default for ConfigBuilder {
    fn default() -> Self {
        Self(Config::default())
    }
}

impl ConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    // ─────────────────────────────────────────────────────────────
    // Printer Parameters
    // ─────────────────────────────────────────────────────────────

    pub fn nozzle_diameter(mut self, value: f32) -> Self {
        self.0.nozzle_diameter = value;
        self
    }

    pub fn layer_height(mut self, value: f32) -> Self {
        self.0.layer_height = value;
        self
    }

    pub fn filament_diameter(mut self, value: f32) -> Self {
        self.0.filament_diameter = value;
        self
    }

    // ─────────────────────────────────────────────────────────────
    // Mesh Quality
    // ─────────────────────────────────────────────────────────────

    pub fn mesh_sides(mut self, value: usize) -> Self {
        self.0.mesh_sides = value;
        self
    }

    // ─────────────────────────────────────────────────────────────
    // Optimizer Settings
    // ─────────────────────────────────────────────────────────────

    pub fn voxel_size_modifier(mut self, value: f32) -> Self {
        self.0.voxel_size_modifier = value;
        self
    }

    pub fn padding_factor(mut self, value: f32) -> Self {
        self.0.padding_factor = value;
        self
    }

    pub fn max_grid_dim(mut self, value: u32) -> Self {
        self.0.max_grid_dim = value;
        self
    }

    pub fn coverage_factor(mut self, value: f32) -> Self {
        self.0.coverage_factor = value;
        self
    }

    pub fn visibility_radius_mult(mut self, value: f32) -> Self {
        self.0.visibility_radius_mult = value;
        self
    }

    // ─────────────────────────────────────────────────────────────
    // Geometry Precision
    // ─────────────────────────────────────────────────────────────

    pub fn direction_epsilon(mut self, value: f32) -> Self {
        self.0.direction_epsilon = value;
        self
    }

    pub fn up_vector_epsilon(mut self, value: f32) -> Self {
        self.0.up_vector_epsilon = value;
        self
    }

    // ─────────────────────────────────────────────────────────────
    // Build
    // ─────────────────────────────────────────────────────────────

    /// Consume the builder and return the configured `Config`
    pub fn build(self) -> Config {
        self.0
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
        assert_eq!(config.mesh_sides, Config::DEFAULT_MESH_SIDES);
        assert_eq!(config.padding_factor, Config::DEFAULT_PADDING_FACTOR);
    }

    #[test]
    fn test_derived_values() {
        let config = Config::new(0.8, 0.3, 1.75, 16);
        
        // Check computed values scale with nozzle diameter
        assert!((config.min_layer_height() - 0.08).abs() < 0.001);  // 0.8 * 0.1
        assert!((config.max_layer_height() - 0.6).abs() < 0.001);   // 0.8 * 0.75
        assert!((config.voxel_size() - 0.4).abs() < 0.001); // 0.8 * 0.5 default modifier
        assert_eq!(config.mesh_sides, 16);
    }

    #[test]
    fn test_builder() {
        let config = Config::builder()
            .nozzle_diameter(0.6)
            .layer_height(0.25)
            .mesh_sides(4)
            .padding_factor(5.0)
            .build();
        
        assert_eq!(config.nozzle_diameter, 0.6);
        assert_eq!(config.layer_height, 0.25);
        assert_eq!(config.mesh_sides, 4);
        assert_eq!(config.padding_factor, 5.0);
    }
}
