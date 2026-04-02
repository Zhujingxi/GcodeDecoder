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
    }

    #[test]
    fn test_derived_values() {
        let config = Config::new(0.8, 0.3, 1.75, 16);

        // Check computed values scale with nozzle diameter
        assert!((config.min_layer_height() - 0.08).abs() < 0.001); // 0.8 * 0.1
        assert!((config.max_layer_height() - 0.6).abs() < 0.001); // 0.8 * 0.75
        assert_eq!(config.mesh_sides, 16);
    }

    #[test]
    fn test_builder() {
        let config = Config::builder()
            .nozzle_diameter(0.6)
            .layer_height(0.25)
            .mesh_sides(4)
            .build();

        assert_eq!(config.nozzle_diameter, 0.6);
        assert_eq!(config.layer_height, 0.25);
        assert_eq!(config.mesh_sides, 4);
    }
}
