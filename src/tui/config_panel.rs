use crate::config::Config;

#[derive(Debug, Clone)]
pub struct ConfigPanel {
    pub config: Config,
    pub original_config: Config,
}

impl ConfigPanel {
    pub fn new(config: Config) -> Self {
        Self {
            original_config: config.clone(),
            config,
        }
    }

    pub fn reset(&mut self) {
        self.config = self.original_config.clone();
    }

    pub fn reset_to_defaults(&mut self) {
        self.config = Config::default();
    }

    pub fn increment_nozzle(&mut self, big: bool) {
        let delta = if big { 0.1 } else { 0.05 };
        self.config.nozzle_diameter = (self.config.nozzle_diameter + delta).min(2.0);
    }

    pub fn decrement_nozzle(&mut self, big: bool) {
        let delta = if big { 0.1 } else { 0.05 };
        self.config.nozzle_diameter = (self.config.nozzle_diameter - delta).max(0.1);
    }

    pub fn increment_layer(&mut self, big: bool) {
        let delta = if big { 0.1 } else { 0.05 };
        self.config.layer_height = (self.config.layer_height + delta).min(1.0);
    }

    pub fn decrement_layer(&mut self, big: bool) {
        let delta = if big { 0.1 } else { 0.05 };
        self.config.layer_height = (self.config.layer_height - delta).max(0.05);
    }

    pub fn increment_filament(&mut self, big: bool) {
        let delta = if big { 0.5 } else { 0.05 };
        self.config.filament_diameter = (self.config.filament_diameter + delta).min(3.0);
    }

    pub fn decrement_filament(&mut self, big: bool) {
        let delta = if big { 0.5 } else { 0.05 };
        self.config.filament_diameter = (self.config.filament_diameter - delta).max(1.0);
    }

    pub fn increment_sides(&mut self, big: bool) {
        let delta = if big { 4 } else { 1 };
        self.config.mesh_sides = (self.config.mesh_sides + delta).min(32);
    }

    pub fn decrement_sides(&mut self, big: bool) {
        let delta = if big { 4 } else { 1 };
        self.config.mesh_sides = self.config.mesh_sides.saturating_sub(delta).max(3);
    }

    pub fn toggle_remove_internal(&mut self) {
        self.config.remove_internal = !self.config.remove_internal;
    }

    pub fn has_changes(&self) -> bool {
        self.config.nozzle_diameter != self.original_config.nozzle_diameter
            || self.config.layer_height != self.original_config.layer_height
            || self.config.filament_diameter != self.original_config.filament_diameter
            || self.config.mesh_sides != self.original_config.mesh_sides
            || self.config.remove_internal != self.original_config.remove_internal
    }

    pub fn validation_warning(&self) -> Option<String> {
        if self.config.nozzle_diameter < self.config.layer_height {
            Some(format!(
                "Warning: Nozzle ({:.2}) < Layer height ({:.2})",
                self.config.nozzle_diameter, self.config.layer_height
            ))
        } else {
            None
        }
    }
}
