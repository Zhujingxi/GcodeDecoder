pub mod app;
pub mod config_panel;
pub mod event;
pub mod preview_panel;
pub mod ui;

use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Field {
    NozzleDiameter,
    LayerHeight,
    FilamentDiameter,
    MeshSides,
    RemoveInternal,
}

impl Field {
    pub fn next(self) -> Self {
        match self {
            Field::NozzleDiameter => Field::LayerHeight,
            Field::LayerHeight => Field::FilamentDiameter,
            Field::FilamentDiameter => Field::MeshSides,
            Field::MeshSides => Field::RemoveInternal,
            Field::RemoveInternal => Field::NozzleDiameter,
        }
    }

    pub fn prev(self) -> Self {
        match self {
            Field::NozzleDiameter => Field::RemoveInternal,
            Field::LayerHeight => Field::NozzleDiameter,
            Field::FilamentDiameter => Field::LayerHeight,
            Field::MeshSides => Field::FilamentDiameter,
            Field::RemoveInternal => Field::MeshSides,
        }
    }
}

#[derive(Debug, Clone)]
pub enum AppState {
    Editing,
    Processing { progress: f32 },
    Complete { output_path: PathBuf },
    Error { message: String },
}

#[derive(Debug, Clone)]
pub struct FileInfo {
    pub path: PathBuf,
    pub size_bytes: u64,
    pub segment_count: usize,
    pub layer_count: usize,
}
