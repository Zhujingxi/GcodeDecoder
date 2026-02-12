use crate::tui::FileInfo;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct PreviewPanel {
    pub file_info: Option<FileInfo>,
}

impl PreviewPanel {
    pub fn new() -> Self {
        Self { file_info: None }
    }

    pub fn load_file(&mut self, path: &Path) -> anyhow::Result<()> {
        let metadata = fs::metadata(path)?;
        let size_bytes = metadata.len();

        let file_info = FileInfo {
            path: path.to_path_buf(),
            size_bytes,
            segment_count: 0,
            layer_count: 0,
        };

        self.file_info = Some(file_info);
        Ok(())
    }

    pub fn update_stats(&mut self, segment_count: usize, layer_count: usize) {
        if let Some(ref mut info) = self.file_info {
            info.segment_count = segment_count;
            info.layer_count = layer_count;
        }
    }

    pub fn format_file_size(bytes: u64) -> String {
        const UNITS: &[&str] = &["B", "KB", "MB", "GB"];
        let mut size = bytes as f64;
        let mut unit_idx = 0;

        while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
            size /= 1024.0;
            unit_idx += 1;
        }

        format!("{:.1} {}", size, UNITS[unit_idx])
    }

    pub fn file_name(&self) -> String {
        self.file_info
            .as_ref()
            .and_then(|info| info.path.file_name())
            .map(|name| name.to_string_lossy().to_string())
            .unwrap_or_else(|| "Unknown".to_string())
    }

    pub fn file_size(&self) -> String {
        self.file_info
            .as_ref()
            .map(|info| Self::format_file_size(info.size_bytes))
            .unwrap_or_else(|| "Unknown".to_string())
    }

    pub fn segment_count(&self) -> String {
        self.file_info
            .as_ref()
            .map(|info| format!("{}", info.segment_count))
            .unwrap_or_else(|| "-".to_string())
    }

    pub fn layer_count(&self) -> String {
        self.file_info
            .as_ref()
            .map(|info| format!("{}", info.layer_count))
            .unwrap_or_else(|| "-".to_string())
    }
}

impl Default for PreviewPanel {
    fn default() -> Self {
        Self::new()
    }
}
