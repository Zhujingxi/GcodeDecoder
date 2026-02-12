# GcodeDecoder TUI Design

## Overview

An interactive terminal user interface for GcodeDecoder that allows users to adjust conversion parameters visually before processing G-code files into STL meshes.

## Architecture

The TUI uses **ratatui** with a component-based architecture:

**Core Components:**
- **App**: Main state container, handles event loop and mode switching
- **ConfigPanel**: Editable form for nozzle diameter, layer height, filament diameter, mesh sides, and "Remove Internal" toggle
- **PreviewPanel**: Shows input file info, segment count, and layer count
- **StatusBar**: Current operation status and keyboard shortcuts

**State Management:**
- App holds `AppState` enum: `Editing`, `Processing`, `Complete`, `Error`
- Config changes trigger real-time validation (e.g., nozzle < layer height warning)
- Processing runs in a background thread with progress updates via channels

**Navigation:**
- Tab/Shift+Tab cycles through form fields
- Arrow keys adjust numeric values (with Ctrl for larger increments)
- Enter starts conversion, 'q' or Esc to quit

## Screen Layout

```
┌─────────────────────────────────────────────────────────────┐
│ GcodeDecoder TUI - Benchy.gcode                    [q]uit  │
├─────────────────────────────────────────────────────────────┤
│ ┌─────────────────────┐  ┌───────────────────────────────┐ │
│ │ Configuration       │  │ Preview                       │ │
│ │                     │  │                               │ │
│ │ Nozzle: 0.40 mm  [▲▼]│  │ Input: Benchy.gcode          │ │
│ │ Layer:  0.20 mm  [▲▼]│  │ Size: 2.4 MB                 │ │
│ │ Filament: 1.75 mm[▲▼]│  │                              │ │
│ │ Sides:  8        [▲▼]│  │ Segments: 15,432             │ │
│ │ Remove Internal: [✓]│  │ Layers: 247                  │ │
│ │                     │  │                              │ │
│ │ [Defaults] [Revert] │  │                              │ │
│ │                     │  │                              │ │
│ └─────────────────────┘  │                              │ │
│                          │                              │ │
│                          │                              │ │
│                          │                              │ │
│                          │                              │ │
│                          │ [Convert]                    │ │
│                          └───────────────────────────────┘ │
├─────────────────────────────────────────────────────────────┤
│ Status: Ready | ↑↓ Adjust | Tab Next | Enter Convert       │
└─────────────────────────────────────────────────────────────┘
```

**Layout Details:**
- Left panel: 35% width for config form with visual controls
- Right panel: 65% width showing input file info
- Bottom bar: Context-aware help text
- Colors: Blue for headers, green for success, yellow for warnings, red for errors

## Data Flow

1. **Startup**: CLI parses args → if `--tui` flag present, launch TUI mode with input file path
2. **Initialize**: Load file, parse G-code to get segment count and layer count
3. **Edit Loop**: User adjusts parameters → validate in real-time
4. **Convert**: Spawn blocking task on thread pool → stream progress updates via `mpsc` channel
5. **Complete**: Show success message with output path

## Error Handling

- **Parse errors**: Show in status bar with red highlighting
- **Validation warnings**: Yellow indicator on problematic field with tooltip
- **Conversion errors**: Modal popup with error details, return to config screen
- **Interrupt handling**: Ctrl+C gracefully exits

## Keyboard Shortcuts

- `Tab` / `Shift+Tab`: Navigate fields
- `↑` / `↓`: Adjust values (hold `Ctrl` for ×10)
- `Enter`: Start conversion
- `r`: Reset to defaults
- `q` / `Esc`: Quit (confirm if unsaved changes)

## Implementation Structure

```
src/
├── main.rs           # CLI entry, --tui flag handling
├── lib.rs            # Existing library code
├── tui/
│   ├── mod.rs        # TUI module exports
│   ├── app.rs        # App struct, event loop, state management
│   ├── ui.rs         # Widget rendering functions
│   ├── config_panel.rs  # Config form component
│   ├── preview_panel.rs # File info display
│   └── event.rs      # Input handling, key mappings
```

## Dependencies

Add to `Cargo.toml`:
```toml
ratatui = "0.29"
crossterm = "0.28"
```

## Key Types

```rust
enum AppState {
    Editing,
    Processing { progress: f32 },
    Complete { output_path: PathBuf },
    Error { message: String },
}

struct App {
    input_path: PathBuf,
    config: Config,
    state: AppState,
    focused_field: Field,
}
```

## Testing Strategy

- Unit tests for validation logic
- Mock input events for navigation testing
- Integration test with test file

## Usage

```bash
# Launch TUI mode
cargo run --release -- --tui ./test.gcode

# Or as default when no output specified
cargo run --release -- ./test.gcode --tui
```
