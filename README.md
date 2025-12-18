# GcodeDecoder

A high-performance G-code to STL converter written in Rust. Converts 3D printer G-code files into accurate 3D mesh representations of the extruded filament paths.

## Features

- **Accurate Geometry** - Creates tube meshes representing actual filament extrusion paths
- **Volumetric Calculation** - Dynamically calculates line width from extrusion amount
- **Mesh Optimization** - Optional internal path culling for smaller files
- **WASM Support** - Can be compiled to WebAssembly for browser use
- **Fully Configurable** - All parameters adapt to your printer setup

## Installation

```bash
cargo build --release
```

## Usage

### Basic

```bash
gcodedecoder input.gcode                    # Outputs input.stl
gcodedecoder input.gcode -o output.stl      # Custom output path
gcodedecoder input.gcode --optimize         # Remove internal paths
```

### Configuration Options

| Option | Default | Description |
|--------|---------|-------------|
| `--nozzle-diameter` | 0.4 | Nozzle diameter in mm |
| `--layer-height` | 0.2 | Layer height in mm |
| `--filament-diameter` | 1.75 | Filament diameter in mm |
| `--mesh-sides` | 8 | Tube resolution (4=low, 8=medium, 16=high) |
| `--optimize` | off | Remove internal/hidden paths |

### Examples

```bash
# Standard 0.4mm nozzle
gcodedecoder Benchy.gcode --optimize

# Large 0.8mm nozzle with thick layers
gcodedecoder input.gcode --nozzle-diameter 0.8 --layer-height 0.3

# High quality mesh (16 sides per tube)
gcodedecoder input.gcode --mesh-sides 16

# 2.85mm filament (Prusa, Ultimaker)
gcodedecoder input.gcode --filament-diameter 2.85
```

## WASM API

```javascript
// Default configuration
const stlBuffer = convert_gcode_to_stl(gcodeString);

// Custom configuration
const stlBuffer = convert_gcode_to_stl_with_config(
  gcodeString,
  0.4,   // nozzle_diameter
  0.2,   // layer_height  
  1.75,  // filament_diameter
  8      // mesh_sides
);
```

## How It Works

1. **Parse** - G-code commands are parsed into move instructions
2. **Process** - Extrusion moves are converted to path segments with width/height
3. **Generate** - Tube mesh geometry is created around each path
4. **Optimize** (optional) - Internal paths are culled via voxel flood-fill
5. **Export** - Binary STL is written

## Derived Parameters

The following values are automatically calculated from nozzle diameter:

| Parameter | Formula | Example (0.4mm nozzle) |
|-----------|---------|------------------------|
| Min layer height | nozzle × 0.1 | 0.04mm |
| Max layer height | nozzle × 0.75 | 0.3mm |
| Voxel size | nozzle | 0.4mm |

## License

MIT
