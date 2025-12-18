# GcodeDecoder

GcodeDecoder is a high-performance Rust tool that converts 3D printer G-code into 3D mesh geometry (STL). It effectively visualizes the extruded filament path, making it useful for verifying slice data, identifying potential print issues, or creating realistic previews.

## Features

- **High Fidelity**: Generates accurate "pillar-like" geometry for every extrusion segment, handling corners and line widths correctly based on volumetric physics.
- **Optimization**: Optional internal path culling to remove geometry that isn't visible from the outside. This significantly reduces mesh complexity and file size while maintaining visual accuracy.
- **Performance**: Built with Rust and `rayon` for multi-threaded processing, ensuring fast conversion speeds even for large G-code files.
- **WASM Support**: Designed to compile to WebAssembly for client-side web usage.
- **Configurable**: Fully customizable printer parameters including nozzle diameter, layer height, filament diameter, and mesh resolution.

## Installation

Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed.

1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/gcodedecoder.git
    cd gcodedecoder
    ```

2.  Build the project:
    ```bash
    cargo build --release
    ```

## Usage

### CLI

Run the tool using `cargo run`. The basic syntax is:

```bash
cargo run --release -- <input_file> [options]
```

#### Examples

**Basic Conversion:**
Convert `test.gcode` to `test.stl` (default output name):
```bash
cargo run --release -- ./test.gcode
```

**Specify Output File:**
```bash
cargo run --release -- ./test.gcode --output ./output_model.stl
```

**With Optimization (Remove Internal Paths):**
This will analyze the mesh and remove voxels/paths that are completely occluded by outer shells.
```bash
cargo run --release -- ./test.gcode --optimize
```

**Custom Printer Settings:**
Adjust parameters to match your slicer settings for the most accurate result:
```bash
cargo run --release -- ./test.gcode \
  --nozzle-diameter 0.6 \
  --layer-height 0.3 \
  --filament-diameter 1.75 \
  --mesh-sides 16
```

### Configuration Options

| Flag | Description | Default |
|------|-------------|---------|
| `--output`, `-o` | Output STL file path | `<input>.stl` |
| `--optimize` | Enable removal of internal paths | Disabled |
| `--nozzle-diameter` | Printer nozzle diameter (mm) | 0.4 |
| `--layer-height` | Layer height (mm) | 0.2 |
| `--filament-diameter` | Filament diameter (mm) | 1.75 |
| `--mesh-sides` | Resolution of the tube mesh (4-16 recommended) | 8 |

## WebAssembly (WASM)

This project can be compiled to WASM for use in web applications.

1.  Install `wasm-pack`: https://rustwasm.github.io/wasm-pack/installer/
2.  Build for web:
    ```bash
    wasm-pack build --target web
    ```
3.  The output will be in the `pkg` directory, ready to be imported into your JavaScript app.

## Development

### Benchmarking

To test the performance of the optimization algorithms:

```bash
cargo run --release --bin benchmark
```


## License

[MIT License](LICENSE)
