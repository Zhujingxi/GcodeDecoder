use eframe::egui;
use gcodedecoder::processor::ExtrusionPath;
use gcodedecoder::{config, export, geometry, parser, processor};
use std::sync::mpsc;
use std::time::Instant;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1100.0, 750.0])
            .with_min_inner_size([800.0, 500.0]),
        ..Default::default()
    };

    eframe::run_native(
        "GcodeDecoder",
        options,
        Box::new(|cc| {
            cc.egui_ctx.set_visuals(egui::Visuals::dark());
            Ok(Box::new(App::default()))
        }),
    )
}

// ---------------------------------------------------------------------------
// Thread messages
// ---------------------------------------------------------------------------

enum Msg {
    Log(String),
    Preview {
        lines: Vec<PreviewLine>,
        min_z: f32,
        max_z: f32,
        center: [f32; 3],
        bounds_min: [f32; 3],
        bounds_max: [f32; 3],
    },
    Done(String),
    Error(String),
}

// ---------------------------------------------------------------------------
// Preview data
// ---------------------------------------------------------------------------

struct PreviewLine {
    start: [f32; 3],
    end: [f32; 3],
    z: f32,
}

// ---------------------------------------------------------------------------
// Application state
// ---------------------------------------------------------------------------

struct App {
    // File paths
    input_path: String,
    output_path: String,

    // Config
    nozzle_diameter: f32,
    layer_height: f32,
    filament_diameter: f32,
    mesh_sides: i32,

    // Status
    log: Vec<String>,
    busy: bool,
    rx: Option<mpsc::Receiver<Msg>>,

    // Preview
    preview_lines: Vec<PreviewLine>,
    min_z: f32,
    max_z: f32,
    layer_cutoff: f32,
    bounds_min: [f32; 3],
    bounds_max: [f32; 3],

    // Camera
    yaw: f32,
    pitch: f32,
    zoom: f32,
    pan: egui::Vec2,
    center: [f32; 3],
}

impl Default for App {
    fn default() -> Self {
        Self {
            input_path: String::new(),
            output_path: String::new(),
            nozzle_diameter: 0.4,
            layer_height: 0.2,
            filament_diameter: 1.75,
            mesh_sides: 8,
            log: vec!["Ready.".into()],
            busy: false,
            rx: None,
            preview_lines: Vec::new(),
            min_z: 0.0,
            max_z: 1.0,
            layer_cutoff: f32::MAX,
            bounds_min: [0.0; 3],
            bounds_max: [0.0; 3],
            yaw: -std::f32::consts::FRAC_PI_4,
            pitch: std::f32::consts::FRAC_PI_6,
            zoom: 1.0,
            pan: egui::Vec2::ZERO,
            center: [0.0; 3],
        }
    }
}

// ---------------------------------------------------------------------------
// eframe::App
// ---------------------------------------------------------------------------

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_messages();

        if self.busy {
            ctx.request_repaint();
        }

        egui::SidePanel::left("controls")
            .resizable(true)
            .default_width(320.0)
            .min_width(280.0)
            .show(ctx, |ui| self.draw_controls(ui));

        egui::CentralPanel::default().show(ctx, |ui| self.draw_preview(ui));
    }
}

// ---------------------------------------------------------------------------
// Message polling
// ---------------------------------------------------------------------------

impl App {
    fn poll_messages(&mut self) {
        // Take the receiver out of self so we don't hold an immutable borrow
        // on self while processing messages that need &mut self.
        let Some(rx) = self.rx.take() else { return };
        let mut finished = false;

        loop {
            match rx.try_recv() {
                Ok(msg) => match msg {
                    Msg::Log(s) => self.log.push(s),
                    Msg::Preview {
                        lines,
                        min_z,
                        max_z,
                        center,
                        bounds_min,
                        bounds_max,
                    } => {
                        self.preview_lines = lines;
                        self.min_z = min_z;
                        self.max_z = max_z;
                        self.center = center;
                        self.bounds_min = bounds_min;
                        self.bounds_max = bounds_max;
                        self.layer_cutoff = max_z;
                        self.auto_zoom();
                    }
                    Msg::Done(s) => {
                        self.log.push(s);
                        finished = true;
                    }
                    Msg::Error(s) => {
                        self.log.push(format!("ERROR: {}", s));
                        finished = true;
                    }
                },
                Err(mpsc::TryRecvError::Empty) => break,
                Err(mpsc::TryRecvError::Disconnected) => {
                    if self.busy {
                        self.log
                            .push("ERROR: Worker thread terminated unexpectedly.".into());
                    }
                    finished = true;
                    break;
                }
            }
        }

        if finished {
            self.busy = false;
            // rx is dropped – don't put it back
        } else {
            self.rx = Some(rx);
        }
    }
}

// ---------------------------------------------------------------------------
// Side-panel controls
// ---------------------------------------------------------------------------

impl App {
    fn draw_controls(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.heading("GcodeDecoder");
            ui.separator();

            // --- Input file ---
            ui.label("Input G-code File");
            ui.horizontal(|ui| {
                let w = (ui.available_width() - 80.0).max(60.0);
                ui.add_sized(
                    [w, 20.0],
                    egui::TextEdit::singleline(&mut self.input_path)
                        .hint_text("Select a .gcode file…"),
                );
                if ui.button("Browse…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("G-code", &["gcode", "gco", "g"])
                        .pick_file()
                    {
                        self.input_path = path.display().to_string();
                        if self.output_path.is_empty() {
                            let mut out = path;
                            out.set_extension("stl");
                            self.output_path = out.display().to_string();
                        }
                    }
                }
            });

            ui.add_space(6.0);

            // --- Output file ---
            ui.label("Output STL File");
            ui.horizontal(|ui| {
                let w = (ui.available_width() - 80.0).max(60.0);
                ui.add_sized(
                    [w, 20.0],
                    egui::TextEdit::singleline(&mut self.output_path)
                        .hint_text("Output path…"),
                );
                if ui.button("Browse…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("STL", &["stl"])
                        .save_file()
                    {
                        self.output_path = path.display().to_string();
                    }
                }
            });

            ui.add_space(10.0);
            ui.separator();

            // --- Printer settings ---
            ui.heading("Printer Settings");
            ui.add_space(4.0);

            egui::Grid::new("printer_settings")
                .num_columns(2)
                .spacing([8.0, 6.0])
                .show(ui, |ui| {
                    ui.label("Nozzle (mm)");
                    ui.add(
                        egui::DragValue::new(&mut self.nozzle_diameter)
                            .range(0.1..=2.0)
                            .speed(0.01)
                            .fixed_decimals(2),
                    );
                    ui.end_row();

                    ui.label("Layer Height (mm)");
                    ui.add(
                        egui::DragValue::new(&mut self.layer_height)
                            .range(0.01..=1.0)
                            .speed(0.01)
                            .fixed_decimals(2),
                    );
                    ui.end_row();

                    ui.label("Filament (mm)");
                    ui.add(
                        egui::DragValue::new(&mut self.filament_diameter)
                            .range(0.5..=3.5)
                            .speed(0.01)
                            .fixed_decimals(2),
                    );
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.separator();

            // --- Mesh quality ---
            ui.heading("Mesh Quality");
            ui.add_space(4.0);
            ui.add(egui::Slider::new(&mut self.mesh_sides, 3..=32).text("sides"));

            ui.add_space(10.0);
            ui.separator();
            ui.add_space(6.0);

            // --- Action buttons ---
            let can_act = !self.busy && !self.input_path.is_empty();

            ui.horizontal(|ui| {
                let half = ui.available_width() / 2.0 - 4.0;
                ui.add_enabled_ui(can_act, |ui| {
                    if ui
                        .add_sized([half, 32.0], egui::Button::new("Preview"))
                        .clicked()
                    {
                        self.start_preview();
                    }
                });
                ui.add_enabled_ui(can_act, |ui| {
                    if ui
                        .add_sized([half, 32.0], egui::Button::new("Export STL"))
                        .clicked()
                    {
                        self.start_export();
                    }
                });
            });

            if self.busy {
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Working…");
                });
            }

            ui.add_space(10.0);
            ui.separator();

            // --- Layer filter (only when data loaded) ---
            if !self.preview_lines.is_empty() {
                ui.heading("View");
                ui.add_space(4.0);
                ui.add(
                    egui::Slider::new(&mut self.layer_cutoff, self.min_z..=self.max_z)
                        .text("mm"),
                );

                let visible = self
                    .preview_lines
                    .iter()
                    .filter(|l| l.z <= self.layer_cutoff)
                    .count();
                ui.label(format!(
                    "{} / {} segments",
                    visible,
                    self.preview_lines.len()
                ));

                if ui.button("Reset View").clicked() {
                    self.yaw = -std::f32::consts::FRAC_PI_4;
                    self.pitch = std::f32::consts::FRAC_PI_6;
                    self.pan = egui::Vec2::ZERO;
                    self.auto_zoom();
                }

                ui.add_space(8.0);
                ui.separator();
            }

            // --- Log ---
            ui.heading("Log");
            ui.add_space(4.0);
            egui::ScrollArea::vertical()
                .max_height(200.0)
                .stick_to_bottom(true)
                .id_salt("log_scroll")
                .show(ui, |ui| {
                    for msg in &self.log {
                        ui.label(msg);
                    }
                });
        });
    }
}

// ---------------------------------------------------------------------------
// 3-D preview
// ---------------------------------------------------------------------------

impl App {
    fn draw_preview(&mut self, ui: &mut egui::Ui) {
        let (response, painter) =
            ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());

        let rect = response.rect;
        let size = rect.size();

        // Dark background
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(25, 25, 35));

        // --- Mouse interaction ---
        if response.dragged_by(egui::PointerButton::Primary) {
            let d = response.drag_delta();
            self.yaw += d.x * 0.005;
            self.pitch = (self.pitch + d.y * 0.005).clamp(
                -std::f32::consts::FRAC_PI_2 + 0.01,
                std::f32::consts::FRAC_PI_2 - 0.01,
            );
        }
        if response.dragged_by(egui::PointerButton::Secondary)
            || response.dragged_by(egui::PointerButton::Middle)
        {
            self.pan += response.drag_delta();
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll.abs() > 0.1 {
                self.zoom *= 1.0 + scroll * 0.003;
                self.zoom = self.zoom.clamp(0.0001, 1000.0);
            }
        }

        // --- Empty state ---
        if self.preview_lines.is_empty() {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Load a G-code file to preview",
                egui::FontId::proportional(18.0),
                egui::Color32::from_rgb(100, 100, 120),
            );
            return;
        }

        let offset = rect.min.to_vec2();

        // --- Build-plate grid ---
        self.draw_grid(&painter, size, offset);

        // --- Extrusion paths ---
        let expanded = rect.expand(100.0);
        for line in &self.preview_lines {
            if line.z > self.layer_cutoff {
                continue;
            }
            let p1 = self.project(line.start, size) + offset;
            let p2 = self.project(line.end, size) + offset;
            if !expanded.contains(p1) && !expanded.contains(p2) {
                continue;
            }
            let color = z_to_color(line.z, self.min_z, self.max_z);
            painter.line_segment([p1, p2], egui::Stroke::new(1.5, color));
        }

        // --- Axes gizmo ---
        self.draw_axes(&painter, rect);

        // --- Help hint ---
        painter.text(
            rect.min + egui::vec2(8.0, 8.0),
            egui::Align2::LEFT_TOP,
            "LMB: Orbit  |  RMB/MMB: Pan  |  Scroll: Zoom",
            egui::FontId::proportional(11.0),
            egui::Color32::from_rgb(80, 80, 100),
        );
    }

    fn draw_grid(&self, painter: &egui::Painter, size: egui::Vec2, offset: egui::Vec2) {
        let color = egui::Color32::from_rgba_premultiplied(50, 50, 70, 80);

        let x0 = (self.bounds_min[0] / 10.0).floor() as i32 * 10;
        let x1 = (self.bounds_max[0] / 10.0).ceil() as i32 * 10;
        let y0 = (self.bounds_min[1] / 10.0).floor() as i32 * 10;
        let y1 = (self.bounds_max[1] / 10.0).ceil() as i32 * 10;

        let mut y = y0;
        while y <= y1 {
            let a = self.project([x0 as f32, y as f32, 0.0], size) + offset;
            let b = self.project([x1 as f32, y as f32, 0.0], size) + offset;
            painter.line_segment([a, b], egui::Stroke::new(0.5, color));
            y += 10;
        }
        let mut x = x0;
        while x <= x1 {
            let a = self.project([x as f32, y0 as f32, 0.0], size) + offset;
            let b = self.project([x as f32, y1 as f32, 0.0], size) + offset;
            painter.line_segment([a, b], egui::Stroke::new(0.5, color));
            x += 10;
        }
    }

    fn draw_axes(&self, painter: &egui::Painter, rect: egui::Rect) {
        let origin = rect.min + egui::vec2(50.0, rect.height() - 50.0);
        let len = 30.0;

        let project_dir = |d: [f32; 3]| -> egui::Vec2 {
            let (sy, cy) = self.yaw.sin_cos();
            let x1 = d[0] * cy - d[1] * sy;
            let y1 = d[0] * sy + d[1] * cy;
            let z1 = d[2];
            let (sp, cp) = self.pitch.sin_cos();
            egui::vec2(x1 * len, -(y1 * sp + z1 * cp) * len)
        };

        let red = egui::Color32::from_rgb(220, 60, 60);
        let green = egui::Color32::from_rgb(60, 220, 60);
        let blue = egui::Color32::from_rgb(80, 80, 255);

        let xe = origin + project_dir([1.0, 0.0, 0.0]);
        let ye = origin + project_dir([0.0, 1.0, 0.0]);
        let ze = origin + project_dir([0.0, 0.0, 1.0]);

        painter.line_segment([origin, xe], egui::Stroke::new(2.0, red));
        painter.line_segment([origin, ye], egui::Stroke::new(2.0, green));
        painter.line_segment([origin, ze], egui::Stroke::new(2.0, blue));

        let f = egui::FontId::proportional(11.0);
        painter.text(xe, egui::Align2::LEFT_CENTER, "X", f.clone(), red);
        painter.text(ye, egui::Align2::LEFT_CENTER, "Y", f.clone(), green);
        painter.text(ze, egui::Align2::LEFT_CENTER, "Z", f, blue);
    }
}

// ---------------------------------------------------------------------------
// Projection helpers
// ---------------------------------------------------------------------------

impl App {
    /// Orthographic projection: world → screen position (relative to viewport origin).
    fn project(&self, p: [f32; 3], size: egui::Vec2) -> egui::Pos2 {
        let x = p[0] - self.center[0];
        let y = p[1] - self.center[1];
        let z = p[2] - self.center[2];

        // Yaw  – rotate around world-Z (the "up" axis in G-code space)
        let (sy, cy) = self.yaw.sin_cos();
        let x1 = x * cy - y * sy;
        let y1 = x * sy + y * cy;
        let z1 = z;

        // Pitch – tilt around the screen-X axis
        let (sp, cp) = self.pitch.sin_cos();
        let screen_x = x1;
        let screen_y = y1 * sp + z1 * cp; // vertical on screen

        let scale = self.zoom * size.x.min(size.y) * 0.5;

        egui::Pos2::new(
            size.x / 2.0 + screen_x * scale + self.pan.x,
            size.y / 2.0 - screen_y * scale + self.pan.y,
        )
    }

    fn auto_zoom(&mut self) {
        if self.preview_lines.is_empty() {
            return;
        }
        let extent = (self.bounds_max[0] - self.bounds_min[0])
            .max(self.bounds_max[1] - self.bounds_min[1])
            .max(self.bounds_max[2] - self.bounds_min[2])
            .max(1.0);
        self.zoom = 1.0 / extent;
        self.pan = egui::Vec2::ZERO;
    }
}

// ---------------------------------------------------------------------------
// Background conversion
// ---------------------------------------------------------------------------

impl App {
    /// Load, parse, and process – show paths in the preview (no mesh / STL).
    fn start_preview(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        self.busy = true;
        self.log.push("Loading preview…".into());

        let input = self.input_path.clone();
        let nozzle = self.nozzle_diameter;
        let layer = self.layer_height;
        let filament = self.filament_diameter;
        let sides = self.mesh_sides as usize;

        std::thread::spawn(move || {
            let start = Instant::now();

            let content = match std::fs::read_to_string(&input) {
                Ok(c) => c,
                Err(e) => {
                    let _ = tx.send(Msg::Error(format!("Read error: {}", e)));
                    return;
                }
            };
            let _ = tx.send(Msg::Log(format!("Read {} bytes", content.len())));

            let lines = match parser::parse_gcode(&content) {
                Ok(l) => l,
                Err(e) => {
                    let _ = tx.send(Msg::Error(format!("Parse error: {}", e)));
                    return;
                }
            };
            let _ = tx.send(Msg::Log(format!("Parsed {} G-code lines", lines.len())));

            let cfg = config::ConfigBuilder::new()
                .nozzle_diameter(nozzle)
                .layer_height(layer)
                .filament_diameter(filament)
                .mesh_sides(sides)
                .build();

            let mut proc = processor::GCodeProcessor::with_config(&cfg);
            let paths = proc.process(&lines);
            let _ = tx.send(Msg::Log(format!("Found {} extrusion paths", paths.len())));

            send_preview(&tx, &paths);

            let elapsed = start.elapsed();
            let _ = tx.send(Msg::Done(format!(
                "Preview loaded ({:.2}s)",
                elapsed.as_secs_f64()
            )));
        });
    }

    /// Full pipeline: parse → process → mesh → write STL.
    fn start_export(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        self.busy = true;
        self.log.push("Starting export…".into());

        let input = self.input_path.clone();
        let output = if self.output_path.is_empty() {
            let mut p = std::path::PathBuf::from(&input);
            p.set_extension("stl");
            p.display().to_string()
        } else {
            self.output_path.clone()
        };
        let nozzle = self.nozzle_diameter;
        let layer = self.layer_height;
        let filament = self.filament_diameter;
        let sides = self.mesh_sides as usize;

        std::thread::spawn(move || {
            let start = Instant::now();

            // Read
            let content = match std::fs::read_to_string(&input) {
                Ok(c) => c,
                Err(e) => {
                    let _ = tx.send(Msg::Error(format!("Read error: {}", e)));
                    return;
                }
            };
            let _ = tx.send(Msg::Log(format!("Read {} bytes", content.len())));

            // Parse
            let lines = match parser::parse_gcode(&content) {
                Ok(l) => l,
                Err(e) => {
                    let _ = tx.send(Msg::Error(format!("Parse error: {}", e)));
                    return;
                }
            };
            let _ = tx.send(Msg::Log(format!("Parsed {} G-code lines", lines.len())));

            // Process
            let cfg = config::ConfigBuilder::new()
                .nozzle_diameter(nozzle)
                .layer_height(layer)
                .filament_diameter(filament)
                .mesh_sides(sides)
                .build();

            let mut proc = processor::GCodeProcessor::with_config(&cfg);
            let paths = proc.process(&lines);
            let _ = tx.send(Msg::Log(format!("Found {} extrusion paths", paths.len())));

            // Send preview so the user sees paths while mesh generates
            send_preview(&tx, &paths);

            // Generate mesh
            let mesh = geometry::generate_mesh(&paths, &cfg);
            let _ = tx.send(Msg::Log(format!("Generated {} triangles", mesh.len())));

            // Write STL
            match std::fs::File::create(&output) {
                Ok(file) => {
                    let mut writer = std::io::BufWriter::new(file);
                    match export::write_stl(&mut writer, &mesh) {
                        Ok(()) => {
                            let elapsed = start.elapsed();
                            let _ = tx.send(Msg::Done(format!(
                                "Exported to {} ({:.2}s)",
                                output,
                                elapsed.as_secs_f64()
                            )));
                        }
                        Err(e) => {
                            let _ = tx.send(Msg::Error(format!("Write error: {}", e)));
                        }
                    }
                }
                Err(e) => {
                    let _ = tx.send(Msg::Error(format!("Create file error: {}", e)));
                }
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Free helpers
// ---------------------------------------------------------------------------

/// Extract lightweight line segments + metadata from extrusion paths.
fn send_preview(tx: &mpsc::Sender<Msg>, paths: &[ExtrusionPath]) {
    let mut lines = Vec::new();
    let mut min_z = f32::MAX;
    let mut max_z = f32::MIN;
    let mut bounds_min = [f32::MAX; 3];
    let mut bounds_max = [f32::MIN; 3];

    for path in paths {
        for w in path.nodes.windows(2) {
            let s = [w[0].pos.x, w[0].pos.y, w[0].pos.z];
            let e = [w[1].pos.x, w[1].pos.y, w[1].pos.z];
            let z = (s[2] + e[2]) * 0.5;

            min_z = min_z.min(z);
            max_z = max_z.max(z);

            for i in 0..3 {
                bounds_min[i] = bounds_min[i].min(s[i]).min(e[i]);
                bounds_max[i] = bounds_max[i].max(s[i]).max(e[i]);
            }

            lines.push(PreviewLine { start: s, end: e, z });
        }
    }

    if lines.is_empty() {
        min_z = 0.0;
        max_z = 1.0;
        bounds_min = [0.0; 3];
        bounds_max = [0.0; 3];
    }

    let center = [
        (bounds_min[0] + bounds_max[0]) * 0.5,
        (bounds_min[1] + bounds_max[1]) * 0.5,
        (bounds_min[2] + bounds_max[2]) * 0.5,
    ];

    let _ = tx.send(Msg::Preview {
        lines,
        min_z,
        max_z,
        center,
        bounds_min,
        bounds_max,
    });
}

/// Map a Z height to a rainbow colour (blue → cyan → green → yellow → red).
fn z_to_color(z: f32, min_z: f32, max_z: f32) -> egui::Color32 {
    let range = (max_z - min_z).max(0.001);
    let t = ((z - min_z) / range).clamp(0.0, 1.0);

    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (0.0, s, 1.0)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (0.0, 1.0, 1.0 - s)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (s, 1.0, 0.0)
    } else {
        let s = (t - 0.75) / 0.25;
        (1.0, 1.0 - s, 0.0)
    };

    egui::Color32::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}
