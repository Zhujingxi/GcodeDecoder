use std::io;
use std::path::PathBuf;
use std::sync::mpsc;

use crossterm::{
    event::{DisableMouseCapture, EnableMouseCapture},
    terminal::{self, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::{Backend, CrosstermBackend},
    Terminal,
};

use crate::config::Config;
use crate::tui::{
    config_panel::ConfigPanel,
    event::{handle_event, Action},
    preview_panel::PreviewPanel,
    ui, AppState, Field,
};
use crate::{export, geometry, parser, processor};

pub struct App {
    pub config_panel: ConfigPanel,
    pub preview_panel: PreviewPanel,
    pub state: AppState,
    pub focused_field: Field,
    pub input_path: PathBuf,
    pub should_quit: bool,
    pub conversion_receiver: Option<mpsc::Receiver<anyhow::Result<PathBuf>>>,
}

impl App {
    pub fn new(input_path: PathBuf, config: Config) -> anyhow::Result<Self> {
        let mut preview_panel = PreviewPanel::new();
        preview_panel.load_file(&input_path)?;

        Ok(Self {
            config_panel: ConfigPanel::new(config),
            preview_panel,
            state: AppState::Editing,
            focused_field: Field::NozzleDiameter,
            input_path,
            should_quit: false,
            conversion_receiver: None,
        })
    }

    pub fn load_file_stats(&mut self, segment_count: usize, layer_count: usize) {
        self.preview_panel.update_stats(segment_count, layer_count);
    }

    pub fn run(
        &mut self,
        terminal: &mut Terminal<impl Backend>,
    ) -> anyhow::Result<Option<PathBuf>> {
        let mut last_tick = std::time::Instant::now();
        let tick_rate = std::time::Duration::from_millis(250);

        while !self.should_quit {
            terminal.draw(|f| {
                ui::render(
                    f,
                    &self.config_panel,
                    &self.preview_panel,
                    &self.state,
                    self.focused_field,
                );
            })?;

            let timeout = tick_rate
                .checked_sub(last_tick.elapsed())
                .unwrap_or_else(|| std::time::Duration::from_secs(0));

            if crossterm::event::poll(timeout)? {
                if let Ok(event) = crossterm::event::read() {
                    let action = handle_event(event, &self.state, self.focused_field);
                    self.handle_action(action)?;
                }
            }

            if last_tick.elapsed() >= tick_rate {
                self.on_tick();
                last_tick = std::time::Instant::now();
            }
        }

        match &self.state {
            AppState::Complete { output_path } => Ok(Some(output_path.clone())),
            _ => Ok(None),
        }
    }

    fn handle_action(&mut self, action: Action) -> anyhow::Result<()> {
        match action {
            Action::Quit => self.should_quit = true,
            Action::NextField => {
                self.focused_field = self.focused_field.next();
            }
            Action::PrevField => {
                self.focused_field = self.focused_field.prev();
            }
            Action::Increment => self.adjust_value(true, false),
            Action::Decrement => self.adjust_value(false, false),
            Action::BigIncrement => self.adjust_value(true, true),
            Action::BigDecrement => self.adjust_value(false, true),
            Action::Toggle => {
                if self.focused_field == Field::RemoveInternal {
                    self.config_panel.toggle_remove_internal();
                }
            }
            Action::Convert => {
                if matches!(self.state, AppState::Editing) {
                    self.start_conversion()?;
                }
            }
            Action::Reset => {
                self.config_panel.reset_to_defaults();
            }
            Action::None => {}
        }
        Ok(())
    }

    fn adjust_value(&mut self, increment: bool, big: bool) {
        match self.focused_field {
            Field::NozzleDiameter => {
                if increment {
                    self.config_panel.increment_nozzle(big);
                } else {
                    self.config_panel.decrement_nozzle(big);
                }
            }
            Field::LayerHeight => {
                if increment {
                    self.config_panel.increment_layer(big);
                } else {
                    self.config_panel.decrement_layer(big);
                }
            }
            Field::FilamentDiameter => {
                if increment {
                    self.config_panel.increment_filament(big);
                } else {
                    self.config_panel.decrement_filament(big);
                }
            }
            Field::MeshSides => {
                if increment {
                    self.config_panel.increment_sides(big);
                } else {
                    self.config_panel.decrement_sides(big);
                }
            }
            Field::RemoveInternal => {
                self.config_panel.toggle_remove_internal();
            }
        }
    }

    fn start_conversion(&mut self) -> anyhow::Result<()> {
        self.state = AppState::Processing { progress: 0.0 };

        let content = std::fs::read_to_string(&self.input_path)?;
        let segments = parser::parse_gcode(&content)?;

        self.preview_panel.update_stats(segments.len(), 0);

        let config = self.config_panel.config.clone();
        let input_path = self.input_path.clone();

        let (tx, rx) = mpsc::channel();
        self.conversion_receiver = Some(rx);

        std::thread::spawn(move || {
            let result = Self::process_conversion(content, config, input_path);
            let _ = tx.send(result);
        });

        self.state = AppState::Processing { progress: 0.5 };

        Ok(())
    }

    fn process_conversion(
        content: String,
        config: Config,
        input_path: PathBuf,
    ) -> anyhow::Result<PathBuf> {
        let segments = parser::parse_gcode(&content)?;

        let mut processor = processor::GCodeProcessor::with_config(&config);
        let mut processed_segments = processor.process(&segments);

        if config.remove_internal {
            processed_segments = crate::optimizer::optimize_paths(processed_segments, &config);
        }

        let mesh = geometry::generate_mesh(&processed_segments, &config);

        let output_path = {
            let mut p = input_path.clone();
            p.set_extension("stl");
            p
        };

        let file = std::fs::File::create(&output_path)?;
        let mut writer = std::io::BufWriter::new(file);
        export::write_stl(&mut writer, &mesh)?;

        Ok(output_path)
    }

    fn on_tick(&mut self) {
        if let AppState::Processing { progress: _ } = self.state {
            if let Some(ref receiver) = self.conversion_receiver {
                match receiver.try_recv() {
                    Ok(Ok(output_path)) => {
                        self.state = AppState::Complete { output_path };
                        self.conversion_receiver = None;
                    }
                    Ok(Err(e)) => {
                        self.state = AppState::Error {
                            message: format!("Conversion failed: {}", e),
                        };
                        self.conversion_receiver = None;
                    }
                    Err(mpsc::TryRecvError::Empty) => {
                        // Still processing, keep waiting
                    }
                    Err(mpsc::TryRecvError::Disconnected) => {
                        self.state = AppState::Error {
                            message: "Conversion thread panicked".to_string(),
                        };
                        self.conversion_receiver = None;
                    }
                }
            }
        }
    }
}

pub fn setup_terminal() -> anyhow::Result<Terminal<impl Backend>> {
    terminal::enable_raw_mode()?;
    let mut stdout = io::stdout();
    crossterm::execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let terminal = Terminal::new(backend)?;
    Ok(terminal)
}

pub fn restore_terminal<B: Backend>(terminal: &mut Terminal<B>) -> anyhow::Result<()> {
    terminal::disable_raw_mode()?;
    crossterm::execute!(io::stdout(), LeaveAlternateScreen, DisableMouseCapture)?;
    terminal.show_cursor()?;
    Ok(())
}

pub fn run_tui(input_path: PathBuf, config: Config) -> anyhow::Result<Option<PathBuf>> {
    let mut terminal = setup_terminal()?;
    let mut app = App::new(input_path, config)?;

    let content = std::fs::read_to_string(&app.input_path)?;
    if let Ok(lines) = parser::parse_gcode(&content) {
        // Count unique Z values as proxy for layer count
        let mut z_values: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for line in &lines {
            if line.command == "G1" || line.command == "G0" {
                if let Some((_, z_val)) = line.params.iter().find(|(k, _)| *k == 'Z') {
                    // Round to 3 decimal places and convert to integer for hashing
                    let z_key = (z_val * 1000.0).round() as u32;
                    z_values.insert(z_key);
                }
            }
        }
        app.load_file_stats(lines.len(), z_values.len());
    }

    let result = app.run(&mut terminal);
    restore_terminal(&mut terminal)?;
    result
}
