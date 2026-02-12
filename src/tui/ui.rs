use ratatui::{
    layout::{Alignment, Constraint, Direction, Layout, Margin, Rect},
    style::{Color, Modifier, Style, Stylize},
    text::{Line, Span, Text},
    widgets::{Block, Borders, Clear, Paragraph, Wrap},
    Frame,
};

use crate::tui::{config_panel::ConfigPanel, preview_panel::PreviewPanel, AppState, Field};

pub fn render(
    frame: &mut Frame,
    config_panel: &ConfigPanel,
    preview_panel: &PreviewPanel,
    state: &AppState,
    focused_field: Field,
) {
    let main_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(10),
            Constraint::Length(1),
        ])
        .split(frame.area());

    render_header(frame, main_layout[0], preview_panel);
    render_body(
        frame,
        main_layout[1],
        config_panel,
        preview_panel,
        state,
        focused_field,
    );
    render_status_bar(frame, main_layout[2], state);

    if let AppState::Error { message } = state {
        render_error_modal(frame, message);
    }
}

fn render_header(frame: &mut Frame, area: Rect, preview_panel: &PreviewPanel) {
    let file_name = preview_panel.file_name();
    let header_text = format!("GcodeDecoder TUI - {}", file_name);

    let header = Paragraph::new(header_text)
        .style(Style::default().fg(Color::White).bg(Color::Blue))
        .alignment(Alignment::Center);

    frame.render_widget(header, area);
}

fn render_body(
    frame: &mut Frame,
    area: Rect,
    config_panel: &ConfigPanel,
    preview_panel: &PreviewPanel,
    state: &AppState,
    focused_field: Field,
) {
    let body_layout = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(35), Constraint::Percentage(65)])
        .split(area);

    render_config_panel(frame, body_layout[0], config_panel, state, focused_field);
    render_preview_panel(frame, body_layout[1], preview_panel, state);
}

fn render_config_panel(
    frame: &mut Frame,
    area: Rect,
    panel: &ConfigPanel,
    state: &AppState,
    focused_field: Field,
) {
    let config = &panel.config;
    let is_editing = matches!(state, AppState::Editing);

    let block = Block::default()
        .title(" Configuration ")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Blue));

    let inner = block.inner(area);
    frame.render_widget(block, area);

    let items = vec![
        (
            "Nozzle",
            format!("{:.2} mm", config.nozzle_diameter),
            Field::NozzleDiameter,
        ),
        (
            "Layer",
            format!("{:.2} mm", config.layer_height),
            Field::LayerHeight,
        ),
        (
            "Filament",
            format!("{:.2} mm", config.filament_diameter),
            Field::FilamentDiameter,
        ),
        ("Sides", format!("{}", config.mesh_sides), Field::MeshSides),
        (
            "Remove Internal",
            if config.remove_internal {
                "[✓]"
            } else {
                "[ ]"
            }
            .to_string(),
            Field::RemoveInternal,
        ),
    ];

    let mut text_lines: Vec<Line> = items
        .into_iter()
        .map(|(label, value, field)| {
            let is_focused = focused_field == field;
            let style = if is_focused && is_editing {
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };

            let arrow = if is_focused && is_editing { " ←" } else { "" };
            Line::from(vec![
                Span::styled(format!("{:15}", label), style),
                Span::styled(format!("{}{}", value, arrow), style),
            ])
        })
        .collect();

    text_lines.push(Line::from(""));
    text_lines.push(Line::from("[Defaults] [Revert]").fg(Color::DarkGray));

    if let Some(warning) = panel.validation_warning() {
        text_lines.push(Line::from(""));
        text_lines.push(Line::from(warning).fg(Color::Yellow));
    }

    let paragraph = Paragraph::new(Text::from(text_lines)).wrap(Wrap { trim: true });

    frame.render_widget(paragraph, inner);
}

fn render_preview_panel(frame: &mut Frame, area: Rect, panel: &PreviewPanel, state: &AppState) {
    let block = Block::default()
        .title(" Preview ")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Blue));

    let inner = block.inner(area);
    frame.render_widget(block, area);

    let mut lines = vec![
        Line::from(format!("Input: {}", panel.file_name())),
        Line::from(format!("Size: {}", panel.file_size())),
        Line::from(""),
        Line::from(format!("Segments: {}", panel.segment_count())),
        Line::from(format!("Layers: {}", panel.layer_count())),
    ];

    match state {
        AppState::Processing { progress } => {
            lines.push(Line::from(""));
            lines.push(Line::from(format!("Processing: {:.0}%", progress * 100.0)));
        }
        AppState::Complete { output_path } => {
            lines.push(Line::from(""));
            lines.push(Line::from("✓ Complete!").fg(Color::Green));
            lines.push(Line::from(format!("Output: {}", output_path.display())).fg(Color::Green));
        }
        _ => {}
    }

    lines.push(Line::from(""));
    lines.push(
        Line::from("[Convert]")
            .fg(Color::Green)
            .add_modifier(Modifier::BOLD),
    );

    let paragraph = Paragraph::new(Text::from(lines));
    frame.render_widget(paragraph, inner);
}

fn render_status_bar(frame: &mut Frame, area: Rect, state: &AppState) {
    let help_text = match state {
        AppState::Editing => "↑↓ Adjust | Tab Next | Enter Convert | r Reset | q Quit",
        AppState::Processing { .. } => "Processing... | q Cancel",
        AppState::Complete { .. } => "Complete! | Enter/q Quit",
        AppState::Error { .. } => "Error! | Enter/q Dismiss",
    };

    let status =
        Paragraph::new(help_text).style(Style::default().fg(Color::White).bg(Color::DarkGray));

    frame.render_widget(status, area);
}

fn render_error_modal(frame: &mut Frame, message: &str) {
    let area = frame.area();
    let popup_area = centered_rect(60, 40, area);

    frame.render_widget(Clear, popup_area);

    let block = Block::default()
        .title(" Error ")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Red));

    let inner = block.inner(popup_area);
    frame.render_widget(block, popup_area);

    let text = Text::from(vec![
        Line::from(message).fg(Color::Red),
        Line::from(""),
        Line::from("Press Enter or q to dismiss").fg(Color::Gray),
    ]);

    let paragraph = Paragraph::new(text).wrap(Wrap { trim: true });
    frame.render_widget(paragraph, inner);
}

fn centered_rect(percent_x: u16, percent_y: u16, r: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}
