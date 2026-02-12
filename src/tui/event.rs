use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyModifiers};
use std::time::Duration;

use crate::tui::{AppState, Field};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    Quit,
    NextField,
    PrevField,
    Increment,
    Decrement,
    BigIncrement,
    BigDecrement,
    Toggle,
    Convert,
    Reset,
    None,
}

pub fn handle_event(state: &AppState, focused_field: Field) -> anyhow::Result<Action> {
    if event::poll(Duration::from_millis(50))? {
        if let Event::Key(key) = event::read()? {
            return Ok(map_key_to_action(key, state, focused_field));
        }
    }
    Ok(Action::None)
}

fn map_key_to_action(key: KeyEvent, state: &AppState, focused_field: Field) -> Action {
    use KeyCode::*;

    match state {
        AppState::Complete { .. } | AppState::Error { .. } => match key.code {
            Char('q') | Esc => Action::Quit,
            Enter => Action::Quit,
            _ => Action::None,
        },
        AppState::Processing { .. } => match key.code {
            Char('q') | Esc => Action::Quit,
            _ => Action::None,
        },
        AppState::Editing => match key.code {
            Char('q') | Esc => Action::Quit,
            Char('r') => Action::Reset,
            Tab => {
                if key.modifiers.contains(KeyModifiers::SHIFT) {
                    Action::PrevField
                } else {
                    Action::NextField
                }
            }
            BackTab => Action::PrevField,
            Enter => Action::Convert,
            Up => {
                if key.modifiers.contains(KeyModifiers::CONTROL) {
                    Action::BigIncrement
                } else {
                    Action::Increment
                }
            }
            Down => {
                if key.modifiers.contains(KeyModifiers::CONTROL) {
                    Action::BigDecrement
                } else {
                    Action::Decrement
                }
            }
            Char(' ') => {
                if focused_field == Field::RemoveInternal {
                    Action::Toggle
                } else {
                    Action::None
                }
            }
            _ => Action::None,
        },
    }
}
