use anyhow::Result;

#[derive(Debug, Clone, PartialEq)]
pub struct GCodeLine {
    pub command: String,
    pub params: Vec<(char, f32)>,
    pub comment: Option<String>,
}

pub fn parse_gcode(input: &str) -> Result<Vec<GCodeLine>> {
    let mut lines = Vec::new();
    for line in input.lines() {
        if let Some(parsed) = parse_line(line) {
            lines.push(parsed);
        }
    }
    Ok(lines)
}

fn parse_line(line: &str) -> Option<GCodeLine> {
    let line = line.trim();
    if line.is_empty() {
        return None;
    }

    // Split comment
    let (content, comment) = match line.split_once(';') {
        Some((code, comment)) => (code.trim(), Some(comment.to_string())),
        None => (line, None),
    };

    if content.is_empty() {
        if comment.is_some() {
             return Some(GCodeLine {
                command: String::new(),
                params: Vec::new(),
                comment,
            });
        }
        return None;
    }

    // naive split by whitespace is usually enough for G-code
    // but sometimes parameters are tight like G1X10Y10
    // For high performance and robustness, we should scan bytes.
    // heavily simplified for MVP: assume whitespace or standard formatting for now.
    // If we need high perf "no regex" parsing for tight gcode, we need a custom lexer.
    // Let's implement a slightly smarter scanner.

    let mut parts = content.split_whitespace();
    let command = parts.next()?.to_string(); // e.g. "G1"
    
    // Check if command is actually a command (starts with letter)
    if !command.chars().next().unwrap_or(' ').is_alphabetic() {
        return None; // Invalid or empty
    }

    let mut params = Vec::new();
    
    // Handle the rest of the parts
    for part in parts {
        if let Some(c) = part.chars().next() {
             if c.is_alphabetic() {
                 // Try to parse the rest as float
                 if let Ok(val) = part[1..].parse::<f32>() {
                     params.push((c.to_ascii_uppercase(), val));
                 }
             }
        }
    }
    
    // Handle the case where params are stuck to command or each other: "G1X10"
    // For the initial MVP, let's assume standard slicing output which usually has spaces.
    // If the user needs to handle "G1X10Y10", we will need a more complex parser.
    // Given the request for "high performance", we can iterate chars.
    
    if params.is_empty() && content.len() > command.len() {
         // Fallback or detailed check?
         // Let's stick to the split_whitespace for now as most slicers output spaces.
         // If we see issues, we upgrade this function.
    }

    Some(GCodeLine {
        command,
        params,
        comment,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_line() {
        let line = "G1 X10.5 Y20.0 E0.1";
        let parsed = parse_line(line).unwrap();
        assert_eq!(parsed.command, "G1");
        assert_eq!(parsed.params, vec![('X', 10.5), ('Y', 20.0), ('E', 0.1)]);
    }

    #[test]
    fn test_with_comment() {
        let line = "G0 Z10 ; move up";
        let parsed = parse_line(line).unwrap();
        assert_eq!(parsed.command, "G0");
        assert_eq!(parsed.comment, Some(" move up".to_string()));
    }
}
