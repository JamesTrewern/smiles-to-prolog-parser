/// Bond type as seen in SMILES syntax
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondKind {
    Single,   // - (usually implicit)
    Double,   // =
    Triple,   // #
    Aromatic, // :
    Up,       // /  (geometric stereo)
    Down,     // \  (geometric stereo)
}

/// An atom from the "organic subset" — written without brackets
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrganicAtom {
    B, C, N, O, P, S, F, Cl, Br, I,
}

/// An aromatic organic atom — lowercase in SMILES
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticAtom {
    B, C, N, O, P, S,
}

/// Contents parsed from a bracket atom like [NH4+], [13C@@H], [Fe+2]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BracketAtom {
    pub isotope: Option<u16>,
    pub symbol: String,
    pub aromatic: bool,
    pub chirality: Option<Chirality>,
    pub hcount: Option<u8>,
    pub charge: i8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Chirality {
    Anticlockwise, // @
    Clockwise,     // @@
}

/// A single SMILES token
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    Organic(OrganicAtom),
    Aromatic(AromaticAtom),
    Bracket(BracketAtom),
    Bond(BondKind),
    RingBond(u8),      // ring closure digit: 0-9
    RingBond2(u8),     // %nn ring closure: 10-99
    BranchOpen,        // (
    BranchClose,       // )
    Dot,               // . (fragment separator)
}

/// Tokenizer error
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TokenError {
    pub pos: usize,
    pub msg: String,
}

impl std::fmt::Display for TokenError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "position {}: {}", self.pos, self.msg)
    }
}

/// Tokenize a SMILES string into a Vec<Token>.
pub fn tokenize(smiles: &str) -> Result<Vec<Token>, TokenError> {
    let chars: Vec<char> = smiles.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < chars.len() {
        let ch = chars[i];
        let start = i;

        match ch {
            // Branch delimiters
            '(' => { tokens.push(Token::BranchOpen);  i += 1; }
            ')' => { tokens.push(Token::BranchClose); i += 1; }

            // Fragment separator
            '.' => { tokens.push(Token::Dot); i += 1; }

            // Explicit bonds
            '-' => { tokens.push(Token::Bond(BondKind::Single));   i += 1; }
            '=' => { tokens.push(Token::Bond(BondKind::Double));   i += 1; }
            '#' => { tokens.push(Token::Bond(BondKind::Triple));   i += 1; }
            ':' => { tokens.push(Token::Bond(BondKind::Aromatic)); i += 1; }
            '/' => { tokens.push(Token::Bond(BondKind::Up));       i += 1; }
            '\\' => { tokens.push(Token::Bond(BondKind::Down));    i += 1; }

            // Ring closure digits
            '0'..='9' => {
                tokens.push(Token::RingBond(ch as u8 - b'0'));
                i += 1;
            }

            // Two-digit ring closure: %nn
            '%' => {
                i += 1;
                let d1 = chars.get(i).copied().ok_or_else(|| TokenError {
                    pos: start, msg: "expected digits after %".into(),
                })?;
                let d2 = chars.get(i + 1).copied().ok_or_else(|| TokenError {
                    pos: start, msg: "expected two digits after %".into(),
                })?;
                if !d1.is_ascii_digit() || !d2.is_ascii_digit() {
                    return Err(TokenError {
                        pos: start, msg: "expected two digits after %".into(),
                    });
                }
                let n = (d1 as u8 - b'0') * 10 + (d2 as u8 - b'0');
                tokens.push(Token::RingBond2(n));
                i += 2;
            }

            // Bracket atom: [...]
            '[' => {
                let bracket = parse_bracket_atom(&chars, &mut i)?;
                tokens.push(Token::Bracket(bracket));
            }

            // Aromatic atoms (lowercase)
            'b' => { tokens.push(Token::Aromatic(AromaticAtom::B)); i += 1; }
            'c' => { tokens.push(Token::Aromatic(AromaticAtom::C)); i += 1; }
            'n' => { tokens.push(Token::Aromatic(AromaticAtom::N)); i += 1; }
            'o' => { tokens.push(Token::Aromatic(AromaticAtom::O)); i += 1; }
            'p' => { tokens.push(Token::Aromatic(AromaticAtom::P)); i += 1; }
            's' => { tokens.push(Token::Aromatic(AromaticAtom::S)); i += 1; }

            // Organic subset atoms (uppercase, with two-char lookahead for Cl/Br)
            'B' => {
                if chars.get(i + 1) == Some(&'r') {
                    tokens.push(Token::Organic(OrganicAtom::Br));
                    i += 2;
                } else {
                    tokens.push(Token::Organic(OrganicAtom::B));
                    i += 1;
                }
            }
            'C' => {
                if chars.get(i + 1) == Some(&'l') {
                    tokens.push(Token::Organic(OrganicAtom::Cl));
                    i += 2;
                } else {
                    tokens.push(Token::Organic(OrganicAtom::C));
                    i += 1;
                }
            }
            'N' => { tokens.push(Token::Organic(OrganicAtom::N)); i += 1; }
            'O' => { tokens.push(Token::Organic(OrganicAtom::O)); i += 1; }
            'P' => { tokens.push(Token::Organic(OrganicAtom::P)); i += 1; }
            'S' => { tokens.push(Token::Organic(OrganicAtom::S)); i += 1; }
            'F' => { tokens.push(Token::Organic(OrganicAtom::F)); i += 1; }
            'I' => { tokens.push(Token::Organic(OrganicAtom::I)); i += 1; }

            _ => {
                return Err(TokenError {
                    pos: start,
                    msg: format!("unexpected character: '{}'", ch),
                });
            }
        }
    }

    Ok(tokens)
}

/// Parse a bracket atom starting at chars[i] == '['.
/// Advances i past the closing ']'.
fn parse_bracket_atom(chars: &[char], i: &mut usize) -> Result<BracketAtom, TokenError> {
    let start = *i;
    *i += 1; // skip '['

    // Isotope (optional leading digits)
    let isotope = parse_number(chars, i);

    // Element symbol
    let (symbol, aromatic) = parse_element_symbol(chars, i, start)?;

    // Chirality (optional @ or @@)
    let chirality = if chars.get(*i) == Some(&'@') {
        *i += 1;
        if chars.get(*i) == Some(&'@') {
            *i += 1;
            Some(Chirality::Clockwise)
        } else {
            Some(Chirality::Anticlockwise)
        }
    } else {
        None
    };

    // Hydrogen count (optional H or Hn)
    let hcount = if chars.get(*i) == Some(&'H') {
        *i += 1;
        let n = parse_number(chars, i);
        Some(n.unwrap_or(1) as u8)
    } else {
        None
    };

    // Charge (optional + or - with optional digit)
    let charge = parse_charge(chars, i);

    // Expect closing ']'
    if chars.get(*i) != Some(&']') {
        return Err(TokenError {
            pos: start,
            msg: format!("unclosed bracket atom (at pos {})", *i),
        });
    }
    *i += 1;

    Ok(BracketAtom {
        isotope: isotope.map(|n| n as u16),
        symbol,
        aromatic,
        chirality,
        hcount,
        charge,
    })
}

/// Parse an optional run of digits, returning the number if present.
fn parse_number(chars: &[char], i: &mut usize) -> Option<u32> {
    let start = *i;
    while *i < chars.len() && chars[*i].is_ascii_digit() {
        *i += 1;
    }
    if *i > start {
        let s: String = chars[start..*i].iter().collect();
        Some(s.parse().unwrap())
    } else {
        None
    }
}

/// Parse element symbol. Returns (symbol_string, is_aromatic).
fn parse_element_symbol(
    chars: &[char], i: &mut usize, start: usize,
) -> Result<(String, bool), TokenError> {
    let ch = *chars.get(*i).ok_or_else(|| TokenError {
        pos: start, msg: "expected element symbol in bracket atom".into(),
    })?;

    if ch.is_ascii_uppercase() {
        // Normal atom: one uppercase + optional one lowercase
        let mut sym = String::new();
        sym.push(ch);
        *i += 1;
        // Take one lowercase letter as part of symbol (Fe, Cl, etc.)
        // but not 'H' which could be hydrogen count
        if *i < chars.len() && chars[*i].is_ascii_lowercase() {
            sym.push(chars[*i]);
            *i += 1;
        }
        Ok((sym, false))
    } else if ch.is_ascii_lowercase() {
        // Aromatic bracket atom: e.g. [se], [as]
        let mut sym = String::new();
        sym.push(ch);
        *i += 1;
        if *i < chars.len() && chars[*i].is_ascii_lowercase() {
            sym.push(chars[*i]);
            *i += 1;
        }
        Ok((sym, true))
    } else if ch == '*' {
        // Wildcard atom
        *i += 1;
        Ok(("*".into(), false))
    } else {
        Err(TokenError {
            pos: *i,
            msg: format!("expected element symbol, got '{}'", ch),
        })
    }
}

/// Parse optional charge: +, -, +2, -3, ++, --, etc.
fn parse_charge(chars: &[char], i: &mut usize) -> i8 {
    match chars.get(*i) {
        Some(&'+') => {
            *i += 1;
            if let Some(n) = parse_number(chars, i) {
                n as i8
            } else {
                let mut count: i8 = 1;
                while chars.get(*i) == Some(&'+') {
                    count += 1;
                    *i += 1;
                }
                count
            }
        }
        Some(&'-') => {
            *i += 1;
            if let Some(n) = parse_number(chars, i) {
                -(n as i8)
            } else {
                let mut count: i8 = -1;
                while chars.get(*i) == Some(&'-') {
                    count -= 1;
                    *i += 1;
                }
                count
            }
        }
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methane() {
        let tokens = tokenize("C").unwrap();
        assert_eq!(tokens, vec![Token::Organic(OrganicAtom::C)]);
    }

    #[test]
    fn test_ethanol() {
        let tokens = tokenize("CCO").unwrap();
        assert_eq!(tokens, vec![
            Token::Organic(OrganicAtom::C),
            Token::Organic(OrganicAtom::C),
            Token::Organic(OrganicAtom::O),
        ]);
    }

    #[test]
    fn test_acetic_acid() {
        let tokens = tokenize("CC(=O)O").unwrap();
        assert_eq!(tokens, vec![
            Token::Organic(OrganicAtom::C),
            Token::Organic(OrganicAtom::C),
            Token::BranchOpen,
            Token::Bond(BondKind::Double),
            Token::Organic(OrganicAtom::O),
            Token::BranchClose,
            Token::Organic(OrganicAtom::O),
        ]);
    }

    #[test]
    fn test_benzene() {
        let tokens = tokenize("c1ccccc1").unwrap();
        assert_eq!(tokens, vec![
            Token::Aromatic(AromaticAtom::C),
            Token::RingBond(1),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::RingBond(1),
        ]);
    }

    #[test]
    fn test_chloroform() {
        let tokens = tokenize("ClC(Cl)Cl").unwrap();
        assert_eq!(tokens, vec![
            Token::Organic(OrganicAtom::Cl),
            Token::Organic(OrganicAtom::C),
            Token::BranchOpen,
            Token::Organic(OrganicAtom::Cl),
            Token::BranchClose,
            Token::Organic(OrganicAtom::Cl),
        ]);
    }

    #[test]
    fn test_bracket_atom_simple() {
        let tokens = tokenize("[NH4+]").unwrap();
        assert_eq!(tokens, vec![
            Token::Bracket(BracketAtom {
                isotope: None,
                symbol: "N".into(),
                aromatic: false,
                chirality: None,
                hcount: Some(4),
                charge: 1,
            }),
        ]);
    }

    #[test]
    fn test_bracket_isotope_chirality() {
        let tokens = tokenize("[13C@@H]").unwrap();
        assert_eq!(tokens, vec![
            Token::Bracket(BracketAtom {
                isotope: Some(13),
                symbol: "C".into(),
                aromatic: false,
                chirality: Some(Chirality::Clockwise),
                hcount: Some(1),
                charge: 0,
            }),
        ]);
    }

    #[test]
    fn test_bracket_iron() {
        let tokens = tokenize("[Fe+2]").unwrap();
        assert_eq!(tokens, vec![
            Token::Bracket(BracketAtom {
                isotope: None,
                symbol: "Fe".into(),
                aromatic: false,
                chirality: None,
                hcount: None,
                charge: 2,
            }),
        ]);
    }

    #[test]
    fn test_ring_two_digit() {
        let tokens = tokenize("C%10CC%10").unwrap();
        assert_eq!(tokens, vec![
            Token::Organic(OrganicAtom::C),
            Token::RingBond2(10),
            Token::Organic(OrganicAtom::C),
            Token::Organic(OrganicAtom::C),
            Token::RingBond2(10),
        ]);
    }

    #[test]
    fn test_pyridine() {
        let tokens = tokenize("c1ccncc1").unwrap();
        assert_eq!(tokens, vec![
            Token::Aromatic(AromaticAtom::C),
            Token::RingBond(1),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::N),
            Token::Aromatic(AromaticAtom::C),
            Token::Aromatic(AromaticAtom::C),
            Token::RingBond(1),
        ]);
    }

    #[test]
    fn test_naphthalene() {
        let tokens = tokenize("c1ccc2ccccc2c1").unwrap();
        assert_eq!(tokens.len(), 14);
        assert_eq!(tokens[5], Token::RingBond(2));
        assert_eq!(tokens[11], Token::RingBond(2));
        assert_eq!(tokens[0], Token::Aromatic(AromaticAtom::C));
    }

    #[test]
    fn test_error_unexpected_char() {
        let result = tokenize("C&C");
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().pos, 1);
    }
}
