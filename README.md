# smiles-to-prolog-parser

Converts [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) molecular notation into Prolog term representations, designed for use in logic programming and inductive learning systems.

## Overview

This crate parses SMILES strings through three stages:

1. **Tokenization** — Breaks SMILES into atoms, bonds, branches, and ring closures
2. **Graph construction** — Builds a molecular graph with implicit hydrogen computation, aromatic bond inference, and ring detection
3. **Term generation** — Converts the graph into nested Prolog terms with canonical ordering

## Example

```rust
use smiles_to_prolog_parser::smiles_to_prolog;

// Ethanol: CCO
let term = smiles_to_prolog("CCO").unwrap();
assert_eq!(term, "('n/a', c, [(1, c, [h, h, (1, o, [h])]), h, h, h])");

// Benzene: c1ccccc1
let term = smiles_to_prolog("c1ccccc1").unwrap();
// mol(ring([...]))  — ring compounds use a ring root
```

## Output Format

Molecules are represented as `(bond, element, [substituents])` where substituents are ordered simplest-first: *Bond at root of molecule will be n/a*

- **Hydrogen** appears as bare `h` (shorthand for `(1,h,[])`)
- **Atom substituents** appear as `(bond, element, [subs])` where bond is `1` (single), `2` (double), `3` (triple), or `a` (aromatic)
- **Ring structures** appear as `ring(bond, [members])` where each member is either `(bond, element, [subs])` or `ring(bond, [members])`

## CLI Usage

```sh
# Convert SMILES strings to Prolog terms
cargo run -- CCO "c1ccccc1" "CC(=O)O"

# Read from stdin (one SMILES per line)
echo -e "CCO\nc1ccccc1" | cargo run
```

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
smiles-to-prolog-parser = "0.1.2"
```

## License

MIT
