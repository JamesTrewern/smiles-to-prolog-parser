use smiles_to_prolog_parser::smiles_to_prolog;
use std::io::{self, BufRead};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.is_empty() {
        // Read SMILES from stdin, one per line
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            match line {
                Ok(smiles) => {
                    let smiles = smiles.trim().to_string();
                    if smiles.is_empty() || smiles.starts_with('#') {
                        continue;
                    }
                    convert(&smiles);
                }
                Err(e) => {
                    eprintln!("Error reading stdin: {}", e);
                    std::process::exit(1);
                }
            }
        }
    } else {
        // Convert each argument as a SMILES string
        for smiles in &args {
            convert(smiles);
        }
    }
}

fn convert(smiles: &str) {
    match smiles_to_prolog(smiles) {
        Ok(term) => println!("{}", term),
        Err(e) => eprintln!("Error [{}]: {}", smiles, e),
    }
}
