pub mod token;
pub mod graph;
pub mod term;

use crate::term::generate_term;

/// Generate a Prolog term string from a SMILES string.
pub fn smiles_to_prolog(smiles: &str) -> Result<String, String> {
    let graph = crate::graph::build_from_smiles(smiles)?;
    let term = generate_term(&graph);
    Ok(term.to_string())
}