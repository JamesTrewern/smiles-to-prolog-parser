mod token;
mod graph;
use smiles_parser::token::tokenize;
use smiles_parser::graph::build_from_smiles;

fn main() {
    let examples = [
        ("methane", "C"),
        ("ethanol", "CCO"),
        ("acetic acid", "CC(=O)O"),
        ("benzene", "c1ccccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("toluene", "Cc1ccccc1"),
        ("biphenyl", "c1ccccc1-c2ccccc2"),
        ("pyridine", "c1ccncc1"),
    ];

    for (name, smiles) in &examples {
        println!("=== {} ({}) ===", name, smiles);
        match tokenize(smiles) {
            Ok(tokens) => println!("  tokens: {:?}", tokens),
            Err(e) => { println!("  tokenize error: {}", e); continue; }
        }
        match build_from_smiles(smiles) {
            Ok(g) => {
                println!("  atoms: {}, ring_bonds: {:?}", g.atom_count(), g.ring_bonds);
                for (i, a) in g.atoms.iter().enumerate() {
                    let nbrs: Vec<String> = g.neighbors(i).iter()
                        .map(|e| format!("{}({:?})", e.target, e.bond_type))
                        .collect();
                    println!("    [{}] {} aro={} h={} nbrs=[{}]",
                        i, a.element, a.aromatic, a.implicit_hcount, nbrs.join(", "));
                }
            }
            Err(e) => println!("  graph error: {}", e),
        }
        println!();
    }
}
