use std::collections::HashMap;
use crate::token::*;

/// Bond type in the molecular graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondType {
    pub fn order(&self) -> u8 {
        match self {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub element: String,
    pub aromatic: bool,
    pub charge: i8,
    pub isotope: Option<u16>,
    pub explicit_hcount: Option<u8>,
    pub implicit_hcount: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub target: usize,
    pub bond_type: BondType,
}

#[derive(Debug)]
pub struct MolGraph {
    pub atoms: Vec<Atom>,
    pub adj: Vec<Vec<Edge>>,
    pub ring_bonds: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GraphError {
    pub msg: String,
}

impl std::fmt::Display for GraphError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "graph error: {}", self.msg)
    }
}

impl MolGraph {
    fn new() -> Self { MolGraph { atoms: Vec::new(), adj: Vec::new(), ring_bonds: Vec::new() } }
    fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len(); self.atoms.push(atom); self.adj.push(Vec::new()); idx
    }
    fn add_bond(&mut self, a: usize, b: usize, bond_type: BondType) {
        self.adj[a].push(Edge { target: b, bond_type });
        self.adj[b].push(Edge { target: a, bond_type });
    }
    pub fn bond_between(&self, a: usize, b: usize) -> Option<BondType> {
        self.adj[a].iter().find(|e| e.target == b).map(|e| e.bond_type)
    }
    pub fn neighbors(&self, a: usize) -> &[Edge] { &self.adj[a] }
    pub fn atom_count(&self) -> usize { self.atoms.len() }
    pub fn bond_order_sum(&self, a: usize) -> u8 { self.adj[a].iter().map(|e| e.bond_type.order()).sum() }
    pub fn degree(&self, a: usize) -> usize { self.adj[a].len() }
}

fn atom_from_organic(oa: &OrganicAtom) -> Atom {
    let element = match oa {
        OrganicAtom::B => "b", OrganicAtom::C => "c", OrganicAtom::N => "n",
        OrganicAtom::O => "o", OrganicAtom::P => "p", OrganicAtom::S => "s",
        OrganicAtom::F => "f", OrganicAtom::Cl => "cl", OrganicAtom::Br => "br", OrganicAtom::I => "i",
    };
    Atom { element: element.to_string(), aromatic: false, charge: 0, isotope: None, explicit_hcount: None, implicit_hcount: 0 }
}

fn atom_from_aromatic(aa: &AromaticAtom) -> Atom {
    let element = match aa {
        AromaticAtom::B => "b", AromaticAtom::C => "c", AromaticAtom::N => "n",
        AromaticAtom::O => "o", AromaticAtom::P => "p", AromaticAtom::S => "s",
    };
    Atom { element: element.to_string(), aromatic: true, charge: 0, isotope: None, explicit_hcount: None, implicit_hcount: 0 }
}

fn atom_from_bracket(ba: &BracketAtom) -> Atom {
    // OpenSMILES: bracket atoms with no H field have 0 implicit hydrogens.
    // If hcount is None (no H specified), treat as explicitly 0.
    let hcount = Some(ba.hcount.unwrap_or(0));
    Atom { element: ba.symbol.to_lowercase(), aromatic: ba.aromatic, charge: ba.charge, isotope: ba.isotope, explicit_hcount: hcount, implicit_hcount: 0 }
}

fn bond_kind_to_type(bk: &BondKind) -> BondType {
    match bk {
        BondKind::Single => BondType::Single, BondKind::Double => BondType::Double,
        BondKind::Triple => BondType::Triple, BondKind::Aromatic => BondType::Aromatic,
        BondKind::Up | BondKind::Down => BondType::Single,
    }
}

fn resolve_bond(explicit: Option<BondKind>, atom_a: &Atom, atom_b: &Atom) -> BondType {
    if let Some(bk) = explicit { bond_kind_to_type(&bk) }
    else if atom_a.aromatic && atom_b.aromatic { BondType::Aromatic }
    else { BondType::Single }
}

fn resolve_ring_bond(open_bond: Option<BondKind>, close_bond: Option<BondKind>, atom_a: &Atom, atom_b: &Atom) -> Result<BondType, GraphError> {
    match (open_bond, close_bond) {
        (Some(a), Some(b)) => {
            let (ta, tb) = (bond_kind_to_type(&a), bond_kind_to_type(&b));
            if ta != tb { Err(GraphError { msg: "conflicting bond types on ring closure".into() }) }
            else { Ok(ta) }
        }
        (Some(bk), None) | (None, Some(bk)) => Ok(bond_kind_to_type(&bk)),
        (None, None) => Ok(resolve_bond(None, atom_a, atom_b)),
    }
}

fn default_valences(element: &str) -> &'static [u8] {
    match element {
        "b" => &[3], "c" => &[4], "n" => &[3, 5], "o" => &[2],
        "p" => &[3, 5], "s" => &[2, 4, 6],
        "f" => &[1], "cl" => &[1], "br" => &[1], "i" => &[1],
        _ => &[],
    }
}

fn compute_implicit_h(atom: &Atom, bond_order_sum: u8) -> u8 {
    if let Some(h) = atom.explicit_hcount { return h; }
    let valences = default_valences(&atom.element);
    if valences.is_empty() { return 0; }
    let valence = valences.iter().copied().find(|&v| v >= bond_order_sum).unwrap_or(*valences.last().unwrap());
    let base_h = valence.saturating_sub(bond_order_sum);
    if atom.aromatic { base_h.saturating_sub(1) } else { base_h }
}

pub fn build_graph(tokens: &[Token]) -> Result<MolGraph, GraphError> {
    let mut graph = MolGraph::new();
    let mut current: Option<usize> = None;
    let mut branch_stack: Vec<usize> = Vec::new();
    let mut pending_bond: Option<BondKind> = None;
    let mut ring_opens: HashMap<u8, (usize, Option<BondKind>)> = HashMap::new();

    for (tok_idx, token) in tokens.iter().enumerate() {
        match token {
            Token::Organic(oa) => {
                let atom = atom_from_organic(oa);
                let idx = graph.add_atom(atom);
                if let Some(curr) = current {
                    let bt = resolve_bond(pending_bond.take(), &graph.atoms[curr], &graph.atoms[idx]);
                    graph.add_bond(curr, idx, bt);
                }
                current = Some(idx);
            }
            Token::Aromatic(aa) => {
                let atom = atom_from_aromatic(aa);
                let idx = graph.add_atom(atom);
                if let Some(curr) = current {
                    let bt = resolve_bond(pending_bond.take(), &graph.atoms[curr], &graph.atoms[idx]);
                    graph.add_bond(curr, idx, bt);
                }
                current = Some(idx);
            }
            Token::Bracket(ba) => {
                let atom = atom_from_bracket(ba);
                let idx = graph.add_atom(atom);
                if let Some(curr) = current {
                    let bt = resolve_bond(pending_bond.take(), &graph.atoms[curr], &graph.atoms[idx]);
                    graph.add_bond(curr, idx, bt);
                }
                current = Some(idx);
            }
            Token::Bond(bk) => { pending_bond = Some(*bk); }
            Token::BranchOpen => {
                let curr = current.ok_or_else(|| GraphError { msg: format!("branch open at token {} with no current atom", tok_idx) })?;
                branch_stack.push(curr);
            }
            Token::BranchClose => {
                current = Some(branch_stack.pop().ok_or_else(|| GraphError { msg: format!("unmatched branch close at token {}", tok_idx) })?);
                pending_bond = None;
            }
            Token::RingBond(d) | Token::RingBond2(d) => {
                let curr = current.ok_or_else(|| GraphError { msg: format!("ring digit at token {} with no current atom", tok_idx) })?;
                let ring_id = *d;
                if let Some((open_atom, open_bond)) = ring_opens.remove(&ring_id) {
                    let bt = resolve_ring_bond(open_bond, pending_bond.take(), &graph.atoms[open_atom], &graph.atoms[curr])?;
                    graph.add_bond(open_atom, curr, bt);
                    graph.ring_bonds.push((open_atom, curr));
                } else {
                    ring_opens.insert(ring_id, (curr, pending_bond.take()));
                }
            }
            Token::Dot => { current = None; pending_bond = None; }
        }
    }

    if !ring_opens.is_empty() {
        let unclosed: Vec<u8> = ring_opens.keys().copied().collect();
        return Err(GraphError { msg: format!("unclosed ring(s): {:?}", unclosed) });
    }
    if !branch_stack.is_empty() {
        return Err(GraphError { msg: format!("{} unclosed branch(es)", branch_stack.len()) });
    }
    for i in 0..graph.atoms.len() {
        let bos = graph.bond_order_sum(i);
        graph.atoms[i].implicit_hcount = compute_implicit_h(&graph.atoms[i], bos);
    }
    Ok(graph)
}

/// Convenience: tokenize a SMILES string and build its molecular graph.
pub fn build_from_smiles(smiles: &str) -> Result<MolGraph, String> {
    let tokens = tokenize(smiles).map_err(|e| format!("tokenize: {}", e))?;
    build_graph(&tokens).map_err(|e| format!("graph: {}", e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(smiles: &str) -> MolGraph {
        build_from_smiles(smiles).unwrap_or_else(|e| panic!("failed on '{}': {}", smiles, e))
    }

    // ── Atom count ──
    #[test] fn methane_atom_count() { assert_eq!(graph("C").atom_count(), 1); }
    #[test] fn ethane_atom_count() { assert_eq!(graph("CC").atom_count(), 2); }
    #[test] fn ethanol_atom_count() { assert_eq!(graph("CCO").atom_count(), 3); }
    #[test] fn benzene_atom_count() { assert_eq!(graph("c1ccccc1").atom_count(), 6); }
    #[test] fn naphthalene_atom_count() { assert_eq!(graph("c1ccc2ccccc2c1").atom_count(), 10); }

    // ── Element identification ──
    #[test] fn methane_element() { let g = graph("C"); assert_eq!(g.atoms[0].element, "c"); assert!(!g.atoms[0].aromatic); }
    #[test] fn benzene_elements_aromatic() { let g = graph("c1ccccc1"); for a in &g.atoms { assert_eq!(a.element, "c"); assert!(a.aromatic); } }
    #[test]
    fn chloroform_elements() {
        let g = graph("ClC(Cl)Cl");
        assert_eq!(g.atom_count(), 4);
        assert_eq!(g.atoms[0].element, "cl"); assert_eq!(g.atoms[1].element, "c");
        assert_eq!(g.atoms[2].element, "cl"); assert_eq!(g.atoms[3].element, "cl");
    }
    #[test]
    fn pyridine_elements() {
        let g = graph("c1ccncc1");
        let elts: Vec<&str> = g.atoms.iter().map(|a| a.element.as_str()).collect();
        assert_eq!(elts, vec!["c","c","c","n","c","c"]);
        assert!(g.atoms.iter().all(|a| a.aromatic));
    }

    // ── Bond types ──
    #[test] fn ethane_single_bond() { assert_eq!(graph("CC").bond_between(0, 1), Some(BondType::Single)); }
    #[test] fn ethene_double_bond() { assert_eq!(graph("C=C").bond_between(0, 1), Some(BondType::Double)); }
    #[test] fn acetylene_triple_bond() { assert_eq!(graph("C#C").bond_between(0, 1), Some(BondType::Triple)); }
    #[test]
    fn benzene_aromatic_bonds() {
        let g = graph("c1ccccc1");
        for i in 0..6 { assert_eq!(g.bond_between(i, (i+1)%6), Some(BondType::Aromatic), "bond {}-{}", i, (i+1)%6); }
    }
    #[test]
    fn explicit_single_between_aromatic_atoms() {
        let g = graph("c1ccccc1-c2ccccc2"); // biphenyl
        assert_eq!(g.atom_count(), 12);
        assert_eq!(g.bond_between(5, 6), Some(BondType::Single));
    }
    #[test]
    fn implicit_bond_aromatic_to_nonaromatic() {
        let g = graph("c1ccccc1C"); // toluene
        assert_eq!(g.bond_between(5, 6), Some(BondType::Single));
    }

    // ── Branching ──
    #[test]
    fn acetic_acid_branching() {
        let g = graph("CC(=O)O");
        assert_eq!(g.atom_count(), 4);
        assert_eq!(g.bond_between(0, 1), Some(BondType::Single));
        assert_eq!(g.bond_between(1, 2), Some(BondType::Double));
        assert_eq!(g.bond_between(1, 3), Some(BondType::Single));
        assert_eq!(g.bond_between(0, 2), None);
    }
    #[test]
    fn nested_branches() {
        let g = graph("C(C)(C)C");
        assert_eq!(g.degree(0), 3);
        for i in 1..=3 { assert_eq!(g.bond_between(0, i), Some(BondType::Single)); }
    }
    #[test]
    fn branch_restores_current_atom() {
        let g = graph("CC(O)N");
        assert_eq!(g.bond_between(1, 2), Some(BondType::Single));
        assert_eq!(g.bond_between(1, 3), Some(BondType::Single));
        assert_eq!(g.bond_between(0, 3), None);
    }

    // ── Ring closures ──
    #[test]
    fn cyclohexane_ring() {
        let g = graph("C1CCCCC1");
        assert_eq!(g.atom_count(), 6);
        assert_eq!(g.bond_between(0, 5), Some(BondType::Single));
        assert_eq!(g.ring_bonds, vec![(0, 5)]);
        for i in 0..6 { assert_eq!(g.degree(i), 2, "atom {}", i); }
    }
    #[test] fn naphthalene_two_rings() { assert_eq!(graph("c1ccc2ccccc2c1").ring_bonds.len(), 2); }
    #[test]
    fn two_digit_ring_closure() {
        let g = graph("C%10CC%10");
        assert_eq!(g.ring_bonds.len(), 1);
        assert_eq!(g.bond_between(0, 2), Some(BondType::Single));
    }
    #[test] fn ring_closure_aromatic_inference() { assert_eq!(graph("c1ccccc1").bond_between(0, 5), Some(BondType::Aromatic)); }


    // ── Implicit hydrogens ──

    #[test]
    fn methane_4h() {
        // C has valence 4, degree 0 → 4 implicit H
        let g = graph("C");
        assert_eq!(g.atoms[0].implicit_hcount, 4);
    }

    #[test]
    fn ethane_3h_each() {
        // Each C has valence 4, 1 bond → 3 implicit H
        let g = graph("CC");
        assert_eq!(g.atoms[0].implicit_hcount, 3);
        assert_eq!(g.atoms[1].implicit_hcount, 3);
    }

    #[test]
    fn ethanol_hydrogen_counts() {
        // C-C-O: C(3H), C(2H), O(1H)
        let g = graph("CCO");
        assert_eq!(g.atoms[0].implicit_hcount, 3); // methyl
        assert_eq!(g.atoms[1].implicit_hcount, 2); // methylene
        assert_eq!(g.atoms[2].implicit_hcount, 1); // hydroxyl
    }

    #[test]
    fn acetic_acid_hydrogen_counts() {
        // CC(=O)O: CH3, C(0H), =O(0H), OH(1H)
        let g = graph("CC(=O)O");
        assert_eq!(g.atoms[0].implicit_hcount, 3); // methyl
        assert_eq!(g.atoms[1].implicit_hcount, 0); // carboxyl C: 1+2+1=4
        assert_eq!(g.atoms[2].implicit_hcount, 0); // =O
        assert_eq!(g.atoms[3].implicit_hcount, 1); // -OH
    }

    #[test]
    fn benzene_aromatic_hydrogen() {
        // Each aromatic C: 2 aromatic bonds (order 1 each) → sum=2
        // Valence 4, base_h = 4-2 = 2, minus 1 for aromatic → 1H
        let g = graph("c1ccccc1");
        for (i, a) in g.atoms.iter().enumerate() {
            assert_eq!(a.implicit_hcount, 1, "atom {} should have 1H", i);
        }
    }

    #[test]
    fn pyridine_nitrogen_0h() {
        // c1ccncc1: N has 2 aromatic bonds, valence 3
        // base_h = 3-2 = 1, minus 1 for aromatic → 0H
        let g = graph("c1ccncc1");
        let n_idx = g.atoms.iter().position(|a| a.element == "n").unwrap();
        assert_eq!(g.atoms[n_idx].implicit_hcount, 0);
    }

    #[test]
    fn ethene_0h_double_bonded() {
        // C=C: each C has 1 double bond (order 2), valence 4 → 2H each
        let g = graph("C=C");
        assert_eq!(g.atoms[0].implicit_hcount, 2);
        assert_eq!(g.atoms[1].implicit_hcount, 2);
    }

    #[test]
    fn hcn_hydrogen_counts() {
        // C#N: C has triple bond (3) → valence 4, 1H; N has triple bond (3) → valence 3, 0H
        let g = graph("C#N");
        assert_eq!(g.atoms[0].implicit_hcount, 1); // HC≡
        assert_eq!(g.atoms[1].implicit_hcount, 0); // ≡N
    }

    #[test]
    fn fluorine_1h() {
        let g = graph("F");
        assert_eq!(g.atoms[0].implicit_hcount, 1);
    }

    #[test]
    fn chlorine_1h() {
        let g = graph("Cl");
        assert_eq!(g.atoms[0].implicit_hcount, 1);
    }

    #[test]
    fn water_2h() {
        let g = graph("O");
        assert_eq!(g.atoms[0].implicit_hcount, 2);
    }

    #[test]
    fn ammonia_3h() {
        let g = graph("N");
        assert_eq!(g.atoms[0].implicit_hcount, 3);
    }

    #[test]
    fn sulfur_multi_valence() {
        // S alone: valence 2 → 2H
        assert_eq!(graph("S").atoms[0].implicit_hcount, 2);
        // S(=O): bond sum = 2, valence 2 → 0H
        // But wait — S has valences [2, 4, 6], bond sum 2 fits valence 2 → 0H
        let g = graph("S=O");
        assert_eq!(g.atoms[0].implicit_hcount, 0);
    }

    #[test]
    fn cyclohexane_2h_each() {
        let g = graph("C1CCCCC1");
        for (i, a) in g.atoms.iter().enumerate() {
            assert_eq!(a.implicit_hcount, 2, "atom {} should have 2H", i);
        }
    }

    // ── Bracket atoms ──

    #[test]
    fn bracket_atom_explicit_hcount() {
        // [NH4+] → explicit hcount = 4, charge = +1
        let g = graph("[NH4+]");
        assert_eq!(g.atoms[0].element, "n");
        assert_eq!(g.atoms[0].explicit_hcount, Some(4));
        assert_eq!(g.atoms[0].implicit_hcount, 4); // uses explicit
        assert_eq!(g.atoms[0].charge, 1);
    }

    #[test]
    fn bracket_atom_isotope() {
        let g = graph("[13C]");
        assert_eq!(g.atoms[0].isotope, Some(13));
        assert_eq!(g.atoms[0].element, "c");
    }

    #[test]
    fn bracket_atom_no_h_means_zero() {
        // [O-] → no H specified in bracket → 0 implicit H per OpenSMILES
        let g = graph("[O-]");
        assert_eq!(g.atoms[0].element, "o");
        assert_eq!(g.atoms[0].charge, -1);
        assert_eq!(g.atoms[0].explicit_hcount, Some(0));
        assert_eq!(g.atoms[0].implicit_hcount, 0);
    }

    #[test]
    fn pyrrole_nh() {
        // c1cc[nH]c1 — pyrrole: 5-membered ring, N has explicit H=1
        let g = graph("c1cc[nH]c1");
        assert_eq!(g.atom_count(), 5, "pyrrole should have 5 atoms, not {}", g.atom_count());
        let n_idx = g.atoms.iter().position(|a| a.element == "n").unwrap();
        assert_eq!(g.atoms[n_idx].explicit_hcount, Some(1));
        assert_eq!(g.atoms[n_idx].implicit_hcount, 1);
        assert!(g.atoms[n_idx].aromatic);
    }

    #[test]
    fn iron_bracket_no_h() {
        // [Fe] — bracket atom, no hcount field → 0 implicit H
        let g = graph("[Fe]");
        assert_eq!(g.atoms[0].element, "fe");
        assert_eq!(g.atoms[0].explicit_hcount, Some(0));
        assert_eq!(g.atoms[0].implicit_hcount, 0);
    }

    // ── Disconnected fragments (dot) ──

    #[test]
    fn dot_separates_fragments() {
        // [Na+].[Cl-] — two disconnected atoms
        let g = graph("[Na+].[Cl-]");
        assert_eq!(g.atom_count(), 2);
        assert_eq!(g.bond_between(0, 1), None);
        assert_eq!(g.atoms[0].charge, 1);
        assert_eq!(g.atoms[1].charge, -1);
    }

    #[test]
    fn dot_two_molecules() {
        // CC.CC — two ethane fragments
        let g = graph("CC.CC");
        assert_eq!(g.atom_count(), 4);
        assert_eq!(g.bond_between(0, 1), Some(BondType::Single));
        assert_eq!(g.bond_between(2, 3), Some(BondType::Single));
        assert_eq!(g.bond_between(1, 2), None);
    }

    // ── Degree / connectivity ──

    #[test]
    fn isobutane_central_degree() {
        // CC(C)C — central carbon bonded to 3 others
        let g = graph("CC(C)C");
        assert_eq!(g.degree(1), 3);
        assert_eq!(g.atoms[1].implicit_hcount, 1);
    }

    #[test]
    fn neopentane_degree_4() {
        // CC(C)(C)C — central C bonded to 4 carbons
        let g = graph("CC(C)(C)C");
        assert_eq!(g.degree(1), 4);
        assert_eq!(g.atoms[1].implicit_hcount, 0);
    }

    // ── Larger molecules ──

    #[test]
    fn toluene_structure() {
        // Cc1ccccc1 — methyl on benzene
        let g = graph("Cc1ccccc1");
        assert_eq!(g.atom_count(), 7);
        assert!(!g.atoms[0].aromatic); // methyl C
        assert!(g.atoms[1].aromatic);  // ring C
        assert_eq!(g.bond_between(0, 1), Some(BondType::Single));
        assert_eq!(g.atoms[0].implicit_hcount, 3); // CH3
    }

    #[test]
    fn benzoic_acid() {
        // OC(=O)c1ccccc1
        let g = graph("OC(=O)c1ccccc1");
        assert_eq!(g.atom_count(), 9); // O, C, =O, 6 ring C
        assert_eq!(g.bond_between(1, 2), Some(BondType::Double));
        assert_eq!(g.bond_between(1, 3), Some(BondType::Single)); // C-c(ring)
    }

    #[test]
    fn furan() {
        // c1ccoc1
        let g = graph("c1ccoc1");
        assert_eq!(g.atom_count(), 5);
        let o_idx = g.atoms.iter().position(|a| a.element == "o").unwrap();
        assert!(g.atoms[o_idx].aromatic);
        assert_eq!(g.atoms[o_idx].implicit_hcount, 0); // O in aromatic ring: val 2, bos 2, -1 aromatic → clamp 0
    }

    // ── Error cases ──

    #[test]
    fn unclosed_ring_error() {
        let result = build_from_smiles("C1CC");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("unclosed ring"));
    }

    #[test]
    fn unmatched_branch_close_error() {
        let result = build_from_smiles("CC)");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("unmatched branch close"));
    }

    #[test]
    fn unclosed_branch_error() {
        let result = build_from_smiles("C(C");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("unclosed branch"));
    }

    #[test]
    fn conflicting_ring_bond_error() {
        // C=1CC-1 → ring open says double, close says single
        let result = build_from_smiles("C=1CC-1");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("conflicting"));
    }
}
