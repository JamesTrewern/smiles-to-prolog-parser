use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt;
use crate::graph::*;

// ── Bond symbol for Prolog output ──

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondSymbol {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondSymbol {
    pub fn from_bond_type(bt: BondType) -> Self {
        match bt {
            BondType::Single => BondSymbol::Single,
            BondType::Double => BondSymbol::Double,
            BondType::Triple => BondSymbol::Triple,
            BondType::Aromatic => BondSymbol::Aromatic,
        }
    }
}

impl fmt::Display for BondSymbol {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            BondSymbol::Single => write!(f, "1"),
            BondSymbol::Double => write!(f, "2"),
            BondSymbol::Triple => write!(f, "3"),
            BondSymbol::Aromatic => write!(f, "a"),
        }
    }
}

// ── Prolog term AST ──

/// A complete molecule term: mol(Element, Subs) or mol(ring(Members))
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolTerm {
    AtomRoot { element: String, subs: Vec<SubTerm> },
    RingRoot { members: Vec<RingMember> },
}

/// A sub-term appearing in a substituent list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubTerm {
    /// Bare `h` shorthand
    H,
    /// (bond, element, subs)
    Atom { bond: BondSymbol, element: String, subs: Vec<SubTerm> },
    /// (bond, ring(members))
    Ring { bond: BondSymbol, members: Vec<RingMember> },
}

/// A member within a ring(...) term.
/// Member 0's bond = closing bond (last→first).
/// Member k's bond (k>0) = bond from member k-1 to k.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingMember {
    pub bond: BondSymbol,
    pub element: String,
    pub subs: Vec<SubTerm>,
}

// ── Subtree complexity for ordering ──

/// Complexity metric: (depth, heavy_atom_count, total_atom_count).
/// Compared lexicographically. Higher = more complex = goes later in lists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Complexity(usize, usize, usize);

impl Complexity {
    fn hydrogen() -> Self { Complexity(0, 0, 1) }
}

fn subterm_complexity(t: &SubTerm) -> Complexity {
    match t {
        SubTerm::H => Complexity::hydrogen(),
        SubTerm::Atom { subs, .. } => {
            let child_max = subs.iter().map(|s| subterm_complexity(s))
                .max().unwrap_or(Complexity(0, 0, 0));
            let (h, t) = count_atoms_subterm_list(subs);
            Complexity(child_max.0 + 1, h + 1, t + 1) // +1 for this atom
        }
        SubTerm::Ring { members, .. } => {
            let mut max_depth = 0usize;
            let mut heavy = 0usize;
            let mut total = 0usize;
            for m in members {
                let mc = ring_member_complexity(m);
                if mc.0 > max_depth { max_depth = mc.0; }
                heavy += mc.1;
                total += mc.2;
            }
            Complexity(max_depth + 1, heavy, total)
        }
    }
}

fn ring_member_complexity(m: &RingMember) -> Complexity {
    let child_max = m.subs.iter().map(|s| subterm_complexity(s))
        .max().unwrap_or(Complexity(0, 0, 0));
    let (h, t) = count_atoms_subterm_list(&m.subs);
    Complexity(child_max.0 + 1, h + 1, t + 1)
}

/// Count (heavy_atoms, total_atoms) in a substituent list.
fn count_atoms_subterm_list(subs: &[SubTerm]) -> (usize, usize) {
    let mut heavy = 0usize;
    let mut total = 0usize;
    for s in subs {
        let c = subterm_complexity(s);
        heavy += c.1;
        total += c.2;
    }
    (heavy, total)
}

/// Atomic number for tiebreaking. Higher = heavier.
fn atomic_number(element: &str) -> u8 {
    match element {
        "h" => 1, "b" => 5, "c" => 6, "n" => 7, "o" => 8,
        "f" => 9, "p" => 15, "s" => 16, "cl" => 17, "br" => 35, "i" => 53,
        "fe" => 26, "na" => 11, "se" => 34, "as" => 33,
        _ => 0, // unknown
    }
}

/// Sort key for substituents: simplest first, most complex last.
/// H < simple atoms < complex subtrees < rings.
fn subterm_sort_key(t: &SubTerm) -> (u8, Complexity, u8) {
    match t {
        SubTerm::H => (0, Complexity::hydrogen(), 1),
        SubTerm::Atom { element, subs: _, .. } => {
            let c = subterm_complexity(t);
            (1, c, atomic_number(element))
        }
        SubTerm::Ring { .. } => {
            let c = subterm_complexity(t);
            (2, c, 0)
        }
    }
}

fn sort_subs(subs: &mut Vec<SubTerm>) {
    subs.sort_by(|a, b| subterm_sort_key(a).cmp(&subterm_sort_key(b)));
}

// ── Ring detection ──

/// A detected ring: ordered list of atom indices forming the cycle.
#[derive(Debug, Clone)]
struct DetectedRing {
    members: Vec<usize>,
    /// The ring closure bond (a, b) that defined this ring.
    #[allow(dead_code)]
    closure: (usize, usize),
}

/// Set of edges that are ring closure bonds (stored as sorted tuples).
fn ring_closure_set(graph: &MolGraph) -> HashSet<(usize, usize)> {
    graph.ring_bonds.iter().map(|&(a, b)| if a < b { (a, b) } else { (b, a) }).collect()
}

/// Is edge (a, b) a ring closure bond?
fn is_ring_closure(a: usize, b: usize, closures: &HashSet<(usize, usize)>) -> bool {
    let key = if a < b { (a, b) } else { (b, a) };
    closures.contains(&key)
}

/// Find the ring cycle for a given ring closure bond (src, dst).
/// BFS from src to dst using only tree edges (non-closure bonds).
/// The ring is: [path from src to dst via tree] which, combined with
/// the closure edge dst→src, forms the cycle.
fn find_ring_path(
    graph: &MolGraph,
    src: usize,
    dst: usize,
    closures: &HashSet<(usize, usize)>,
) -> Vec<usize> {
    // BFS from src, avoiding closure edges, to find shortest path to dst.
    let mut visited = HashSet::new();
    let mut parent: HashMap<usize, usize> = HashMap::new();
    let mut queue = VecDeque::new();
    visited.insert(src);
    queue.push_back(src);

    while let Some(node) = queue.pop_front() {
        for edge in graph.neighbors(node) {
            let nbr = edge.target;
            // Skip the direct closure edge between src and dst
            if (node == src && nbr == dst) || (node == dst && nbr == src) {
                if is_ring_closure(node, nbr, closures) {
                    continue;
                }
            }
            // Skip other closure edges
            if is_ring_closure(node, nbr, closures) && !(node == src || node == dst || nbr == src || nbr == dst) {
                continue;
            }
            if visited.contains(&nbr) { continue; }
            visited.insert(nbr);
            parent.insert(nbr, node);
            if nbr == dst {
                // Reconstruct path
                let mut path = vec![dst];
                let mut cur = dst;
                while cur != src {
                    cur = parent[&cur];
                    path.push(cur);
                }
                path.reverse();
                return path;
            }
            queue.push_back(nbr);
        }
    }
    // Should not happen for a valid ring closure
    vec![src, dst]
}

/// Detect all rings in the graph.
fn detect_rings(graph: &MolGraph) -> Vec<DetectedRing> {
    let closures = ring_closure_set(graph);
    graph.ring_bonds.iter().map(|&(a, b)| {
        let members = find_ring_path(graph, a, b, &closures);
        DetectedRing { members, closure: (a, b) }
    }).collect()
}

/// Build a map: atom_index → list of ring indices it belongs to.
fn atom_ring_membership(rings: &[DetectedRing]) -> HashMap<usize, Vec<usize>> {
    let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ri, ring) in rings.iter().enumerate() {
        for &atom in &ring.members {
            map.entry(atom).or_default().push(ri);
        }
    }
    map
}

// ── Root selection ──

/// Compute subtree complexity for an atom in the graph, used for root selection.
/// This is a rough metric: count heavy atoms reachable from `start`
/// without crossing into ring members (except through start itself).
fn graph_subtree_size(
    graph: &MolGraph,
    start: usize,
    exclude_ring_members: &HashSet<usize>,
) -> usize {
    let mut visited = HashSet::new();
    let mut stack = vec![start];
    let mut count = 0;
    while let Some(node) = stack.pop() {
        if !visited.insert(node) { continue; }
        count += 1;
        for edge in graph.neighbors(node) {
            if visited.contains(&edge.target) { continue; }
            if node != start && exclude_ring_members.contains(&edge.target) { continue; }
            stack.push(edge.target);
        }
    }
    count
}

/// Select which ring should be the root (if any rings exist).
/// Pick the ring whose members have the most complex external subtrees.
fn select_root_ring(graph: &MolGraph, rings: &[DetectedRing]) -> usize {
    if rings.len() == 1 { return 0; }
    let mut best = 0;
    let mut best_score = 0usize;
    for (ri, ring) in rings.iter().enumerate() {
        let ring_set: HashSet<usize> = ring.members.iter().copied().collect();
        let mut score = 0;
        for &atom in &ring.members {
            score += graph_subtree_size(graph, atom, &ring_set);
        }
        if score > best_score {
            best_score = score;
            best = ri;
        }
    }
    best
}

/// Select root atom for a molecule with no rings.
/// Highest non-H degree, then atomic number descending.
fn select_root_atom(graph: &MolGraph) -> usize {
    (0..graph.atom_count())
        .max_by_key(|&i| {
            let deg = graph.degree(i);
            let an = atomic_number(&graph.atoms[i].element);
            (deg, an)
        })
        .unwrap_or(0)
}

// ── Term generation (DFS) ──

/// Generate a MolTerm from a MolGraph.
pub fn generate_term(graph: &MolGraph) -> MolTerm {
    let rings = detect_rings(graph);
    let membership = atom_ring_membership(&rings);
    let closures = ring_closure_set(graph);

    // Track which atoms have been consumed
    let mut visited: HashSet<usize> = HashSet::new();

    if rings.is_empty() {
        // No rings — atom root
        let root = select_root_atom(graph);
        let (element, subs) = build_atom_term(
            graph, root, &rings, &membership, &closures, &mut visited,
        );
        MolTerm::AtomRoot { element, subs }
    } else {
        // Ring root
        let root_ring_idx = select_root_ring(graph, &rings);
        let members = build_ring_term(
            graph, root_ring_idx, &rings, &membership, &closures, &mut visited,
        );
        MolTerm::RingRoot { members }
    }
}

/// Build the sub-terms for an atom (non-ring context).
/// Returns (element, sorted substituent list).
fn build_atom_term(
    graph: &MolGraph,
    atom_idx: usize,
    rings: &[DetectedRing],
    membership: &HashMap<usize, Vec<usize>>,
    closures: &HashSet<(usize, usize)>,
    visited: &mut HashSet<usize>,
) -> (String, Vec<SubTerm>) {
    visited.insert(atom_idx);
    let atom = &graph.atoms[atom_idx];
    let mut subs = Vec::new();

    // Add implicit hydrogens
    for _ in 0..atom.implicit_hcount {
        subs.push(SubTerm::H);
    }

    // Process neighbors
    for edge in graph.neighbors(atom_idx) {
        let nbr = edge.target;
        if visited.contains(&nbr) { continue; }

        let bond = BondSymbol::from_bond_type(edge.bond_type);

        // Is neighbor part of a ring we haven't yet consumed?
        if let Some(ring_indices) = membership.get(&nbr) {
            // Find a ring that contains nbr and hasn't been fully visited
            let ring_to_build = ring_indices.iter().find(|&&ri| {
                rings[ri].members.iter().any(|&m| !visited.contains(&m))
            });
            if let Some(&ri) = ring_to_build {
                let members = build_ring_term(
                    graph, ri, rings, membership, closures, visited,
                );
                subs.push(SubTerm::Ring { bond, members });
                continue;
            }
        }

        // Regular atom neighbor
        let (el, child_subs) = build_atom_term(
            graph, nbr, rings, membership, closures, visited,
        );
        subs.push(SubTerm::Atom { bond, element: el, subs: child_subs });
    }

    sort_subs(&mut subs);
    (atom.element.clone(), subs)
}

/// Build the members list for a ring.
fn build_ring_term(
    graph: &MolGraph,
    ring_idx: usize,
    rings: &[DetectedRing],
    membership: &HashMap<usize, Vec<usize>>,
    closures: &HashSet<(usize, usize)>,
    visited: &mut HashSet<usize>,
) -> Vec<RingMember> {
    let ring = &rings[ring_idx];
    let ring_set: HashSet<usize> = ring.members.iter().copied().collect();

    // Mark all ring members as visited
    for &m in &ring.members {
        visited.insert(m);
    }

    // Build raw members with their non-ring substituents
    let mut raw_members: Vec<(usize, RingMember)> = Vec::new();
    let n = ring.members.len();

    for (pos, &atom_idx) in ring.members.iter().enumerate() {
        let atom = &graph.atoms[atom_idx];

        // Determine the bond for this member
        let bond = if pos == 0 {
            // Closing bond: from last member back to first
            let last = ring.members[n - 1];
            let bt = graph.bond_between(last, atom_idx).unwrap_or(BondType::Single);
            BondSymbol::from_bond_type(bt)
        } else {
            let prev = ring.members[pos - 1];
            let bt = graph.bond_between(prev, atom_idx).unwrap_or(BondType::Single);
            BondSymbol::from_bond_type(bt)
        };

        // Build non-ring substituents for this ring member
        let mut subs = Vec::new();
        for _ in 0..atom.implicit_hcount {
            subs.push(SubTerm::H);
        }
        for edge in graph.neighbors(atom_idx) {
            let nbr = edge.target;
            if ring_set.contains(&nbr) { continue; } // skip ring neighbors
            if visited.contains(&nbr) { continue; }

            let bond_sym = BondSymbol::from_bond_type(edge.bond_type);

            // Check if neighbor starts a new ring
            if let Some(ring_indices) = membership.get(&nbr) {
                let ring_to_build = ring_indices.iter().find(|&&ri| {
                    ri != ring_idx && rings[ri].members.iter().any(|&m| !visited.contains(&m))
                });
                if let Some(&ri) = ring_to_build {
                    let members = build_ring_term(
                        graph, ri, rings, membership, closures, visited,
                    );
                    subs.push(SubTerm::Ring { bond: bond_sym, members });
                    continue;
                }
            }

            let (el, child_subs) = build_atom_term(
                graph, nbr, rings, membership, closures, visited,
            );
            subs.push(SubTerm::Atom { bond: bond_sym, element: el, subs: child_subs });
        }
        sort_subs(&mut subs);
        raw_members.push((atom_idx, RingMember { bond, element: atom.element.clone(), subs }));
    }

    // ── Ring member ordering ──
    // Find the member with the most complex non-ring subtree → goes first.
    // Complexity of a member's subs (excluding H for this comparison).
    let member_complexity: Vec<Complexity> = raw_members.iter()
        .map(|(_, m)| {
            m.subs.iter()
                .map(|s| subterm_complexity(s))
                .max()
                .unwrap_or(Complexity(0, 0, 0))
        })
        .collect();

    // Find the starting position (most complex subtree)
    let start_pos = member_complexity.iter().enumerate()
        .max_by(|(ia, ca), (ib, cb)| {
            ca.cmp(cb).then_with(|| {
                let an_a = atomic_number(&raw_members[*ia].1.element);
                let an_b = atomic_number(&raw_members[*ib].1.element);
                an_a.cmp(&an_b)
            })
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Determine traversal direction: prefer direction where the next
    // member (position start+1) is more complex than the previous (start-1).
    let n = raw_members.len();
    let fwd_next = (start_pos + 1) % n;
    let bwd_next = (start_pos + n - 1) % n;
    let fwd_c = &member_complexity[fwd_next];
    let bwd_c = &member_complexity[bwd_next];

    let go_forward = match fwd_c.cmp(bwd_c) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Equal => {
            // Tiebreak by atomic number
            let fwd_an = atomic_number(&raw_members[fwd_next].1.element);
            let bwd_an = atomic_number(&raw_members[bwd_next].1.element);
            fwd_an >= bwd_an
        }
    };

    // Build final ordered members
    let mut ordered = Vec::with_capacity(n);
    for i in 0..n {
        let idx = if go_forward {
            (start_pos + i) % n
        } else {
            (start_pos + n - i) % n
        };
        ordered.push(raw_members[idx].1.clone());
    }

    // When we reorder, we need to fix up the bond symbols.
    // Member 0's bond = closing bond (from last member to first).
    // Member k's bond = bond from member k-1 to member k.
    // After reordering, recalculate from the graph.
    let reordered_atoms: Vec<usize> = {
        let mut v = Vec::with_capacity(n);
        for i in 0..n {
            let idx = if go_forward {
                (start_pos + i) % n
            } else {
                (start_pos + n - i) % n
            };
            v.push(raw_members[idx].0); // atom index
        }
        v
    };

    for i in 0..n {
        let bond = if i == 0 {
            // Closing bond: last → first
            let last = reordered_atoms[n - 1];
            let first = reordered_atoms[0];
            let bt = graph.bond_between(last, first).unwrap_or(BondType::Single);
            BondSymbol::from_bond_type(bt)
        } else {
            let prev = reordered_atoms[i - 1];
            let curr = reordered_atoms[i];
            let bt = graph.bond_between(prev, curr).unwrap_or(BondType::Single);
            BondSymbol::from_bond_type(bt)
        };
        ordered[i].bond = bond;
    }

    ordered
}

// ── Display (Prolog syntax) ──

/// Format a substituent list as Prolog list syntax: [a, b, c]
fn fmt_subs(subs: &[SubTerm], f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "[")?;
    for (i, s) in subs.iter().enumerate() {
        if i > 0 { write!(f, ", ")?; }
        write!(f, "{}", s)?;
    }
    write!(f, "]")
}

impl fmt::Display for MolTerm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            MolTerm::AtomRoot { element, subs } => {
                write!(f, "mol({}, ", element)?;
                fmt_subs(subs, f)?;
                write!(f, ")")
            }
            MolTerm::RingRoot { members } => {
                write!(f, "mol(ring([")?;
                for (i, m) in members.iter().enumerate() {
                    if i > 0 { write!(f, ", ")?; }
                    write!(f, "{}", m)?;
                }
                write!(f, "]))")
            }
        }
    }
}

impl fmt::Display for SubTerm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SubTerm::H => write!(f, "h"),
            SubTerm::Atom { bond, element, subs } => {
                write!(f, "({}, {}, ", bond, element)?;
                fmt_subs(subs, f)?;
                write!(f, ")")
            }
            SubTerm::Ring { bond, members } => {
                write!(f, "({}, ring([", bond)?;
                for (i, m) in members.iter().enumerate() {
                    if i > 0 { write!(f, ", ")?; }
                    write!(f, "{}", m)?;
                }
                write!(f, "]))")
            }
        }
    }
}

impl fmt::Display for RingMember {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, ", self.bond, self.element)?;
        fmt_subs(&self.subs, f)?;
        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::super::*;

    fn prolog(smiles: &str) -> String {
        smiles_to_prolog(smiles).unwrap_or_else(|e| panic!("failed on '{}': {}", smiles, e))
    }

    // ── Simple molecules (no rings) ──

    #[test]
    fn methane() {
        assert_eq!(prolog("C"), "mol(c, [h, h, h, h])");
    }

    #[test]
    fn water() {
        assert_eq!(prolog("O"), "mol(o, [h, h])");
    }

    #[test]
    fn hf() {
        // H-F: F has 1 implicit H. Single atom.
        assert_eq!(prolog("F"), "mol(f, [h])");
    }

    #[test]
    fn ethane() {
        // CC: root should be either C (both degree 1). Both identical,
        // so either is fine. Each has 3H + 1 bond to other.
        let s = prolog("CC");
        // Root C bonded to other C. Both have 3H.
        assert_eq!(s, "mol(c, [h, h, h, (1, c, [h, h, h])])");
    }

    #[test]
    fn ethanol() {
        // CCO: root should be middle C (degree 2), or C bonded to O
        // Actually root = highest degree. All are degree 1 or 2.
        // C(idx0) deg=1, C(idx1) deg=2, O(idx2) deg=1
        // Root = C(idx1) with degree 2
        let s = prolog("CCO");
        // Root is middle C: bonded to CH3 and OH
        // Subs ordered simplest first: h, h before the complex ones
        // CH3 complexity: depth=1, heavy=1, total=4
        // OH complexity: depth=1, heavy=1, total=2
        // OH is simpler (fewer total atoms), so OH before CH3
        assert_eq!(s, "mol(c, [h, h, (1, o, [h]), (1, c, [h, h, h])])");
    }

    #[test]
    fn acetic_acid() {
        // CC(=O)O: carboxyl C (idx1) has degree 3, should be root
        let s = prolog("CC(=O)O");
        // Root = C(idx1), subs: CH3(single), =O(double), -OH(single)
        // =O: depth=1, heavy=1, total=1 (no H on =O)
        // -OH: depth=1, heavy=1, total=2
        // CH3: depth=1, heavy=1, total=4
        // Order: =O(simplest), -OH, CH3(most complex)
        assert_eq!(s, "mol(c, [(2, o, []), (1, o, [h]), (1, c, [h, h, h])])");
    }

    #[test]
    fn formaldehyde() {
        // C=O: C has degree 1 (to O), O has degree 1 (to C)
        // Tiebreak by atomic number: O(8) > C(6), so O is root
        let s = prolog("C=O");
        assert_eq!(s, "mol(o, [(2, c, [h, h])])");
    }

    #[test]
    fn hcn() {
        // C#N: both degree 1. N(7) > C(6) by atomic number → N is root
        let s = prolog("C#N");
        assert_eq!(s, "mol(n, [(3, c, [h])])");
    }

    // ── Ring molecules ──

    #[test]
    fn benzene() {
        // c1ccccc1: all symmetric, all aromatic C with 1H
        let s = prolog("c1ccccc1");
        // All members identical: 6 × (a, c, [h])
        assert_eq!(s, "mol(ring([(a, c, [h]), (a, c, [h]), (a, c, [h]), (a, c, [h]), (a, c, [h]), (a, c, [h])]))");
    }

    #[test]
    fn cyclohexane() {
        let s = prolog("C1CCCCC1");
        assert_eq!(s, "mol(ring([(1, c, [h, h]), (1, c, [h, h]), (1, c, [h, h]), (1, c, [h, h]), (1, c, [h, h]), (1, c, [h, h])]))");
    }

    #[test]
    fn pyridine() {
        // c1ccncc1: ring with nitrogen
        let s = prolog("c1ccncc1");
        // N has 0H and highest atomic number (7) but has simplest subtree (empty).
        // All C have [h] subtrees (depth=0,heavy=0,total=1).
        // N has [] subtree (depth=0,heavy=0,total=0) — simplest.
        // Most complex member = any C with [h], all tied.
        // So ordering by atomic number: C(6) < N(7).
        // The member with most complex subtree is any C (they all have h).
        // N goes... depends on implementation.
        assert!(s.contains("ring("));
        assert!(s.contains("(a, n, [])"));
    }
}
