# Representing Organic Chemistry as Prolog Terms

## A Design for SMILES-to-Prolog Conversion for Meta-Interpretive Learning

---

## 1. Motivation

Inductive Logic Programming (ILP) has a long history of application to molecular structure-activity relationship (SAR) problems, dating back to the work of King, Muggleton, Srinivasan and Sternberg in the 1990s. In that tradition, molecules are represented as collections of ground Prolog facts — `atom(mol1, a1, c, 22, 0.015)`, `bond(mol1, a1, a2, 7)` — which describe the molecular graph explicitly. While effective, this approach is verbose (hundreds of facts per molecule) and requires the learner to reconstruct higher-level chemical concepts like ring systems and functional groups from individual atom-bond relations.

SMILES (Simplified Molecular Input Line Entry System) is the dominant string encoding for molecular structure in cheminformatics. It encodes the molecular graph as a depth-first traversal, using parentheses for branching and numeric labels for ring closures. This document describes a design for converting SMILES strings into nested Prolog terms that preserve the hierarchical structure of a molecule in a form amenable to Meta-Interpretive Learning (MIL).

The goal is a representation where:

- Chemical substructures are directly visible as subterms.
- A MIL system can learn structural patterns by matching on term structure and inventing predicates that compose features.
- The representation is compact and uniform, without requiring per-molecule ground fact databases.
- Substituent lists enable uniform traversal predicates without pattern-matching on every possible arity.

---

## 2. Core Representation

### 2.1 The Term Forms

The representation uses a small set of term shapes:

| Form | Structure | Meaning |
|------|-----------|---------|
| **Molecule root (atom)** | `mol(Element, Subs)` | Complete molecule rooted at an atom |
| **Molecule root (ring)** | `mol(ring(Members))` | Complete molecule whose root is a ring |
| **Sub-atom** | `(BondType, Element, Subs)` | Atom within a molecule, bonded to parent |
| **Sub-ring** | `(BondType, ring(Members))` | Ring within a molecule, bonded to parent |
| **Hydrogen shorthand** | `h` | Equivalent to `(1, h, [])` |

**`mol`** marks a complete molecule. It has no bond type because there is no parent — this is the top-level entry point. When the root is a plain atom, it is `mol(Element, Subs)`. When the root is a ring system, it is `mol(ring(Members))`.

**Sub-atoms** are 3-element tuples `(BondType, Element, Subs)` where:

- **BondType** describes the bond from the parent to this atom:
  - `1` — single bond
  - `2` — double bond
  - `3` — triple bond
  - `a` — aromatic (delocalised) bond
- **Element** is the atomic symbol: `c`, `n`, `o`, `s`, `p`, `f`, `cl`, `br`, `i`, etc.
- **Subs** is a **list** of substituents (child atoms, hydrogen shorthand, or ring terms).

**Sub-rings** are 2-element tuples `(BondType, ring(Members))`. The bond type describes the bond from the parent atom to the **first** member of the ring. Sub-rings have no substituent list — all substituents on ring atoms live inside the individual ring members' own lists.

**Hydrogen shorthand**: Bare `h` appearing in a substituent list is shorthand for `(1, h, [])`.

### 2.2 Ring Structure

All atoms of a ring live inside the `ring(Members)` term. The ring replaces what would normally be an element symbol, appearing either at the root (`mol(ring(...))`) or as a child (`(Bond, ring(...))`).

```
ring([(B1, E1, Subs1), (B2, E2, Subs2), ..., (BN, EN, SubsN)])
```

Where:

- **Member 1's BondType (B1)** describes the **closing bond** — the bond from the last member back to the first, completing the cycle.
- **Member K's BondType (K > 1)** describes the bond from member K-1 to member K.
- Each member's **Subs** list contains that atom's non-ring substituents (hydrogens, functional groups, or further rings).

This makes the bond description fully cyclic: B1 closes the ring, B2 connects member 1→2, B3 connects 2→3, and so on.

When a ring appears as a sub-term `(BondType, ring(Members))`, the outer BondType describes the bond from the **parent atom** to the **first member** of the ring. This bond is separate from the closing bond B1 — the first member may have two different bonds connecting it to the rest of the molecule vs closing the ring.

### 2.3 Examples — Linear and Branched Molecules

**Methane** (CH₄) — `C` in SMILES:

```prolog
mol(c, [h, h, h, h])
```

**Ethanol** (CH₃CH₂OH) — `CCO` in SMILES:

```prolog
mol(c, [h, h, h, (1, c, [h, h, (1, o, [h])])])
```

**Acetic acid** (CH₃COOH) — `CC(=O)O` in SMILES:

```prolog
mol(c, [h, h, h, (1, c, [(2, o), (1, o, [h])])])
```

Here the carboxyl carbon has two substituents: a double-bonded oxygen `(2, o)` (shorthand for `(2, o, [])`) and a hydroxyl group `(1, o, [h])`.

**Isobutane** — `CC(C)C` in SMILES:

```prolog
mol(c, [h, h, h, (1, c, [h, (1, c, [h, h, h]), (1, c, [h, h, h])])])
```

### 2.4 Examples — Rings

**Benzene** (`c1ccccc1` in SMILES):

All six carbons are inside the ring. The ring is the root of the molecule:

```prolog
mol(ring([
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))
```

Every bond type is `a` (aromatic). The first member's `a` describes the closing bond (member 6 → member 1).

**Cyclohexane** (`C1CCCCC1` in SMILES):

```prolog
mol(ring([
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h])
]))
```

**Pyridine** (`c1ccncc1` in SMILES):

An aromatic ring with a nitrogen. The nitrogen has no hydrogen (its lone pair participates in the aromatic system):

```prolog
mol(ring([
    (a, c, [h]),
    (a, c, [h]),
    (a, n, []),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))
```

### 2.5 Examples — Substituted Rings

**Benzoic acid** (`OC(=O)c1ccccc1` in SMILES):

The carboxyl group hangs off one ring carbon. Rooting at the ring:

```prolog
mol(ring([
    (a, c, [(1, c, [(2, o), (1, o, [h])])]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))
```

The carboxyl carbon `(1, c, [(2, o), (1, o, [h])])` sits in the first ring member's substituent list.

Alternatively, rooting at the carboxyl carbon with the ring as a child:

```prolog
mol(c, [(2, o), (1, o, [h]), (a, ring([
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))])
```

Here `(a, ring([...]))` means "aromatic bond from the carboxyl carbon to the first ring member."

**Cyclohexanol** — hydroxyl on cyclohexane:

```prolog
mol(ring([
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, h]),
    (1, c, [h, (1, o, [h])])
]))
```

**4-hydroxybenzoic acid** — carboxyl at position 1, hydroxyl at position 4:

```prolog
mol(ring([
    (a, c, [(1, c, [(2, o), (1, o, [h])])]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [(1, o, [h])]),
    (a, c, [h]),
    (a, c, [h])
]))
```

### 2.6 Examples — Multiple Rings

**Biphenyl** (two benzene rings connected by a single bond):

```prolog
mol(ring([
    (a, c, [(1, ring([
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h])
    ]))]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))
```

The inner ring `(1, ring([...]))` sits in the first ring member's substituent list. The `1` means a single bond connects the two rings.

### 2.7 Leaf Atom Convention

When a sub-atom has no substituents, its list is empty: `(2, o, [])`. For brevity in inline examples, a leaf atom may be written as `(2, o)` where context is unambiguous, but the canonical form always includes the list.

---

## 3. Fused Rings

Fused ring systems (e.g., naphthalene, indole) share an edge between two rings. This is represented using an **`edge` functor** within a ring member's substituent list, which marks the shared atom and introduces the fused ring:

```
edge(SharedAtom, ring(Members))
```

The `edge` functor sits in the substituent list of the ring member at the junction point. The SharedAtom is the other atom shared between the two rings (the junction atom carrying the edge is implicit). The fused ring's member list describes the path through the second ring.

**Naphthalene** (`c1ccc2ccccc2c1` in SMILES):

```prolog
mol(ring([
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [edge(c, ring([
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h]),
        (a, c, [h])
    ]))]),
    (a, c, [h]),
    (a, c, [h]),
    (a, c, [h])
]))
```

*Note: The fused ring design needs further validation against complex polycyclic systems.*

---

## 4. Bond Types

| Symbol | Meaning | Bond Order |
|--------|---------|------------|
| `1` | Single bond | 1 |
| `2` | Double bond | 2 |
| `3` | Triple bond | 3 |
| `a` | Aromatic / delocalised | ~1.5 |

The aromatic bond type `a` represents delocalised electrons shared across the ring system. This is physically accurate — aromatic bonds genuinely have a bond order of approximately 1.5, rather than being strictly alternating single and double. Using `a` avoids the canonicalisation problem that arises with Kekulé structures, where the same molecule could be represented with two different alternating patterns.

Non-aromatic rings that have genuine alternating bonds (e.g., cyclooctatetraene, which fails Hückel's rule) use the appropriate integer bond types.

---

## 5. Canonical Ordering

A given molecule can have multiple valid SMILES representations depending on which atom the traversal starts from and which direction it proceeds. For the Prolog term representation to be unique per molecule, a **canonical ordering** is required.

### 5.1 Root Selection

The root node should be chosen deterministically:

1. If the molecule contains ring(s), the primary ring system becomes the root (`mol(ring(...))`).
2. If no rings, select the atom with the highest degree of branching.
3. Break ties by atomic number (heavier elements first), then by substituent complexity.

### 5.2 Ring Member Ordering

Within a ring, choose the starting member and traversal direction deterministically:

1. Start at the most-substituted member (most non-hydrogen substituents).
2. Break ties by atomic number, then substituent complexity.
3. Traverse in the direction that yields the lexicographically smaller sequence.

### 5.3 Substituent Ordering Within Lists

When an atom has multiple substituents in its list, order them by:

1. Hydrogens first (simplest).
2. Non-ring, non-hydrogen substituents next, ordered by atomic number then complexity.
3. Ring substituents last.

---

## 6. Background Knowledge for MIL

The list-based structure enables clean, uniform background predicates. Because substituents are always in lists, predicates can use standard list operations (`member/2`) rather than needing separate clauses for every possible arity. Rings require their own clauses since they have a different term shape.

### 6.1 Core Accessors

```prolog
% Element of a molecule root (atom only)
element(mol(E, _), E).
element((_, E, _), E).

% Bond type (sub-atoms and sub-rings)
bond_type((B, _, _), B).
bond_type((B, ring(_)), B).

% Substituent list (atoms only — rings don't have one)
subs(mol(_, S), S).
subs((_, _, S), S).

% Ring detection
is_ring(mol(ring(_))).
is_ring((_, ring(_))).

% Ring member access
ring_members(mol(ring(M)), M).
ring_members((_, ring(M)), M).
```

### 6.2 Recursive Traversal

```prolog
% Does molecule/subterm contain element E anywhere?
contains(mol(E, _), E).
contains(mol(_, Subs), E) :- member(S, Subs), contains(S, E).
contains(mol(ring(Members)), E) :- member((_, E, _), Members).
contains(mol(ring(Members)), E) :-
    member((_, _, Subs), Members), member(S, Subs), contains(S, E).
contains((_, E, _), E).
contains((_, _, Subs), E) :- member(S, Subs), contains(S, E).
contains((_, ring(Members)), E) :- member((_, E, _), Members).
contains((_, ring(Members)), E) :-
    member((_, _, Subs), Members), member(S, Subs), contains(S, E).
```

### 6.3 Functional Group Patterns

```prolog
% Carbonyl: carbon double-bonded to oxygen
carbonyl(mol(c, Subs)) :- member((2, o, _), Subs).
carbonyl((_, c, Subs)) :- member((2, o, _), Subs).

% Carboxyl: carbon with both C=O and C-OH
carboxyl(X) :-
    element(X, c),
    subs(X, S),
    member((2, o, _), S),
    member((1, o, [h]), S).

% Hydroxyl group
hydroxyl((_, o, [h])).

% Primary amine
amine((_, n, [h, h])).
```

### 6.4 Ring Predicates

```prolog
% Get ring members from any context
get_ring(mol(ring(M)), M).
get_ring((_, ring(M)), M).
get_ring(mol(_, Subs), M) :- member(S, Subs), get_ring(S, M).
get_ring((_, _, Subs), M) :- member(S, Subs), get_ring(S, M).

% Ring contains a specific element
ring_contains(ring(Members), E) :-
    member((_, E, _), Members).

% Aromatic ring (all bonds are 'a')
aromatic_ring(ring(Members)) :-
    forall(member((B, _, _), Members), B = a).

% Ring size
ring_size(ring(Members), N) :- length(Members, N).

% Heterocyclic ring (contains a non-carbon atom)
heterocyclic(ring(Members)) :-
    member((_, E, _), Members), E \= c.
```

### 6.5 Advantages of This Representation for MIL

1. **Uniform traversal**: `member/2` works on substituent lists regardless of atom valence.

2. **Ring-as-element**: Rings occupy the element position, so pattern matching naturally separates ring vs non-ring cases. A predicate can match `(_, ring(_))` to detect any ring.

3. **Self-contained rings**: All ring atoms are inside the ring term, so ring-level predicates (`aromatic_ring`, `ring_contains`, `ring_size`) operate on the ring term alone without needing context from a parent.

4. **Cyclic bond description**: The first member's bond type closing the ring means `ring(Members)` fully describes its own bonding — no external closing bond argument needed.

5. **Predicate invention**: MIL can compose `member`, `element`, `ring_contains`, and functional group predicates to discover patterns like "molecules with a hydroxyl on a nitrogen-containing aromatic ring."

6. **mol vs sub-atom distinction**: Lets the system distinguish complete molecules from fragments when relevant, while sharing traversal logic through `element/2` and `subs/2`.

---

## 7. Parsing: SMILES to Prolog Terms

### 7.1 Overview

The parser converts a SMILES string into the nested term representation described above. The process has three stages:

1. **Tokenisation**: Convert the SMILES string into a sequence of tokens (atoms, bonds, branch-open, branch-close, ring-open/close digits).
2. **Graph Construction**: Build an explicit molecular graph (atoms + bonds + ring membership) by interpreting the token sequence.
3. **Term Generation**: Convert the molecular graph into the canonical nested `mol`/tuple/ring representation.

### 7.2 Stage 1: Tokenisation

SMILES tokens are:

| Token Type | Pattern | Examples |
|-----------|---------|----------|
| Organic atom | `[BCNOPSFI]` or `Cl`, `Br` | `C`, `N`, `Cl` |
| Aromatic atom | `[bcnops]` | `c`, `n` |
| Bracketed atom | `[...]` | `[NH4+]`, `[Fe]` |
| Bond | `-`, `=`, `#`, `:` | `=`, `#` |
| Branch open | `(` | |
| Branch close | `)` | |
| Ring digit | `0-9` or `%NN` | `1`, `%12` |
| Dot | `.` | (disconnected fragments) |

The tokeniser produces a flat list of typed tokens. Multi-character atoms (`Cl`, `Br`) and bracket expressions must be handled as single tokens.

### 7.3 Stage 2: Graph Construction

Process the token stream with a stack-based algorithm:

1. Maintain a **current atom** pointer and a **branch stack**.
2. For each atom token: create a new atom node, connect it to the current atom with the pending bond type (default: single, or aromatic if both atoms are aromatic).
3. For `(`: push the current atom onto the branch stack.
4. For `)`: pop the branch stack and set current atom to the popped value. Clear any pending bond.
5. For bond tokens (`=`, `#`, `:`, `/`, `\`): set the pending bond type for the next atom. Stereo bonds (`/`, `\`) are treated as single bonds in the graph.
6. For ring digits:
   - If this digit is **unseen**: record the current atom and pending bond as an open ring.
   - If this digit was **previously opened**: create a bond between the current atom and the stored atom, closing the ring. Both sides may carry an explicit bond type; if both are specified, they must agree.
7. For `.` (dot): reset the current atom to None, starting a disconnected fragment.
8. After processing all tokens: verify no unclosed rings or branches remain.

**Hydrogen filling**: After graph construction, compute implicit hydrogens for each atom based on its element type and current bond order sum, following standard SMILES valence rules:

- **Organic subset atoms** (no brackets): find the lowest valid valence ≥ bond order sum, then `implicit_h = valence - bond_order_sum`. For aromatic atoms, subtract 1 further (minimum 0) since one electron participates in the π system.
- **Bracket atoms**: if an explicit H count is given (e.g., `[NH4+]`), use it directly. If no H field is present (e.g., `[Fe]`, `[O-]`), implicit hydrogen count is **0** per the OpenSMILES specification.

**Aromatic bond inference**: When two adjacent atoms are both aromatic (lowercase in SMILES) and no explicit bond symbol is given, the bond type is `a` (aromatic). This applies to both chain bonds and ring-closure bonds. An explicit `-` between two aromatic atoms forces a single bond (e.g., biphenyl `c1ccccc1-c2ccccc2`). When one atom is aromatic and the other is not, the default bond is single.

**Default valences** used for implicit hydrogen computation:

| Element | Valid valences |
|---------|---------------|
| B | 3 |
| C | 4 |
| N | 3, 5 |
| O | 2 |
| P | 3, 5 |
| S | 2, 4, 6 |
| F, Cl, Br, I | 1 |

### 7.4 Stage 3: Term Generation

Convert the molecular graph into the nested term representation:

1. **Identify ring systems**: Group ring bonds into ring systems. For each ring, determine the ordered set of member atoms.
2. **Select root**: If rings exist, the primary ring system is the root. Otherwise, select the highest-branching atom.
3. **Generate ring terms**: For each ring system, produce `ring(Members)` with the first member's bond type encoding the closing bond. Each member's substituent list includes non-ring children.
4. **Wrap root**: If the root is a ring, produce `mol(ring(Members))`. If an atom, produce `mol(Element, Subs)`.
5. **Tree traversal**: DFS from the root. Non-ring children become `(BondType, Element, Subs)` tuples. Ring children become `(BondType, ring(Members))` tuples.
6. **Canonicalise** member ordering and substituent list ordering per Section 5.
7. **Fill hydrogens** as `h` entries in substituent lists.

### 7.5 Implementation

The parser is implemented in **Rust** as a standalone crate (`smiles-parser`) for eventual integration with the Prolog2 system. Zero external dependencies.

**Stage 1 (Tokenisation)** — implemented in `src/token.rs`:

- Hand-rolled character-by-character scanner with one-character lookahead.
- Handles two-character atoms (`Cl`, `Br`), two-digit ring closures (`%10`–`%99`), bracket atoms with isotope/chirality/H-count/charge, aromatic bracket atoms (`[se]`, `[as]`), wildcard (`[*]`), and old-style charges (`++`, `---`).
- 12 unit tests covering all token types and error reporting.

**Stage 2 (Graph Construction)** — implemented in `src/graph.rs`:

- Adjacency-list molecular graph (`MolGraph`) with typed atoms and bonds.
- Stack-based algorithm processing tokens left-to-right with branch stack and ring-open map.
- Aromatic bond inference: implicit aromatic bonds between two lowercase atoms; explicit `-` overrides to single.
- Ring closure bonds: both open and close sides may carry explicit bond types, which must agree.
- Implicit hydrogen computation using standard SMILES valence rules with multi-valence support (N: 3/5, S: 2/4/6, P: 3/5).
- Bracket atom H convention: no H field → 0 implicit hydrogens (per OpenSMILES).
- Disconnected fragments via `.` (dot) operator.
- Error detection: unclosed rings, unmatched branches, conflicting ring bond types.
- 52 unit tests covering atom counts, element identification, bond types, branching, ring closures, implicit hydrogens (all element types), bracket atoms (charges, isotopes, explicit H), dot fragments, connectivity, larger molecules (toluene, benzoic acid, furan, pyrrole, biphenyl), and error cases.

**Stage 3 (Term Generation)** — not yet implemented. Will convert the molecular graph into nested Prolog terms per the representation design in Sections 2–3.

### 7.6 SELFIES as Alternative Input

SELFIES (Self-Referencing Embedded Strings) is a more recent molecular notation designed for machine learning applications. Every valid SELFIES string encodes a valid molecule, unlike SMILES where syntax errors are possible. The parsing pipeline for SELFIES would follow the same three-stage structure, with Stage 1 adapted for SELFIES token syntax. Stages 2 and 3 remain identical.

---

## 8. Open Questions and Future Work

1. **Fused ring encoding**: The `edge` functor design needs validation against complex polycyclic systems (steroids, porphyrins). Bridged bicyclic systems (e.g., norbornane) may require additional constructs.

2. **Stereochemistry**: SMILES can encode E/Z isomerism (`/`, `\`) and chirality (`@`, `@@`). These are currently not represented. They could be added as wrapper functors or additional list elements where relevant.

3. **Charges and radicals**: Charged atoms (`[NH4+]`, `[O-]`) and radical species need representation. A `charge/1` wrapper in the substituent list is one option: `mol(n, [h, h, h, h, charge(1)])`.

4. **Leaf atom shorthand**: Whether `(2, o)` should be permitted as sugar for `(2, o, [])`. The canonical form always includes the list, but parser output could accept either.

5. **Ring starting position**: The convention for which member starts the ring list needs thorough specification, particularly for symmetric rings and rings with multiple substituents.

6. **Canonicalisation testing**: The canonical ordering rules need thorough testing against molecular databases to ensure uniqueness and consistency.

7. **Performance at scale**: For SAR datasets with thousands of molecules, the term representation should be benchmarked against the traditional atom-bond fact approach for both memory usage and learning speed.

8. **Integration with Prolog2**: Meta-rules for molecular pattern learning would operate over the `mol`/tuple/ring structure, potentially inventing predicates that recognise novel functional groups or ring-substituent patterns.