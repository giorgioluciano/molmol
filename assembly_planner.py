from __future__ import annotations

from dataclasses import dataclass
from rdkit import Chem
from .fragment_recognizer import FragmentMatch, find_fragment_matches
from .fragment_library import get_fragment_template


@dataclass
class FragmentInstance:
    id: str
    template_id: str
    source_atom_indices: list[int]


@dataclass
class LinkerAtom:
    id: str
    element: str
    source_atom_index: int


@dataclass
class AssemblyNode:
    id: str
    kind: str
    fragment_instance: FragmentInstance | None = None
    linker_atom: LinkerAtom | None = None


@dataclass
class AssemblyEdge:
    id: str
    node_a_id: str
    node_b_id: str
    bond_order: int = 1


@dataclass
class AssemblyGraph:
    nodes: dict[str, AssemblyNode]
    edges: list[AssemblyEdge]
    source_name: str | None = None


def select_non_overlapping_matches(matches: list[FragmentMatch]) -> list[FragmentMatch]:
    selected = []
    occupied = set()
    for match in matches:
        atoms = set(match.atom_indices)
        if atoms & occupied:
            continue
        selected.append(match)
        occupied.update(atoms)
    return selected


def reorder_match_atom_indices_by_template(mol: Chem.Mol, match: FragmentMatch) -> tuple[int, ...]:
    atom_indices = tuple(int(i) for i in match.atom_indices)
    template = get_fragment_template(match.fragment_id)

    if template.smarts is not None:
        pattern = Chem.MolFromSmarts(template.smarts)
        if pattern is not None:
            for rdkit_match in mol.GetSubstructMatches(pattern):
                if set(rdkit_match) == set(atom_indices):
                    return tuple(rdkit_match)

    # fallback: ordine originale dal recognizer
    return atom_indices


def build_assembly_graph(mol: Chem.Mol) -> AssemblyGraph:
    matches = select_non_overlapping_matches(find_fragment_matches(mol))
    nodes: dict[str, AssemblyNode] = {}
    atom_to_node: dict[int, str] = {}
    used_atoms: set[int] = set()

    for i, match in enumerate(matches):
        node_id = f"frag_{i:03d}"
        normalized = reorder_match_atom_indices_by_template(mol, match)
        fi = FragmentInstance(id=node_id, template_id=match.fragment_id, source_atom_indices=list(normalized))
        nodes[node_id] = AssemblyNode(id=node_id, kind="fragment", fragment_instance=fi)
        for atom_idx in match.atom_indices:
            atom_to_node[int(atom_idx)] = node_id
            used_atoms.add(int(atom_idx))

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        if atom_idx in used_atoms:
            continue
        node_id = f"atom_{atom_idx:03d}"
        nodes[node_id] = AssemblyNode(
            id=node_id,
            kind="linker_atom",
            linker_atom=LinkerAtom(id=node_id, element=atom.GetSymbol(), source_atom_index=atom_idx),
        )
        atom_to_node[atom_idx] = node_id

    edges = []
    seen_pairs: set[frozenset[str]] = set()
    edge_count = 0

    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        node_a = atom_to_node[a]
        node_b = atom_to_node[b]
        if node_a == node_b:
            continue

        pair = frozenset((node_a, node_b))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)

        order = max(1, min(3, int(round(bond.GetBondTypeAsDouble()))))
        edges.append(
            AssemblyEdge(
                id=f"edge_{edge_count:03d}",
                node_a_id=node_a,
                node_b_id=node_b,
                bond_order=order,
            )
        )
        edge_count += 1

    return AssemblyGraph(
        nodes=nodes,
        edges=edges,
        source_name=mol.GetProp("_Name") if mol.HasProp("_Name") else None,
    )
