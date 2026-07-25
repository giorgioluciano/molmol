from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from rdkit import Chem
from .fragment_library import list_fragment_templates


@dataclass(frozen=True)
class FragmentMatch:
    fragment_id: str
    fragment_name: str
    atom_indices: tuple[int, ...]
    smarts: str
    size: int


def load_molecule(path: str | Path):
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
        if mol is None:
            raise ValueError(f"Could not read MOL file: {path}")
        return mol

    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        for mol in supplier:
            if mol is not None:
                return mol
        raise ValueError(f"No valid molecules found in SDF file: {path}")

    raise ValueError(f"Unsupported file format: {suffix}")



def _reorder_by_expected_elements(mol, template, atom_indices: tuple[int, ...]) -> tuple[int, ...]:
    expected = list(template.atom_elements)
    actual = [mol.GetAtomWithIdx(i).GetSymbol() for i in atom_indices]
    if len(expected) != len(actual) or expected == actual:
        return tuple(int(i) for i in atom_indices)
    used = [False] * len(actual)
    reordered = []
    for symbol in expected:
        found = None
        for idx, real in enumerate(actual):
            if not used[idx] and real == symbol:
                found = idx
                break
        if found is None:
            return tuple(int(i) for i in atom_indices)
        reordered.append(int(atom_indices[found]))
        used[found] = True
    return tuple(reordered)

def find_fragment_matches(mol) -> list[FragmentMatch]:
    matches: list[FragmentMatch] = []

    for template in list_fragment_templates():
        query = Chem.MolFromSmarts(template.smarts)
        if query is None:
            raise ValueError(
                f"Invalid SMARTS in fragment library: {template.id} -> {template.smarts}"
            )

        for atom_indices in mol.GetSubstructMatches(query, uniquify=True):
            normalized = tuple(int(i) for i in atom_indices)
            matches.append(
                FragmentMatch(
                    fragment_id=template.id,
                    fragment_name=template.name,
                    atom_indices=normalized,
                    smarts=template.smarts,
                    size=len(normalized),
                )
            )

    matches.sort(key=lambda match: (-match.size, match.fragment_id, match.atom_indices))
    return matches



def find_fragment_matches_from_file(path: str | Path) -> list[FragmentMatch]:
    return find_fragment_matches(load_molecule(path))
