from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np



@dataclass
class FragmentTemplate:
    id: str
    name: str
    smarts: str
    atom_names: List[str]
    atom_elements: List[str]
    atom_local_positions: np.ndarray
    anchor_atom_indices: List[int] = field(default_factory=list)
    symmetry_label: Optional[str] = None
    blender_template_name: Optional[str] = None   # root object name in .blend
    blender_collection_name: Optional[str] = None # collection name in .blend
    tags: List[str] = field(default_factory=list)
    molmol_recipe_id: Optional[str] = None
    primary_attachment_atom_local_index: Optional[int] = None
    atom_marker_names: List[str] = field(default_factory=list)
    internal_bond_pairs: List[Tuple[int, int]] = field(default_factory=list)



def _v(x: float, y: float, z: float) -> np.ndarray:
    return np.array([x, y, z], dtype=float)


# ---------------------------------------------------------------------------
# Fragment templates — names aligned to current Blender dump
# ---------------------------------------------------------------------------

CARBONYL = FragmentTemplate(
    id="carbonyl",
    name="Carbonyl",
    smarts="[C:1]=[O:2]",
    atom_names=["C", "O"],
    atom_elements=["C", "O"],
    atom_local_positions=np.array([
        [0.00, 0.00, 0.00],
        [1.23, 0.00, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1],
    symmetry_label="planar_sp2",
    blender_template_name="c_o_root",
    blender_collection_name="frag_c_o",
    atom_marker_names=[
        "c_o_c_mk_atom_1",
        "c_o_o_mk_atom_1",
    ],
    internal_bond_pairs=[
        (0, 1),
    ],
    tags=["planar", "sp2", "functional_group"],
)



ALKENE = FragmentTemplate(
    id="alkene",
    name="Alkene",
    smarts="[C:1]=[C:2]",
    atom_names=["C1", "C2"],
    atom_elements=["C", "C"],
    atom_local_positions=np.array([
        [-0.67, 0.00, 0.00],
        [0.67, 0.00, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1],
    symmetry_label="double_bond_planar",
    blender_template_name="c_c_double_root",
    blender_collection_name="frag_c_c_double",
    atom_marker_names=[
        "c_c_double_mesh_c1_mk_atom",
        "c_c_double_mesh_c2_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
    ],
    tags=["planar", "double_bond"],
)



ALKYNE = FragmentTemplate(
    id="alkyne",
    name="Alkyne",
    smarts="[C:1]#[C:2]",
    atom_names=["C1", "C2"],
    atom_elements=["C", "C"],
    atom_local_positions=np.array([
        [-0.60, 0.00, 0.00],
        [0.60, 0.00, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1],
    symmetry_label="linear",
    blender_template_name="c_c_triple_root",
    blender_collection_name="frag_c_c_triple",
    atom_marker_names=[
        "c_c_triple_c1_mk_atom",
        "c_c_triple_c2_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
    ],
    tags=["linear", "triple_bond"],
)



NITRILE = FragmentTemplate(
    id="nitrile",
    name="Nitrile",
    smarts="[C:1]#[N:2]",
    atom_names=["C", "N"],
    atom_elements=["C", "N"],
    atom_local_positions=np.array([
        [0.00, 0.00, 0.00],
        [1.16, 0.00, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1],
    symmetry_label="linear_terminal",
    blender_template_name="c_n_triple_root",
    blender_collection_name="frag_c_n_triple",
    atom_marker_names=[
        "c_n_triple_c_mk_atom",
        "c_n_triple_n_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
    ],
    tags=["linear", "terminal", "triple_bond_like"],
    molmol_recipe_id="nitrile_linear_v1",
)



BENZENE = FragmentTemplate(
    id="benzene",
    name="Benzene",
    smarts="c1ccccc1",
    atom_names=["C1", "C2", "C3", "C4", "C5", "C6"],
    atom_elements=["C", "C", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [1.40, 0.00, 0.00],
        [0.70, 1.21, 0.00],
        [-0.70, 1.21, 0.00],
        [-1.40, 0.00, 0.00],
        [-0.70, -1.21, 0.00],
        [0.70, -1.21, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4, 5],
    symmetry_label="ring_c6",
    blender_template_name="benzene_root",
    blender_collection_name="frag_benzene",
    atom_marker_names=[
        "benzene_c1_mk_atom",
        "benzene_c2_mk_atom",
        "benzene_c3_mk_atom",
        "benzene_c4_mk_atom",
        "benzene_c5_mk_atom",
        "benzene_c6_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
    ],
    tags=["aromatic", "planar", "ring"],
)

PHENYL = FragmentTemplate(
    id="phenyl",
    name="Phenyl",
    smarts="c1ccccc1",
    atom_names=["C1", "C2", "C3", "C4", "C5", "C6"],
    atom_elements=["C", "C", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [1.40, 0.00, 0.00],
        [0.70, 1.21, 0.00],
        [-0.70, 1.21, 0.00],
        [-1.40, 0.00, 0.00],
        [-0.70, -1.21, 0.00],
        [0.70, -1.21, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4, 5],
    symmetry_label="ring_c6_substitutable",
    blender_template_name="phenyl_root",
    blender_collection_name="frag_phenyl",
    atom_marker_names=[
        "phenyl_mesh_c1_mk_atom",
        "phenyl_mesh_c2_mk_atom",
        "phenyl_mesh_c3_mk_atom",
        "phenyl_mesh_c4_mk_atom",
        "phenyl_mesh_c5_mk_atom",
        "phenyl_mesh_c6_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
    ],
    tags=["aromatic", "planar", "ring", "substitutable"],
)


COO = FragmentTemplate(
    id="coo",
    name="COO Group",
    smarts="[C:1](=[O:2])[O:3]",
    atom_names=["C", "O_double", "O_single"],
    atom_elements=["C", "O", "O"],
    atom_local_positions=np.array([
        [0.00, 0.00, 0.00],
        [1.23, 0.00, 0.00],
        [-0.67, 1.16, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2],
    symmetry_label="planar_sp2_coo",
    blender_template_name="c_o_o_root",
    blender_collection_name="frag_c_o_o",
    atom_marker_names=[
        "c_o_o_c_mk_atom",
        "c_o_o_o1_mk_atom",
        "c_o_o_o2_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (0, 2),
    ],
    tags=["planar", "sp2", "functional_group", "ester", "carboxylic_acid"],
)


CH3 = FragmentTemplate(
    id="ch3",
    name="Methyl",
    smarts="[CH3:1]",
    atom_names=["C"],
    atom_elements=["C"],
    atom_local_positions=np.array([
        [0.00, 0.00, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0],
    symmetry_label="tetrahedral_methyl",
    blender_template_name="ch3_root",
    blender_collection_name="frag_ch3",
    atom_marker_names=[],
    internal_bond_pairs=[],
    tags=["methyl", "terminal"],
)


NAPHTHALENE = FragmentTemplate(
    id="naphthalene",
    name="Naphthalene",
    smarts="c1ccc2ccccc2c1",
    atom_names=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"],
    atom_elements=["C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [2.10, 0.00, 0.00],
        [1.40, 1.21, 0.00],
        [0.00, 1.21, 0.00],
        [-0.70, 0.00, 0.00],
        [0.00, -1.21, 0.00],
        [1.40, -1.21, 0.00],
        [-2.10, 0.00, 0.00],
        [-1.40, -1.21, 0.00],
        [-2.80, -1.21, 0.00],
        [-2.80, 1.21, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    symmetry_label="fused_ring_c10",
    blender_template_name="naphthalene_root",
    blender_collection_name="frag_naphthalene",
    atom_marker_names=[
        "naphthalene_c1_mk_atom",
        "naphthalene_c2_mk_atom",
        "naphthalene_c3_mk_atom",
        "naphthalene_c4_mk_atom",
        "naphthalene_c5_mk_atom",
        "naphthalene_c6_mk_atom",
        "naphthalene_c7_mk_atom",
        "naphthalene_c8_mk_atom",
        "naphthalene_c9_mk_atom",
        "naphthalene_c10_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (3, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 2),
    ],
    tags=["aromatic", "planar", "ring", "fused_ring"],
)


PYRIDINE = FragmentTemplate(
    id="pyridine",
    name="Pyridine",
    smarts="n1ccccc1",
    atom_names=["N1", "C1", "C2", "C3", "C4", "C5"],
    atom_elements=["N", "C", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [1.40, 0.00, 0.00],
        [0.70, 1.21, 0.00],
        [-0.70, 1.21, 0.00],
        [-1.40, 0.00, 0.00],
        [-0.70, -1.21, 0.00],
        [0.70, -1.21, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4, 5],
    symmetry_label="ring_c6_hetero",
    blender_template_name="pyridine_root",
    blender_collection_name="frag_pyridine",
    atom_marker_names=[
        "pyridine_n1_mk_atom",
        "pyridine_c1_mk_atom",
        "pyridine_c2_mk_atom",
        "pyridine_c3_mk_atom",
        "pyridine_c4_mk_atom",
        "pyridine_c5_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
    ],
    tags=["aromatic", "planar", "ring", "heterocycle"],
)

QUINOLINE = FragmentTemplate(
    id="quinoline",
    name="Quinoline",
    smarts="c1ccc2ncccc2c1",
    atom_names=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "N1"],
    atom_elements=["C", "C", "C", "C", "C", "C", "C", "C", "C", "N"],
    atom_local_positions=np.array([
        [2.10, 0.00, 0.00],
        [1.40, 1.21, 0.00],
        [0.00, 1.21, 0.00],
        [-0.70, 0.00, 0.00],
        [0.00, -1.21, 0.00],
        [1.40, -1.21, 0.00],
        [-2.10, 0.00, 0.00],
        [-1.40, -1.21, 0.00],
        [-2.80, -1.21, 0.00],
        [-1.40, 1.21, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    symmetry_label="fused_ring_quinoline",
    blender_template_name="quinoline_root",
    blender_collection_name="frag_quinoline",
    atom_marker_names=[
        "quinoline_c1_mk_atom",
        "quinoline_c2_mk_atom",
        "quinoline_c3_mk_atom",
        "quinoline_c4_mk_atom",
        "quinoline_c5_mk_atom",
        "quinoline_c6_mk_atom",
        "quinoline_c7_mk_atom",
        "quinoline_c8_mk_atom",
        "quinoline_c9_mk_atom",
        "quinoline_n1_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (3, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 2),
    ],
    tags=["aromatic", "planar", "ring", "fused_ring", "heterocycle"],
)


THIOPHENE = FragmentTemplate(
    id="thiophene",
    name="Thiophene",
    smarts="s1cccc1",
    atom_names=["S1", "C1", "C2", "C3", "C4"],
    atom_elements=["S", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [0.00, 1.25, 0.00],
        [1.19, 0.39, 0.00],
        [0.74, -1.01, 0.00],
        [-0.74, -1.01, 0.00],
        [-1.19, 0.39, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4],
    symmetry_label="ring_c5_hetero",
    blender_template_name="thiophene_root",
    blender_collection_name="frag_thiophene",
    atom_marker_names=[
        "thiophene_s_mk_atom_1",
        "thiophene_c1_mk_atom_1",
        "thiophene_c2_mk_atom_1",
        "thiophene_c3_mk_atom_1",
        "thiophene_c4_mk_atom_1",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
    ],
    tags=["aromatic", "planar", "ring", "heterocycle"],
)


FURAN = FragmentTemplate(
    id="furan",
    name="Furan",
    smarts="o1cccc1",
    atom_names=["O1", "C1", "C2", "C3", "C4"],
    atom_elements=["O", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [0.00, 1.25, 0.00],
        [1.19, 0.39, 0.00],
        [0.74, -1.01, 0.00],
        [-0.74, -1.01, 0.00],
        [-1.19, 0.39, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4],
    symmetry_label="ring_c5_hetero",
    blender_template_name="furane_root",
    blender_collection_name="frag_furane",
    atom_marker_names=[
        "furane_o1_mk_atom",
        "furane_c1_mk_atom",
        "furane_c2_mk_atom",
        "furane_c3_mk_atom",
        "furane_c4_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
    ],
    tags=["aromatic", "planar", "ring", "heterocycle"],
)

IMIDAZOLE = FragmentTemplate(
    id="imidazole",
    name="Imidazole",
    smarts="n1c[nH]cc1",
    atom_names=["N1", "C1", "C2", "N2", "C4"],
    atom_elements=["N", "C", "C", "N", "C"],
    atom_local_positions=np.array([
        [0.00, 1.25, 0.00],
        [1.19, 0.39, 0.00],
        [0.74, -1.01, 0.00],
        [-0.74, -1.01, 0.00],
        [-1.19, 0.39, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4],
    symmetry_label="ring_c5_di_hetero",
    blender_template_name="imidazole_root",
    blender_collection_name="frag_imidazole",
    atom_marker_names=[
        "imidazole_n1_mk_atom",
        "imidazole_c1_mk_atom",
        "imidazole_c2_mk_atom",
        "imidazole_n2_mk_atom",
        "imidazole_c4_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
    ],
    tags=["aromatic", "planar", "ring", "heterocycle", "di_hetero"],
)


PYRROLE = FragmentTemplate(
    id="pyrrole",
    name="Pyrrole",
    smarts="[nH]1cccc1",
    atom_names=["N1", "C1", "C2", "C3", "C4"],
    atom_elements=["N", "C", "C", "C", "C"],
    atom_local_positions=np.array([
        [0.00, 1.25, 0.00],
        [1.19, 0.39, 0.00],
        [0.74, -1.01, 0.00],
        [-0.74, -1.01, 0.00],
        [-1.19, 0.39, 0.00],
    ], dtype=float),
    anchor_atom_indices=[0, 1, 2, 3, 4],
    symmetry_label="ring_c5_hetero",
    blender_template_name="pyrrole_root",
    blender_collection_name="frag_pyrrole",
    atom_marker_names=[
        "pyrrole_n_mk_atom",
        "pyrrole_c1_mk_atom",
        "pyrrole_c2_mk_atom",
        "pyrrole_c3_mk_atom",
        "pyrrole_c4_mk_atom",
    ],
    internal_bond_pairs=[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
    ],
    tags=["aromatic", "planar", "ring", "heterocycle"],
)


FRAGMENT_LIBRARY: Dict[str, FragmentTemplate] = {
    CARBONYL.id: CARBONYL,
    ALKENE.id: ALKENE,
    ALKYNE.id: ALKYNE,
    NITRILE.id: NITRILE,
    BENZENE.id: BENZENE,
    PHENYL.id: PHENYL,
    COO.id: COO,
    CH3.id: CH3,
    NAPHTHALENE.id: NAPHTHALENE,
    PYRIDINE.id: PYRIDINE,
    QUINOLINE.id: QUINOLINE,
    THIOPHENE.id: THIOPHENE,
    FURAN.id: FURAN,
    IMIDAZOLE.id: IMIDAZOLE,
    PYRROLE.id: PYRROLE,
}


def get_fragment_template(fragment_id: str) -> FragmentTemplate:
    return FRAGMENT_LIBRARY[fragment_id]


def list_fragment_templates() -> List[FragmentTemplate]:
    return list(FRAGMENT_LIBRARY.values())


FragmentTemplate = FragmentTemplate
fragment_template = FragmentTemplate

# ---------------------------------------------------------------------------
# Atom template names — names aligned to current Blender dump
# ---------------------------------------------------------------------------

ATOM_TEMPLATE_NAMES: Dict[str, str] = {
    "tpl_tetra_sp3": "tetra_root",
    "tpl_single_17": "single_17_root",
    "tpl_single_23": "single_23_root",
    "tpl_linear": "linear_root",
    "tpl_trigonal_planar": "trigonal_planar_root",
    "tpl_trigonal_pyramidal": "trigonal_pyramidal_root",
    "tpl_trigonal_pyramidal_lp": "trigonal_pyramidal_with_lp_root",
    "tpl_trigonal_bipyramidal": "trigonal_bipyramidal_root",
    "tpl_octahedral": "octahedral_23_root",
    "tpl_bent_lp": "tpl_bent_with_lp_root",
}


ATOM_TEMPLATE_COLLECTIONS: Dict[str, str] = {
    "tpl_tetra_sp3": "tpl_tetra",
    "tpl_single_17": "tpl_single_17",
    "tpl_single_23": "tpl_single_23",
    "tpl_linear": "tpl_linear",
    "tpl_trigonal_planar": "tpl_trigonal_planar",
    "tpl_trigonal_pyramidal": "tpl_trigonal_pyramidal",
    "tpl_trigonal_pyramidal_lp": "tpl_trigonal_pyramidal_with_lp",
    "tpl_trigonal_bipyramidal": "tpl_trigonal_bipyramidal",
    "tpl_octahedral": "tpl_octahedral_23",
    "tpl_bent_lp": "tpl_bent_with_lp",
}


def get_atom_template_name(element: str, residual_contact_count: int) -> str:
    """
    Returns the Blender object name for a linker atom template,
    based on element and number of remaining connections.
    """
    element = (element or "").strip().upper()

    if element == "H":
        return ATOM_TEMPLATE_NAMES["tpl_single_17"]
    if element == "C":
        return ATOM_TEMPLATE_NAMES["tpl_tetra_sp3"]
    if residual_contact_count == 1:
        return ATOM_TEMPLATE_NAMES["tpl_single_23"]
    if residual_contact_count == 2:
        return ATOM_TEMPLATE_NAMES["tpl_linear"]
    if residual_contact_count == 3:
        return ATOM_TEMPLATE_NAMES["tpl_trigonal_planar"]
    return ATOM_TEMPLATE_NAMES["tpl_tetra_sp3"]


def get_atom_template_collection(template_key: str) -> Optional[str]:
    """Returns the Blender collection name for an atom template key."""
    return ATOM_TEMPLATE_COLLECTIONS.get(template_key)


# ---------------------------------------------------------------------------
# Bond template registry (canonical naming)
# ---------------------------------------------------------------------------

BOND_TEMPLATE_NAMES: Dict[str, str] = {
    "single": "bond_single_medium_root",
    "single_short": "bond_single_short_root",
    "single_half": "half_bond_root",
    "opaque": "bond_opaque_root",
}

BOND_TEMPLATE_COLLECTIONS: Dict[str, str] = {
    "single": "tpl_bond_single_medium",
    "single_short": "tpl_bond_single_short",
    "single_half": "tpl_half_bond",
    "opaque": "tpl_bond_opaque",
}

BOND_TEMPLATE_LEGACY_NAMES: Dict[str, List[str]] = {
    "single": [
        "bond_single_medium",
        "bond_single",
    ],
    "single_short": [
        "bond_single_short",
        "bond_single_short_root",
    ],
    "single_half": [
        "half_bond",
        "half_bond_root",
    ],
    "opaque": [
        "bond_opaque",
        "bond_opaque_root",
    ],
}


def get_bond_template_key(bond_order: int) -> str:
    """
    Return the canonical logical key for an independent bond template.

    Only single independent bonds are supported by the current pipeline.
    Double and triple bonds are embedded in fragment templates.
    """
    bond_order = int(bond_order)
    if bond_order != 1:
        raise ValueError(
            f"Independent bond templates are only supported for bond_order=1, got {bond_order}"
        )
    return "single"


def get_bond_template_name(bond_order: int) -> str:
    """
    Return the canonical Blender root object name for an independent bond template.
    """
    template_key = get_bond_template_key(bond_order)
    return BOND_TEMPLATE_NAMES[template_key]


def get_bond_template_collection_name(bond_order: int) -> str:
    """
    Return the Blender collection name containing the requested bond template.
    """
    template_key = get_bond_template_key(bond_order)
    return BOND_TEMPLATE_COLLECTIONS[template_key]


def get_bond_template_legacy_names(bond_order: int) -> List[str]:
    """
    Return legacy-compatible root object names that may still exist in old libraries.
    """
    template_key = get_bond_template_key(bond_order)
    return list(BOND_TEMPLATE_LEGACY_NAMES.get(template_key, []))