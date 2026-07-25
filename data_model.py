"""
molmol/datamodel.py
-------------------
Strutture dati centrali del nuovo geometry/contact core.
Nessuna dipendenza da Blender o RDKit: testabile in Python puro.

Convenzione vertex group:
    contact_<lato>_<indice>
    es: contact_a_0, contact_b_0, contact_b_1
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Livello 1 - superfici di contatto
# ---------------------------------------------------------------------------

@dataclass
class contact_patch:
    """
    Superficie di contatto reale derivata da un vertex group Blender.
    Tutti i vettori sono in world space al momento dell'estrazione.
    """
    vertex_group_name: str      # es. 'contact_a_0'
    centroid: np.ndarray        # shape (3,)  - baricentro dei vertici del vg
    normal: np.ndarray          # shape (3,)  - normale media outward, normalizzata
    sample_points: np.ndarray   # shape (N,3) - vertici world del vg

    def __post_init__(self) -> None:
        self.centroid      = np.asarray(self.centroid,      dtype=float)
        self.normal        = np.asarray(self.normal,        dtype=float)
        self.sample_points = np.asarray(self.sample_points, dtype=float)
        n = np.linalg.norm(self.normal)
        if n > 1e-8:
            self.normal = self.normal / n

    @property
    def side(self) -> str:
        """
        Lato dal nome del vertex group.
        'contact_a_0' -> 'a'
        Ritorna '' se il nome non segue la convenzione.
        """
        parts = self.vertex_group_name.split("_")
        return parts[1] if len(parts) >= 3 else ""

    @property
    def index(self) -> int:
        """
        Indice numerico dal nome del vertex group.
        'contact_a_0' -> 0
        Ritorna -1 se non parsabile.
        """
        parts = self.vertex_group_name.split("_")
        try:
            return int(parts[2]) if len(parts) >= 3 else -1
        except ValueError:
            return -1

    def __repr__(self) -> str:
        return (
            f"contact_patch(vg='{self.vertex_group_name}', "
            f"centroid={self.centroid.round(4).tolist()}, "
            f"normal={self.normal.round(4).tolist()}, "
            f"n_points={len(self.sample_points)})"
        )


# ---------------------------------------------------------------------------
# Livello 2 - frammenti piazzati in scena
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Livello 3 - legami tra frammenti
# ---------------------------------------------------------------------------

@dataclass
class bond_edge:
    """
    Collegamento tra due placed_fragment tramite le rispettive
    superfici di contatto.

    transform_b_to_a: matrice 4x4 numpy che porta l'oggetto B
    in posizione allineata rispetto ad A; calcolata da contact_aligner.
    """
    node_a: str                                 # node_id del frammento A
    node_b: str                                 # node_id del frammento B
    patch_a: str                                # vg_name su A
    patch_b: str                                # vg_name su B
    bond_order: int = 1                         # 1=singolo, 2=doppio, 3=triplo
    transform_b_to_a: Optional[np.ndarray] = None  # 4x4, popolato da contact_aligner
    resolved: bool = False                      # True dopo che l'allineamento e applicato

    def __repr__(self) -> str:
        return (
            f"bond_edge({self.node_a}[{self.patch_a}] "
            f"<-{self.bond_order}-> "
            f"{self.node_b}[{self.patch_b}], "
            f"resolved={self.resolved})"
        )


# ---------------------------------------------------------------------------
# Livello 4 - grafo di assemblaggio arricchito
# ---------------------------------------------------------------------------

@dataclass
class contact_assembly_graph:
    """
    Grafo completo: nodi = placed_fragment, archi = bond_edge.

    Affianca AssemblyGraph di assembly_planner aggiungendo
    la dimensione geometrica reale (contact_patch).
    """
    fragments: Dict[str, placed_fragment] = field(default_factory=dict)
    edges: List[bond_edge] = field(default_factory=list)

    def add_fragment(self, frag: placed_fragment) -> None:
        self.fragments[frag.node_id] = frag

    def add_edge(self, edge: bond_edge) -> None:
        self.edges.append(edge)

    def get_fragment(self, node_id: str) -> Optional[placed_fragment]:
        return self.fragments.get(node_id)

    def edges_for(self, node_id: str) -> List[bond_edge]:
        """Tutti gli archi che coinvolgono node_id."""
        return [e for e in self.edges if node_id in (e.node_a, e.node_b)]

    def unresolved_edges(self) -> List[bond_edge]:
        return [e for e in self.edges if not e.resolved]

    def __repr__(self) -> str:
        return (
            f"contact_assembly_graph("
            f"{len(self.fragments)} fragments, "
            f"{len(self.edges)} edges, "
            f"{len(self.unresolved_edges())} unresolved)"
        )
