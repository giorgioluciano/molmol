from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Optional

import mathutils
import numpy as np
from rdkit import Chem

from .fragment_library import get_fragment_template


PIPELINE_EPS = 1e-8


# ---------------------------------------------------------------------------
# Runtime models
# ---------------------------------------------------------------------------

@dataclass
class FragmentInstance:
    id: str
    template_id: str
    source_atom_indices: list[int]


# ---------------------------------------------------------------------------
# RDKit coordinate utilities
# ---------------------------------------------------------------------------

def rdkit_atom_position(
    mol: Chem.Mol,
    atom_idx: int,
    center_offset: np.ndarray | None = None,
) -> np.ndarray:
    """
    Return RDKit atom 3D coordinates as a NumPy array, optionally centered.
    """
    conf = mol.GetConformer()
    pos = conf.GetAtomPosition(int(atom_idx))
    arr = np.array([pos.x, pos.y, pos.z], dtype=float)
    if center_offset is not None:
        arr -= center_offset
    return arr


def _ensure_conformer(mol: Chem.Mol):
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer")
    return mol.GetConformer()


def _bond_length(mol: Chem.Mol, a: int, b: int) -> float:
    pa = rdkit_atom_position(mol, int(a))
    pb = rdkit_atom_position(mol, int(b))
    return float(np.linalg.norm(pb - pa))


def compute_molecule_centroid(mol: Chem.Mol) -> np.ndarray:
    """
    Compute the centroid of all atoms in the molecule in RDKit coordinates.
    """
    conf = _ensure_conformer(mol)
    coords: list[list[float]] = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    return np.mean(coords, axis=0) if coords else np.array([0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# Rigid alignment / Kabsch
# ---------------------------------------------------------------------------

def _kabsch(
    source_points: np.ndarray,
    target_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Kabsch algorithm: find optimal rotation + translation that aligns source_points to target_points.
    Points are Nx3 arrays with one-to-one correspondence.
    """
    source_points = np.asarray(source_points, dtype=float)
    target_points = np.asarray(target_points, dtype=float)

    if (
        source_points.shape != target_points.shape
        or source_points.ndim != 2
        or source_points.shape[1] != 3
    ):
        raise ValueError("Expected matching Nx3 point arrays")

    source_centroid = source_points.mean(axis=0)
    target_centroid = target_points.mean(axis=0)

    source_centered = source_points - source_centroid
    target_centered = target_points - target_centroid

    covariance = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T

    # Reflection fix: enforce proper rotation (det(R) = +1)
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T

    translation = target_centroid - rotation @ source_centroid
    return rotation, translation


def rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """
    Root-mean-square deviation after optimal rigid alignment.
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape:
        raise ValueError("RMSD requires same shape")
    rotation, translation = _kabsch(P, Q)
    Pp = (rotation @ P.T).T + translation
    return float(np.sqrt(np.mean(np.sum((Pp - Q) ** 2, axis=1))))


def make_transform_matrix(rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """
    Build a 4x4 homogeneous transform matrix from rotation (3x3) and translation (3,).
    """
    transform = np.eye(4, dtype=float)
    transform[:3, :3] = np.asarray(rotation, dtype=float)
    transform[:3, 3] = np.asarray(translation, dtype=float).reshape(3)
    return transform


def best_fit_rigid(
    source_points: np.ndarray,
    target_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return optimal rotation + translation aligning source_points to target_points.
    """
    return _kabsch(source_points, target_points)


def best_fit_transform_matrix(
    source_points: np.ndarray,
    target_points: np.ndarray,
) -> np.ndarray:
    """
    Return a 4x4 transform matrix that best aligns source_points to target_points.
    """
    rotation, translation = _kabsch(source_points, target_points)
    return make_transform_matrix(rotation, translation)


def apply_transform_matrix(points: np.ndarray, transform: np.ndarray) -> np.ndarray:
    """
    Apply a 4x4 homogeneous transform matrix to a set of 3D points (Nx3).
    """
    points = np.asarray(points, dtype=float)
    transform = np.asarray(transform, dtype=float)
    ones = np.ones((points.shape[0], 1), dtype=float)
    hp = np.concatenate([points, ones], axis=1)
    out = (transform @ hp.T).T
    return out[:, :3]


def matrix4x4_to_blender(transform: np.ndarray) -> mathutils.Matrix:
    """
    Convert a NumPy 4x4 float array into a Blender mathutils.Matrix.
    """
    T = np.asarray(transform, dtype=float)
    return mathutils.Matrix(
        [
            [float(T[0, 0]), float(T[0, 1]), float(T[0, 2]), float(T[0, 3])],
            [float(T[1, 0]), float(T[1, 1]), float(T[1, 2]), float(T[1, 3])],
            [float(T[2, 0]), float(T[2, 1]), float(T[2, 2]), float(T[2, 3])],
            [float(T[3, 0]), float(T[3, 1]), float(T[3, 2]), float(T[3, 3])],
        ]
    )


# ---------------------------------------------------------------------------
# Anchor selection (if still used)
# ---------------------------------------------------------------------------

def choose_anchor_indices(
    local_positions: np.ndarray,
    anchor_indices: list[int],
) -> list[int]:
    """
    Choose up to 3 anchor atom indices that define a stable frame.

    If more than 2 anchor indices are available, pick the first triplet that is non-collinear.
    """
    if len(anchor_indices) < 2:
        raise ValueError("Need at least 2 anchor atom indices")
    if len(anchor_indices) == 2:
        return anchor_indices[:]

    pts = np.asarray(local_positions, dtype=float)
    first, second = anchor_indices[0], anchor_indices[1]
    p0, p1 = pts[first], pts[second]
    v01 = p1 - p0

    for idx in anchor_indices[2:]:
        v02 = pts[idx] - p0
        if np.linalg.norm(np.cross(v01, v02)) > 1e-6:
            return [first, second, idx]

    return anchor_indices[:2]


# ---------------------------------------------------------------------------
# Fragment placement: RDKit-driven fit onto Blender templates
# ---------------------------------------------------------------------------

def fit_fragment_pose(mol: Chem.Mol, fi: FragmentInstance, scale: float = 1.0) -> np.ndarray:
    """
    Compute a 4x4 transform matrix that aligns the Blender fragment template markers
    to the RDKit atom positions specified in FragmentInstance.source_atom_indices.

    RDKit controls the global pose; Blender markers provide the local template geometry.
    """
    import bpy

    def normalize_marker_name(name: str) -> str:
        if "." in name:
            base, suffix = name.rsplit(".", 1)
            if suffix.isdigit():
                return base
        return name

    def marker_world_position(obj) -> np.ndarray:
        translation = obj.matrix_world.translation
        return np.array([translation.x, translation.y, translation.z], dtype=float)

    template = get_fragment_template(fi.template_id)
    if template is None:
        raise ValueError(f"Template '{fi.template_id}' not found")

    rdkit_positions = np.asarray(
        [rdkit_atom_position(mol, int(atom_index)) for atom_index in fi.source_atom_indices],
        dtype=float,
    ) * float(scale)

    if rdkit_positions.size == 0:
        return make_transform_matrix(np.eye(3, dtype=float), np.zeros(3, dtype=float))

    template_name = template.blender_template_name or fi.template_id
    source_obj = bpy.data.objects.get(template_name)
    if source_obj is None:
        raise ValueError(f"Blender template root object '{template_name}' not found")

    all_marker_children = [
        child
        for child in source_obj.children_recursive
        if "mk_atom" in child.name
    ]
    if not all_marker_children:
        raise ValueError(f"No atom markers found under template '{template_name}'")

    marker_by_exact_name = {child.name: child for child in all_marker_children}
    marker_by_normalized_name = {
        normalize_marker_name(child.name): child for child in all_marker_children
    }

    requested_marker_names = list(getattr(template, "atom_marker_names", []) or [])
    ordered_markers: list[object] = []

    if requested_marker_names:
        for marker_name in requested_marker_names:
            marker_obj = marker_by_exact_name.get(marker_name)
            if marker_obj is None:
                marker_obj = marker_by_normalized_name.get(
                    normalize_marker_name(marker_name)
                )
            if marker_obj is None:
                raise ValueError(
                    f"Marker '{marker_name}' not found in template '{template_name}'"
                )
            ordered_markers.append(marker_obj)
    else:
        ordered_markers = sorted(
            all_marker_children,
            key=lambda obj: normalize_marker_name(obj.name),
        )

    if len(ordered_markers) != len(rdkit_positions):
        raise ValueError(
            f"Marker count mismatch for template '{template_name}': "
            f"{len(ordered_markers)} markers vs {len(rdkit_positions)} RDKit positions"
        )

    template_positions = np.asarray(
        [marker_world_position(obj) for obj in ordered_markers],
        dtype=float,
    )

    return best_fit_transform_matrix(template_positions, rdkit_positions)
