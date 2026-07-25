"""
molmol/site_assignment.py
-------------------------
Legacy contact-patch extraction utilities based on Blender vertex groups.

This module is retained only for backward compatibility and debugging of
older patch-based assets. It must not drive fragment placement, topology
resolution, or independent bond generation in the current marker-based pipeline.
"""

from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
from mathutils import Matrix

from .data_model import contact_patch


CONTACT_VERTEX_GROUP_PREFIX = "contact_"
_MIN_NORMAL_EPS = 1e-8


class ContactAssignmentError(RuntimeError):
    """Explicit error raised during contact patch extraction."""


def _ensure_mesh_object(obj) -> None:
    if obj is None:
        raise ContactAssignmentError("obj is None")
    if getattr(obj, "type", None) != "MESH":
        raise ContactAssignmentError(
            f"Object '{getattr(obj, 'name', '<unnamed>')}' is not a MESH"
        )
    if getattr(obj, "data", None) is None:
        raise ContactAssignmentError(f"Object '{obj.name}' has no mesh data")


def _world_vertex_positions(obj) -> np.ndarray:
    mesh = obj.data
    world = obj.matrix_world
    positions = np.empty((len(mesh.vertices), 3), dtype=float)
    for i, vertex in enumerate(mesh.vertices):
        world_position = world @ vertex.co
        positions[i] = (world_position.x, world_position.y, world_position.z)
    return positions


def _normal_matrix_3x3(obj) -> Matrix:
    world_3x3 = obj.matrix_world.to_3x3()
    try:
        return world_3x3.inverted().transposed()
    except Exception as exc:
        raise ContactAssignmentError(
            f"Cannot build normal matrix for object '{obj.name}': {exc}"
        )


def _vertex_group_index_by_name(obj, vertex_group_name: str) -> int:
    vertex_group = obj.vertex_groups.get(vertex_group_name)
    if vertex_group is None:
        raise ContactAssignmentError(
            f"Object '{obj.name}' has no vertex group '{vertex_group_name}'"
        )
    return int(vertex_group.index)


def _vertex_indices_for_group(obj, vertex_group_name: str) -> List[int]:
    group_index = _vertex_group_index_by_name(obj, vertex_group_name)
    indices: List[int] = []

    for vertex in obj.data.vertices:
        for group in vertex.groups:
            if group.group == group_index and group.weight > 0.0:
                indices.append(int(vertex.index))
                break

    if not indices:
        raise ContactAssignmentError(
            f"Vertex group '{vertex_group_name}' on object '{obj.name}' has no vertices"
        )

    return indices


def _world_vertex_normals(obj, vertex_indices: Iterable[int]) -> np.ndarray:
    mesh = obj.data
    normal_matrix = _normal_matrix_3x3(obj)
    normals = []

    for index in vertex_indices:
        local_normal = mesh.vertices[int(index)].normal
        world_normal = normal_matrix @ local_normal
        if world_normal.length > _MIN_NORMAL_EPS:
            world_normal.normalize()
        normals.append((world_normal.x, world_normal.y, world_normal.z))

    return np.asarray(normals, dtype=float)


def _safe_normalize(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= _MIN_NORMAL_EPS:
        raise ContactAssignmentError("Cannot normalize zero-length vector")
    return vector / norm





def _covariance_normal(sample_points: np.ndarray) -> np.ndarray:
    """
    Robust fallback: estimate patch normal by PCA.
    The selected direction is the eigenvector associated with the smallest variance.
    """
    if sample_points.shape[0] < 3:
        raise ContactAssignmentError(
            "At least 3 points are required to estimate a patch normal by covariance"
        )

    centered = sample_points - sample_points.mean(axis=0)
    covariance = centered.T @ centered
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    normal = eigenvectors[:, int(np.argmin(eigenvalues))]
    return _safe_normalize(normal)


def _mean_normal_with_fallback(
    obj,
    vertex_indices: List[int],
    sample_points: np.ndarray,
) -> np.ndarray:
    normals = _world_vertex_normals(obj, vertex_indices)
    if len(normals) == 0:
        return _covariance_normal(sample_points)

    mean_normal = normals.mean(axis=0)
    norm = float(np.linalg.norm(mean_normal))
    if norm > _MIN_NORMAL_EPS:
        return mean_normal / norm

    return _covariance_normal(sample_points)


def _disambiguate_normal_direction(sample_points: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """
    Make the normal orientation deterministic.

    Local geometry cannot always determine a true outward direction.
    To avoid random flips between sessions, enforce a stable convention:
    the component with the largest absolute value must be positive.
    """
    dominant_index = int(np.argmax(np.abs(normal)))
    if normal[dominant_index] < 0.0:
        return -normal
    return normal


def iter_contact_vertex_group_names(obj) -> List[str]:
    """
    Return all vertex groups following the 'contact_<side>_<index>' naming convention,
    sorted lexicographically.
    """
    _ensure_mesh_object(obj)
    names = [
        vertex_group.name
        for vertex_group in obj.vertex_groups
        if vertex_group.name.startswith(CONTACT_VERTEX_GROUP_PREFIX)
    ]
    names.sort()
    return names


def extract_contact_patch(obj, vertex_group_name: str) -> contact_patch:
    """
    Extract a complete contact_patch from a Blender mesh vertex group.

    Parameters
    ----------
    obj:
        Blender mesh object.
    vertex_group_name:
        Vertex group name, for example 'contact_a_0'.

    Returns
    -------
    contact_patch
        Patch with centroid, normal, and sample_points in world space.
    """
    _ensure_mesh_object(obj)

    vertex_indices = _vertex_indices_for_group(obj, vertex_group_name)
    world_positions = _world_vertex_positions(obj)
    sample_points = world_positions[vertex_indices]

    centroid = sample_points.mean(axis=0)
    normal = _mean_normal_with_fallback(obj, vertex_indices, sample_points)
    normal = _disambiguate_normal_direction(sample_points, normal)

    return contact_patch(
        vertex_group_name=vertex_group_name,
        centroid=centroid,
        normal=normal,
        sample_points=sample_points,
    )


def assign_contact_patches(obj) -> Dict[str, contact_patch]:
    """
    Scan all 'contact_*' vertex groups on the object and build
    the mapping {vertex_group_name: contact_patch}.

    Raises ContactAssignmentError for invalid objects or empty patches.
    """
    _ensure_mesh_object(obj)
    result: Dict[str, contact_patch] = {}
    for vertex_group_name in iter_contact_vertex_group_names(obj):
        result[vertex_group_name] = extract_contact_patch(obj, vertex_group_name)
    return result


def has_contact_patches(obj) -> bool:
    """Return True if the object contains at least one 'contact_*' vertex group."""
    return len(iter_contact_vertex_group_names(obj)) > 0


def patch_debug_dict(patch: contact_patch) -> dict:
    """
    Utility ready for logging/debug/UI without extra dependencies.
    """
    return {
        "vertex_group_name": patch.vertex_group_name,
        "side": patch.side,
        "index": patch.index,
        "centroid": patch.centroid.tolist(),
        "normal": patch.normal.tolist(),
        "sample_count": int(len(patch.sample_points)),
    }