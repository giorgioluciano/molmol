from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional,Set, Tuple

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector
from rdkit import Chem

from .assembly_planner import build_assembly_graph


from .fragment_library import (
    get_fragment_template,
    get_atom_template_name,
    get_bond_template_name,
    get_bond_template_collection_name,
)

from .fragment_recognizer import load_molecule
from .geometry_utils import (
    apply_transform_matrix,
    fit_fragment_pose,
    make_transform_matrix,
    rdkit_atom_position,
    rmsd, matrix4x4_to_blender,
)

from .smiles_builder import (
    duplicate_root_with_children,
    ensure_fragment_collection,
    get_template_source_object,
)
   
PIPELINE_EPS = 1e-8


class contact_pipeline_error(RuntimeError):
    """Explicit contact-pipeline error."""


ContactPipelineError = contact_pipeline_error


@dataclass
class placed_fragment:
    node_id: str
    blender_obj_name: str
    fragment_id: str
    node_kind: str = "fragment"
    connection_markers: Dict[str, "connection_marker"] = field(default_factory=dict)


@dataclass
class peg_hole_slot:
    slot_key: str
    kind: str  # "atom", "hole", "peg", "anchor"
    source_object_name: str
    parent_root_name: str
    location: np.ndarray


@dataclass
class bond_edge:
    node_a: str
    node_b: str
    marker_a: str | None = None
    marker_b: str | None = None
    bond_order: int = 1
    resolved: bool = False


@dataclass
class contact_assembly_graph:
    fragments: Dict[str, placed_fragment] = field(default_factory=dict)
    edges: List[bond_edge] = field(default_factory=list)

    def add_fragment(self, fragment: placed_fragment) -> None:
        self.fragments[fragment.node_id] = fragment

    def get_fragment(self, node_id: str) -> Optional[placed_fragment]:
        return self.fragments.get(node_id)

    def add_edge(self, edge: bond_edge) -> None:
        self.edges.append(edge)

@dataclass
class planner_node_binding:
    node_id: str
    blender_obj_name: str
    fragment_id: str | None = None


@dataclass
class unresolved_connection:
    node_a: str
    node_b: str
    bond_order: int = 1


@dataclass
class match_resolution_result:
    resolved_edges: List[bond_edge]
    rejected_connections: List[unresolved_connection]


@dataclass
class build_execution_result:
    planner_graph: object
    contact_graph: contact_assembly_graph
    match_result: match_resolution_result
    built_collection_name: str
    instantiated_objects: List[str]


@dataclass
class connection_marker:
    marker_name: str
    marker_type: str
    source_object_name: str
    local_position: np.ndarray
    world_position: np.ndarray

def _detect_marker_type(object_name: str) -> Optional[str]:
    """
    Detect the semantic marker type from a Blender object name.
    """
    name = str(object_name or "").strip().lower()
    if not name:
        return None

    marker_tokens = (
        "mk_atom",
        "mk_hole",
        "mk_peg",
        "mk_axis",
    )
    for token in marker_tokens:
        if token in name:
            return token

    return None

def _vector_to_numpy(value) -> np.ndarray:
    """
    Convert a Blender vector-like value to a float numpy array with shape (3,).
    """
    return np.array([float(value[0]), float(value[1]), float(value[2])], dtype=float)

def _blender_matrix_to_numpy(matrix) -> np.ndarray:
    """
    Convert a Blender 4x4 matrix to a numpy array.
    """
    return np.array([[float(matrix[r][c]) for c in range(4)] for r in range(4)], dtype=float)

def extract_connection_markers_from_hierarchy(root_obj) -> Dict[str, connection_marker]:
    """
    Extract all supported connection markers from a fragment hierarchy.

    Returned dictionary is keyed by the full Blender object name so marker names
    remain stable in logs and resolver diagnostics.
    """
    markers: Dict[str, connection_marker] = {}

    if root_obj is None:
        return markers

    supported_marker_types = {
        "mk_atom",
        "mk_hole",
        "mk_peg",
        "mk_axis",
    }

    for obj in _iter_object_hierarchy(root_obj):
        if obj is None or obj == root_obj:
            continue

        marker_type = _detect_marker_type(getattr(obj, "name", ""))
        if marker_type is None:
            continue
        if marker_type not in supported_marker_types:
            continue

        local_position = _vector_to_numpy(obj.location)
        world_position = _vector_to_numpy(obj.matrix_world.to_translation())

        marker = connection_marker(
            marker_name=str(obj.name),
            marker_type=marker_type,
            source_object_name=str(root_obj.name),
            local_position=local_position,
            world_position=world_position,
        )
        markers[marker.marker_name] = marker

    return markers

def debug_print_connection_markers(root_obj) -> None:
    """
    Print all detected connection markers for a fragment root.
    Useful during marker-based resolver debugging.
    """
    root_name = getattr(root_obj, "name", None)
    markers = extract_connection_markers_from_hierarchy(root_obj)

    print("MOLMOL_MARKERS_BEGIN", root_name, len(markers))

    for marker_name in sorted(markers.keys()):
        marker = markers[marker_name]
        print(
            "MOLMOL_MARKER",
            marker.source_object_name,
            marker.marker_name,
            marker.marker_type,
            tuple(float(x) for x in marker.local_position),
            tuple(float(x) for x in marker.world_position),
        )

    print("MOLMOL_MARKERS_END", root_name, len(markers))


def _get_attr(obj, name: str, default=None):
    return getattr(obj, name, default)

def _normalize_library_path(value) -> str:
    """
    Normalize library_path values that may arrive as:
    - a plain string
    - a tuple/list where the first item is the .blend path
    - None
    """
    if value is None:
        return ""

    if isinstance(value, (tuple, list)):
        if len(value) == 0:
            return ""
        value = value[0]

    return str(value).strip()



def _extract_graph_nodes(planner_graph) -> Dict[str, object]:
    nodes = _get_attr(planner_graph, "nodes", None)
    if nodes is None:
        raise contact_pipeline_error("Planner graph has no nodes attribute")

    if isinstance(nodes, dict):
        return {str(k): v for k, v in nodes.items()}

    result: Dict[str, object] = {}
    for node in nodes:
        node_id = _get_attr(node, "id", None)
        if not node_id:
            raise contact_pipeline_error("Planner node has no id attribute")
        result[str(node_id)] = node
    return result


def _extract_graph_edges(planner_graph) -> List[object]:
    edges = _get_attr(planner_graph, "edges", None)
    if edges is None:
        raise contact_pipeline_error("Planner graph has no edges attribute")
    return list(edges)


def _edge_endpoints(edge) -> Tuple[str, str]:
    for a_name, b_name in (
        ("node_a_id", "node_b_id"),
        ("node_a", "node_b"),
        ("source", "target"),
    ):
        a = _get_attr(edge, a_name, None)
        b = _get_attr(edge, b_name, None)
        if a is not None and b is not None:
            return str(a), str(b)

    raise contact_pipeline_error(
        f"Cannot determine endpoints for edge type {type(edge).__name__}"
    )



def _edge_bond_order(edge) -> int:
    for attr_name in ("bond_order", "order", "bond_type"):
        value = _get_attr(edge, attr_name, None)
        if value is not None:
            try:
                return int(value)
            except Exception:
                pass
    return 1


def _planner_node_fragment_id(planner_node) -> str:
    fragment_instance = _get_attr(planner_node, "fragment_instance", None)
    if fragment_instance is not None:
        value = _get_attr(fragment_instance, "template_id", None)
        if value:
            return str(value)

    linker_atom = _get_attr(planner_node, "linker_atom", None)
    if linker_atom is not None:
        value = _get_attr(linker_atom, "element", None)
        if value:
            return str(value)

    for attr_name in ("fragment_id", "template_id", "element"):
        value = _get_attr(planner_node, attr_name, None)
        if value:
            return str(value)

    return str(_get_attr(planner_node, "id", "unknown"))


def _planner_node_kind(planner_node) -> str:
    kind = _get_attr(planner_node, "kind", None)
    if kind == "fragment":
        return "fragment"
    if kind == "linker_atom":
        return "atom"

    node_type = _get_attr(planner_node, "type", None)
    if node_type == "fragment":
        return "fragment"
    if node_type == "atom":
        return "atom"

    linker_atom = _get_attr(planner_node, "linker_atom", None)
    if linker_atom is not None and _get_attr(linker_atom, "element", None) is not None:
        return "atom"

    return "fragment"


def _default_object_name(node_id: str, fragment_id: str) -> str:
    return f"{fragment_id}__{node_id}"


def _strip_slot_marker_tokens(slot_key: str) -> str:
    value = str(slot_key or "").strip().lower()
    if not value:
        return ""

    if ":" in value:
        _, value = value.split(":", 1)

    replacements = [
        "__atom__",
        "__hole__",
        "__peg__",
        "mk_atom",
        "mkatom",
        "mk_hole",
        "mkhole",
        "mk_peg",
        "mkpeg",
    ]
    for token in replacements:
        value = value.replace(token, "_")

    while "__" in value:
        value = value.replace("__", "_")

    return value.strip("_")


def _slot_family_signature(slot: peg_hole_slot) -> str:
    return _strip_slot_marker_tokens(slot.slot_key)


def _slot_is_compatible(slot_a: peg_hole_slot, slot_b: peg_hole_slot) -> bool:
    pair = {slot_a.kind, slot_b.kind}
    return pair in (
        {"peg", "hole"},
        {"atom", "hole"},
        {"atom", "peg"},
    )


def _slot_pair_priority(slot_a: peg_hole_slot, slot_b: peg_hole_slot) -> float:
    if not _slot_is_compatible(slot_a, slot_b):
        return -1.0

    kinds = {slot_a.kind, slot_b.kind}
    if kinds == {"peg", "hole"}:
        return 3.0
    if kinds == {"atom", "hole"}:
        return 2.0
    if kinds == {"atom", "peg"}:
        return 1.5
    return -1.0


def _count_residual_contacts(planner_graph, node_id: str) -> int:
    count = 0
    for edge in _extract_graph_edges(planner_graph):
        node_a, node_b = _edge_endpoints(edge)
        if node_a == node_id or node_b == node_id:
            count += 1
    return count



def append_object_from_blend_library(library_path: str, object_name: str):
    library_path = _normalize_library_path(library_path)
    object_name = str(object_name or "").strip()

    if not library_path or not object_name:
        return None

    existing = bpy.data.objects.get(object_name)
    if existing is not None:
        return existing

    try:
        with bpy.data.libraries.load(library_path, link=False) as (data_from, data_to):
            if object_name not in data_from.objects:
                return None
            data_to.objects = [object_name]
    except Exception as exc:
        raise ContactPipelineError(
            f"Cannot load object '{object_name}' from library '{library_path}': {exc}"
        )

    imported = bpy.data.objects.get(object_name)
    return imported

def _find_bond_template_pegs(root_obj) -> tuple[object, object]:
    peg_objects = []

    for obj in _iter_object_hierarchy(root_obj):
        if obj is None or obj == root_obj:
            continue
        marker_kind = _detect_marker_kind_from_name(getattr(obj, "name", ""))
        if marker_kind == "peg":
            peg_objects.append(obj)

    if len(peg_objects) != 2:
        raise ContactPipelineError(
            f"Bond template '{root_obj.name}' must contain exactly 2 peg markers, found {len(peg_objects)}"
        )

    peg_objects.sort(key=lambda obj: obj.name.lower())
    return peg_objects[0], peg_objects[1]


def _compute_root_transform_from_two_pegs(
    root_obj,
    peg_a_obj,
    peg_b_obj,
    target_a: Vector,
    target_b: Vector,
) -> Matrix:
    local_a = peg_a_obj.matrix_local.translation.copy()
    local_b = peg_b_obj.matrix_local.translation.copy()

    local_axis = local_b - local_a
    target_axis = target_b - target_a

    if local_axis.length <= PIPELINE_EPS:
        raise ContactPipelineError(
            f"Bond template '{root_obj.name}' has coincident local peg positions"
        )
    if target_axis.length <= PIPELINE_EPS:
        raise ContactPipelineError(
            f"Cannot place bond template '{root_obj.name}' on zero-length target axis"
        )

    local_dir = local_axis.normalized()
    target_dir = target_axis.normalized()

    rotation = _rotation_matrix_between_vectors(
        np.asarray(local_dir, dtype=float),
        np.asarray(target_dir, dtype=float),
    )
    rotation_blender = matrix4x4_to_blender(make_transform_matrix(rotation, np.zeros(3, dtype=float))).to_3x3()

    rotated_local_a = rotation_blender @ local_a
    translation = target_a - rotated_local_a

    transform = rotation_blender.to_4x4()
    transform.translation = translation
    return transform


def _duplicate_bond_template_instance(
    template_obj,
    target_collection,
    new_name: str,
):
    duplicate_result = duplicate_root_with_children(
        source_obj=template_obj,
        target_collection=target_collection,
        new_name=new_name,
    )
    if duplicate_result is None:
        raise ContactPipelineError(
            f"duplicate_root_with_children returned None for bond template '{template_obj.name}'"
        )

    root_obj = duplicate_result[0] if isinstance(duplicate_result, tuple) else duplicate_result
    if root_obj is None:
        raise ContactPipelineError(
            f"duplicate_root_with_children returned no root object for bond template '{template_obj.name}'"
        )
    return root_obj


def _build_internal_template_bond_pairs(
    mol,
    planner_graph,
    binding_map: Dict[str, planner_node_binding],
) -> set[tuple[int, int]]:
    """
    Return the set of RDKit atom-index bond pairs already covered by placed fragment templates.

    Only bonds explicitly declared as internal_bond_pairs in the corresponding FragmentTemplate
    are considered covered. All other bonds remain candidates for independent cylinders.
    """
    covered_pairs: set[tuple[int, int]] = set()
    nodes = _extract_graph_nodes(planner_graph)

    for node_id, binding in binding_map.items():
        planner_node = nodes.get(node_id)
        if planner_node is None:
            continue
        if _planner_node_kind(planner_node) != "fragment":
            continue

        fragment_instance = getattr(planner_node, "fragment_instance", None)
        if fragment_instance is None:
            continue

        source_atom_indices = getattr(fragment_instance, "source_atom_indices", None)
        if not source_atom_indices:
            continue

        fragment_id = _planner_node_fragment_id(planner_node)
        try:
            fragment_template = get_fragment_template(fragment_id)
        except KeyError:
            continue

        if not fragment_template.internal_bond_pairs:
            # Nessun bond interno dichiarato: questo frammento non copre bond via template mesh
            continue

        mapping = _build_template_local_to_rdkit_atom_map(
            mol=mol,
            fragment_template=fragment_template,
            source_atom_indices=list(source_atom_indices),
        )

        for local_a, local_b in fragment_template.internal_bond_pairs:
            rdkit_a = mapping.get(local_a)
            rdkit_b = mapping.get(local_b)
            if rdkit_a is None or rdkit_b is None:
                continue
            covered_pairs.add(tuple(sorted((int(rdkit_a), int(rdkit_b)))))

    return covered_pairs


def _append_fragment_template_collection_if_needed(
    fragment_id: str,
    library_path: str,
) -> None:
    """
    Ensure the Blender collection containing the fragment template hierarchy
    is loaded into bpy.data before resolving the root object.
    """
    library_path = _normalize_library_path(library_path)
    if not library_path:
        return

    try:
        template = get_fragment_template(fragment_id)
    except Exception as exc:
        raise ContactPipelineError(
            f"Cannot resolve fragment template for fragment_id '{fragment_id}': {exc}"
        )

    collection_name = str(template.blender_collection_name or "").strip()
    if not collection_name:
        return

    imported = append_fragment_collection_from_library(
        library_path=library_path,
        collection_name=collection_name,
    )
    if imported is None:
        raise ContactPipelineError(
            f"Cannot import fragment collection '{collection_name}' "
            f"for fragment_id '{fragment_id}' from '{library_path}'"
        )

    print(
        "MOLMOL_FRAGMENT_COLLECTION_READY",
        fragment_id,
        collection_name,
        template.blender_template_name,
    )

def resolve_template_source_object(source_name: str, library_path: str):
    """
    Resolve a template source object from the current file or from the library.

    Priority:
    1. Object already present in current file
    2. Import owning collection and search inside it
    3. Fallback single-object append by exact object name
    """
    import bpy
    from .fragment_library import FRAGMENT_LIBRARY, ATOM_TEMPLATE_COLLECTIONS, ATOM_TEMPLATE_NAMES

    library_path = _normalize_library_path(library_path)
    source_name = str(source_name or "").strip()

    if not source_name:
        raise ContactPipelineError("resolve_template_source_object: source_name is empty")

    source_obj = get_template_source_object(source_name)
    if source_obj is not None:
        print("MOLMOL_TEMPLATE_LOOKUP", source_name, "FROM_CURRENT_FILE")
        return source_obj

    print("MOLMOL_TEMPLATE_LOOKUP_FAILED", source_name)

    if not library_path:
        raise ContactPipelineError(
            f"Template {source_name!r} not in current file and no library_path provided"
        )

    collection_name = None

    for frag in FRAGMENT_LIBRARY.values():
        if frag.blender_template_name == source_name:
            collection_name = frag.blender_collection_name
            break

    if collection_name is None:
        for key, object_name in ATOM_TEMPLATE_NAMES.items():
            if object_name == source_name:
                collection_name = ATOM_TEMPLATE_COLLECTIONS.get(key)
                break

    if collection_name:
        print("MOLMOL_LIBRARY_IMPORT_COLLECTION", collection_name, "for", source_name)
        imported_collection = append_fragment_collection_from_library(
            library_path=library_path,
            collection_name=collection_name,
        )

        if imported_collection is not None:
            bpy.context.view_layer.update()

            for obj in imported_collection.all_objects:
                if obj.name == source_name:
                    print("MOLMOL_TEMPLATE_LOOKUP", source_name, "FROM_IMPORTED_COLLECTION_EXACT")
                    return obj

            for obj in imported_collection.all_objects:
                if obj.name.startswith(source_name):
                    print(
                        "MOLMOL_TEMPLATE_LOOKUP",
                        source_name,
                        "FROM_IMPORTED_COLLECTION_PREFIX",
                        obj.name,
                    )
                    return obj

    source_obj = get_template_source_object(source_name)
    if source_obj is not None:
        print("MOLMOL_TEMPLATE_LOOKUP", source_name, "FROM_LIBRARY_COLLECTION_POSTCHECK")
        return source_obj

    source_obj = append_object_from_blend_library(library_path, source_name)
    if source_obj is not None:
        print("MOLMOL_TEMPLATE_LOOKUP", source_name, "FROM_LIBRARY_OBJECT_FALLBACK")
        return source_obj

    available_candidates = [
        obj.name for obj in bpy.data.objects if source_name.lower() in obj.name.lower()
    ]
    print("MOLMOL_TEMPLATE_LOOKUP_FAILED_FINAL", source_name, f"candidates={available_candidates}")

    raise ContactPipelineError(
        f"Template {source_name!r} not found in current file or library {library_path!r}"
    )

def instantiate_planner_graph(
    planner_graph,
    collection_name: str = "MolMol_ContactBuild",
    library_path: str = "",
) -> Tuple[List[planner_node_binding], List[str], str]:
    library_path = _normalize_library_path(library_path)

    nodes = _extract_graph_nodes(planner_graph)
    collection = ensure_fragment_collection(collection_name)
    bindings: List[planner_node_binding] = []
    instantiated_objects: List[str] = []

    for node_id, planner_node in nodes.items():
        fragment_id = _planner_node_fragment_id(planner_node)
        object_kind = _planner_node_kind(planner_node)
        blender_obj_name = _default_object_name(node_id, fragment_id)

        existing = bpy.data.objects.get(blender_obj_name)
        if existing is not None:
            bpy.data.objects.remove(existing, do_unlink=True)

        if object_kind == "atom":
            linker_atom = _get_attr(planner_node, "linker_atom", None)
            if linker_atom is None:
                raise ContactPipelineError(
                    f"Planner node '{node_id}' is classified as atom but has no linker_atom"
                )

            element = str(_get_attr(linker_atom, "element", "C"))
            residual_contact_count = _count_residual_contacts(planner_graph, node_id)
            source_name = get_atom_template_name(element, residual_contact_count)
        else:
            if fragment_id == "bond_single":
                source_name = get_bond_template_name(1)
            else:
                template = get_fragment_template(fragment_id)
                _append_fragment_template_collection_if_needed(fragment_id, library_path)
                source_name = template.blender_template_name or fragment_id

        source_obj = resolve_template_source_object(
            source_name=source_name,
            library_path=library_path,
        )

        print("MOLMOL_SOURCE_OBJECT_FOR_DUPLICATION", node_id, fragment_id, source_obj.name)
        debug_print_connection_markers(source_obj)

        if source_obj is None:
            raise ContactPipelineError(
                f"No template source object for fragment_id '{fragment_id}' "
                f"(source_name='{source_name}')"
            )

        print("MOLMOL_TEMPLATE_LOOKUP", node_id, fragment_id, source_name)

        duplicate_result = duplicate_root_with_children(
            source_obj=source_obj,
            target_collection=collection,
            new_name=blender_obj_name,
        )

        if duplicate_result is None:
            raise ContactPipelineError(
                f"duplicate_root_with_children returned None for node '{node_id}'"
            )

        root_obj = duplicate_result[0] if isinstance(duplicate_result, tuple) else duplicate_result
        if root_obj is None:
            raise ContactPipelineError(
                f"duplicate_root_with_children returned no root object for node '{node_id}'"
            )


        bindings.append(
            planner_node_binding(
                node_id=node_id,
                blender_obj_name=root_obj.name,
                fragment_id=fragment_id,
            )
        )
        instantiated_objects.append(root_obj.name)
        print("MOLMOL_DUPLICATED_ROOT", node_id, fragment_id, root_obj.name)
        debug_print_connection_markers(root_obj)

    return bindings, instantiated_objects, collection.name


def build_contact_assembly_graph(
    planner_graph,
    bindings: Iterable[planner_node_binding],
) -> Tuple[contact_assembly_graph, List[unresolved_connection]]:
    nodes = _extract_graph_nodes(planner_graph)
    edges = _extract_graph_edges(planner_graph)

    binding_map = {binding.node_id: binding for binding in bindings}
    missing = [node_id for node_id in nodes if node_id not in binding_map]
    if missing:
        raise contact_pipeline_error(f"Missing bindings for node ids: {', '.join(sorted(missing))}")

    graph = contact_assembly_graph()
    unresolved_connections: List[unresolved_connection] = []

    for node_id, planner_node in nodes.items():
        binding = binding_map[node_id]
        blender_obj = bpy.data.objects.get(binding.blender_obj_name)
        if blender_obj is None:
            raise contact_pipeline_error(f"Object '{binding.blender_obj_name}' not found")

        connection_markers = extract_connection_markers_from_hierarchy(blender_obj)
        debug_print_connection_markers(blender_obj)

        graph.add_fragment(
            placed_fragment(
                node_id=node_id,
                blender_obj_name=blender_obj.name,
                fragment_id=binding.fragment_id or _planner_node_fragment_id(planner_node),
                node_kind=_planner_node_kind(planner_node),
                connection_markers=connection_markers,
            )
        )

    for edge in edges:
        node_a, node_b = _edge_endpoints(edge)
        unresolved_connections.append(
            unresolved_connection(
                node_a=node_a,
                node_b=node_b,
                bond_order=_edge_bond_order(edge),
            )
        )

    return graph, unresolved_connections

def create_resolved_bond_edge(
    graph: contact_assembly_graph,
    node_a: str,
    marker_a: str,
    node_b: str,
    marker_b: str,
    bond_order: int = 1,
) -> bond_edge:
    frag_a = graph.get_fragment(node_a)
    frag_b = graph.get_fragment(node_b)

    if frag_a is None or frag_b is None:
        raise ContactPipelineError("Unknown fragment reference in resolved edge creation")
    if marker_a not in frag_a.connection_markers:
        raise ContactPipelineError(f"Fragment '{node_a}' has no marker '{marker_a}'")
    if marker_b not in frag_b.connection_markers:
        raise ContactPipelineError(f"Fragment '{node_b}' has no marker '{marker_b}'")

    return bond_edge(
        node_a=node_a,
        node_b=node_b,
        marker_a=marker_a,
        marker_b=marker_b,
        bond_order=int(bond_order),
        resolved=True,
    )


def _marker_atom_key(marker_name: str) -> str:
    """
    Extract a stable atom-side key from a marker object name.

    Examples:
    - c_o_o_c_mk_atom      -> c
    - c_o_o_o2_mk_hole_1   -> o2
    - benzene_mk_atom_C_0  -> c_0
    """
    name = str(marker_name or "").strip().lower()
    if not name:
        return ""

    for token in ("_mk_atom", "_mk_hole", "_mk_peg", "_mk_axis"):
        if token in name:
            prefix = name.split(token, 1)[0]
            parts = [part for part in prefix.split("_") if part]
            if not parts:
                return ""
            if len(parts) >= 2 and parts[-1].isdigit():
                return f"{parts[-2]}_{parts[-1]}"
            return parts[-1]

    return ""


def _marker_pair_priority(marker_a: connection_marker, marker_b: connection_marker) -> int:
    """
    Lower is better.
    """
    pair = {marker_a.marker_type, marker_b.marker_type}

    if pair == {"mk_hole", "mk_atom"}:
        return 0
    if pair == {"mk_peg", "mk_hole"}:
        return 1
    if pair == {"mk_atom", "mk_atom"}:
        return 2

    return 100


def _marker_types_compatible(marker_a: connection_marker, marker_b: connection_marker) -> bool:
    return _marker_pair_priority(marker_a, marker_b) < 100


def _choose_connection_markers(
    frag_a: placed_fragment,
    frag_b: placed_fragment,
) -> Optional[Tuple[str, str]]:
    """
    Select the best marker pair between two fragments.

    Strategy:
    1. Match by atom-side key extracted from marker names.
    2. Prefer mk_hole<->mk_atom.
    3. Fallback to any compatible pair.
    """
    candidates: List[Tuple[int, str, str]] = []

    markers_a = list(frag_a.connection_markers.values())
    markers_b = list(frag_b.connection_markers.values())

    for marker_a in markers_a:
        for marker_b in markers_b:
            if not _marker_types_compatible(marker_a, marker_b):
                continue

            key_a = _marker_atom_key(marker_a.marker_name)
            key_b = _marker_atom_key(marker_b.marker_name)

            same_atom_key = bool(key_a) and bool(key_b) and key_a == key_b
            priority = _marker_pair_priority(marker_a, marker_b)

            if same_atom_key:
                candidates.append((priority, marker_a.marker_name, marker_b.marker_name))
            else:
                candidates.append((priority + 10, marker_a.marker_name, marker_b.marker_name))

    if not candidates:
        return None

    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    _, marker_name_a, marker_name_b = candidates[0]
    return marker_name_a, marker_name_b


def resolve_all_connections_from_markers(
    graph: contact_assembly_graph,
    unresolved_connections: Iterable[unresolved_connection],
) -> match_resolution_result:
    """
    Resolve graph connections using marker compatibility only.

    Important:
    - this function resolves topology, not fragment placement
    - fragment transforms from the initial RDKit / balls+stick layout must remain untouched
    """
    unresolved_connections = list(unresolved_connections)

    resolved_edges: List[bond_edge] = []
    rejected_connections: List[unresolved_connection] = []
    used_marker_names: set[str] = set()

    print("MOLMOL_MARKER_RESOLUTION_BEGIN", len(unresolved_connections))

    for connection in unresolved_connections:
        frag_a = graph.get_fragment(connection.node_a)
        frag_b = graph.get_fragment(connection.node_b)

        if frag_a is None or frag_b is None:
            rejected_connections.append(connection)
            print(
                "MOLMOL_MARKER_CONNECTION_REJECTED",
                connection.node_a,
                connection.node_b,
                "missing_fragment",
            )
            continue

        chosen_pair = _choose_connection_markers(frag_a, frag_b)
        if chosen_pair is None:
            rejected_connections.append(connection)
            print(
                "MOLMOL_MARKER_CONNECTION_REJECTED",
                connection.node_a,
                connection.node_b,
                "no_compatible_markers",
            )
            continue

        marker_a, marker_b = chosen_pair

        if marker_a in used_marker_names or marker_b in used_marker_names:
            rejected_connections.append(connection)
            print(
                "MOLMOL_MARKER_CONNECTION_REJECTED",
                connection.node_a,
                connection.node_b,
                "marker_already_used",
                marker_a,
                marker_b,
            )
            continue

        edge = create_resolved_bond_edge(
            graph=graph,
            node_a=connection.node_a,
            marker_a=marker_a,
            node_b=connection.node_b,
            marker_b=marker_b,
            bond_order=connection.bond_order,
        )

        graph.add_edge(edge)
        resolved_edges.append(edge)
        used_marker_names.add(marker_a)
        used_marker_names.add(marker_b)

        print(
            "MOLMOL_MARKER_CONNECTION_RESOLVED",
            connection.node_a,
            marker_a,
            "->",
            connection.node_b,
            marker_b,
            "bond_order",
            int(connection.bond_order),
        )

    print(
        "MOLMOL_MARKER_RESOLUTION_END",
        "resolved",
        len(resolved_edges),
        "rejected",
        len(rejected_connections),
    )

    return match_resolution_result(
        resolved_edges=resolved_edges,
        rejected_connections=rejected_connections,
    )


def resolve_all_connections(
    graph: contact_assembly_graph,
    unresolved_connections: Iterable[unresolved_connection],
) -> match_resolution_result:
    """
    Resolve all planner connections using the marker-based topology resolver only.
    """
    print("MOLMOL_CONNECTION_RESOLVER", "marker_only")
    return resolve_all_connections_from_markers(
        graph=graph,
        unresolved_connections=unresolved_connections,
    )


def _compute_mean_rdkit_single_bond_length(mol) -> float:
    distances: List[float] = []

    for bond in mol.GetBonds():
        bond_order = int(round(float(bond.GetBondTypeAsDouble())))
        bond_order = max(1, min(3, bond_order))
        if bond_order != 1:
            continue

        atom_idx_a = int(bond.GetBeginAtomIdx())
        atom_idx_b = int(bond.GetEndAtomIdx())
        pos_a = np.asarray(rdkit_atom_position(mol, atom_idx_a), dtype=float)
        pos_b = np.asarray(rdkit_atom_position(mol, atom_idx_b), dtype=float)
        distance = float(np.linalg.norm(pos_b - pos_a))

        if distance > PIPELINE_EPS:
            distances.append(distance)

    if not distances:
        raise ContactPipelineError("No RDKit single bonds found for layout scaling")

    mean_distance = float(np.mean(distances))
    if mean_distance <= PIPELINE_EPS:
        raise ContactPipelineError("Mean RDKit single-bond length is zero")

    return mean_distance

def compute_global_layout_scale(
    mol,
    template_single_bond_length: float = 1.62,
    extra_offset_factor: float = 0.10,
    fallback_rdkit_single_bond_length: float = 1.40,
) -> float:
    target_single_bond_length = float(template_single_bond_length) * (
        1.0 + float(extra_offset_factor)
    )
    if target_single_bond_length <= PIPELINE_EPS:
        raise ContactPipelineError("Target single-bond length must be positive")

    try:
        mean_rdkit_single_bond_length = _compute_mean_rdkit_single_bond_length(mol)
    except ContactPipelineError as exc:
        if "No RDKit single bonds found for layout scaling" not in str(exc):
            raise
        mean_rdkit_single_bond_length = float(fallback_rdkit_single_bond_length)

    if mean_rdkit_single_bond_length <= PIPELINE_EPS:
        raise ContactPipelineError("Fallback RDKit single-bond length must be positive")

    return float(target_single_bond_length / mean_rdkit_single_bond_length)


def _representative_rdkit_position(mol, planner_node, scale: float) -> Optional[np.ndarray]:
    linker_atom = _get_attr(planner_node, "linker_atom", None)
    if linker_atom is not None:
        source_atom_index = _get_attr(linker_atom, "source_atom_index", None)
        if source_atom_index is not None:
            return np.asarray(rdkit_atom_position(mol, int(source_atom_index)), dtype=float) * scale

    fragment_instance = _get_attr(planner_node, "fragment_instance", None)
    if fragment_instance is not None:
        source_atom_indices = _get_attr(fragment_instance, "source_atom_indices", None)
        if source_atom_indices:
            points = [
                np.asarray(rdkit_atom_position(mol, int(atom_index)), dtype=float)
                for atom_index in source_atom_indices
            ]
            centroid = np.mean(np.asarray(points, dtype=float), axis=0)
            return centroid * scale

    return None

def _iter_object_hierarchy(root_obj):
    if root_obj is None:
        return

    stack = [root_obj]
    while stack:
        current_obj = stack.pop()
        yield current_obj
        direct_children = list(getattr(current_obj, "children", []) or [])
        stack.extend(reversed(direct_children))


def _safe_normalize(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm <= PIPELINE_EPS:
        raise ContactPipelineError("Cannot normalize a near-zero vector")
    return vector / norm


def _rotation_matrix_between_vectors(
    source_vector: np.ndarray,
    target_vector: np.ndarray,
) -> np.ndarray:
    source = _safe_normalize(source_vector)
    target = _safe_normalize(target_vector)

    dot_value = float(np.clip(np.dot(source, target), -1.0, 1.0))

    if dot_value >= 1.0 - PIPELINE_EPS:
        return np.eye(3, dtype=float)

    if dot_value <= -1.0 + PIPELINE_EPS:
        fallback_axis = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(float(np.dot(source, fallback_axis))) >= 0.9:
            fallback_axis = np.array([0.0, 1.0, 0.0], dtype=float)

        rotation_axis = np.cross(source, fallback_axis)
        rotation_axis = _safe_normalize(rotation_axis)

        x, y, z = rotation_axis
        angle = np.pi
        c = float(np.cos(angle))
        s = float(np.sin(angle))
        t = 1.0 - c

        return np.array(
            [
                [t * x * x + c,     t * x * y - s * z, t * x * z + s * y],
                [t * x * y + s * z, t * y * y + c,     t * y * z - s * x],
                [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
            ],
            dtype=float,
        )

    rotation_axis = np.cross(source, target)
    axis_norm = float(np.linalg.norm(rotation_axis))
    if axis_norm <= PIPELINE_EPS:
        return np.eye(3, dtype=float)

    rotation_axis = rotation_axis / axis_norm
    x, y, z = rotation_axis

    angle = float(np.arccos(dot_value))
    c = float(np.cos(angle))
    s = float(np.sin(angle))
    t = 1.0 - c

    return np.array(
        [
            [t * x * x + c,     t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c,     t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ],
        dtype=float,
    )


def _neighbor_node_ids(planner_graph, node_id: str) -> List[str]:
    normalized_node_id = str(node_id)
    neighbors: List[str] = []
    seen: set[str] = set()

    for edge in _extract_graph_edges(planner_graph):
        node_a, node_b = _edge_endpoints(edge)

        if node_a == normalized_node_id and node_b not in seen:
            neighbors.append(node_b)
            seen.add(node_b)
        elif node_b == normalized_node_id and node_a not in seen:
            neighbors.append(node_a)
            seen.add(node_a)

    return neighbors


def _place_single_atom_fragment_with_axis(
    blender_obj,
    planner_graph,
    planner_node,
    mol,
    scale: float,
    axis_marker_name: str,
) -> None:
    fragment_instance = _get_attr(planner_node, "fragment_instance", None)
    if fragment_instance is None:
        raise ContactPipelineError("Single-atom fragment has no fragment_instance")

    source_atom_indices = list(_get_attr(fragment_instance, "source_atom_indices", []) or [])
    if len(source_atom_indices) != 1:
        raise ContactPipelineError(
            "_place_single_atom_fragment_with_axis requires exactly one source atom"
        )

    atom_index = int(source_atom_indices[0])
    atom_position = np.asarray(rdkit_atom_position(mol, atom_index), dtype=float) * scale

    node_id = str(_get_attr(planner_node, "id", ""))
    neighbor_ids = _neighbor_node_ids(planner_graph, node_id)

    if not neighbor_ids:
        transform = make_transform_matrix(np.eye(3, dtype=float), atom_position)
        blender_obj.matrix_world = matrix4x4_to_blender(transform)
        return

    nodes = _extract_graph_nodes(planner_graph)
    target_direction = None

    for neighbor_node_id in neighbor_ids:
        neighbor_node = nodes.get(neighbor_node_id)
        if neighbor_node is None:
            continue

        neighbor_position = _representative_rdkit_position(mol, neighbor_node, scale)
        if neighbor_position is None:
            continue

        direction = np.asarray(neighbor_position, dtype=float) - atom_position
        if np.linalg.norm(direction) > PIPELINE_EPS:
            target_direction = _safe_normalize(direction)
            break

    if target_direction is None:
        transform = make_transform_matrix(np.eye(3, dtype=float), atom_position)
        blender_obj.matrix_world = matrix4x4_to_blender(transform)
        return

    axis_marker_obj = _find_child_object_by_name(blender_obj, axis_marker_name)
    if axis_marker_obj is None:
        raise ContactPipelineError(
            f"Axis marker {axis_marker_name!r} not found under fragment root {blender_obj.name!r}"
        )

    local_axis = np.asarray(axis_marker_obj.location, dtype=float)
    if np.linalg.norm(local_axis) <= PIPELINE_EPS:
        raise ContactPipelineError(
            f"Axis marker {axis_marker_name!r} has zero local offset in {blender_obj.name!r}"
        )

    attachment_axis_local = _safe_normalize(local_axis)
    rotation = _rotation_matrix_between_vectors(attachment_axis_local, target_direction)
    transform = make_transform_matrix(rotation, atom_position)
    blender_obj.matrix_world = matrix4x4_to_blender(transform)


def _get_root_object_from_name(object_name: str):
    normalized_name = str(object_name or "").strip()
    if not normalized_name:
        return None
    return bpy.data.objects.get(normalized_name)


def _collect_marker_objects_by_type(root_obj, marker_type: str) -> List[object]:
    result: List[object] = []

    if root_obj is None:
        return result

    normalized_marker_type = str(marker_type or "").strip().lower()
    if not normalized_marker_type:
        return result

    for child_obj in _iter_object_hierarchy(root_obj):
        if child_obj == root_obj:
            continue

        child_name = str(getattr(child_obj, "name", "") or "").strip()
        detected_type = _detect_marker_type(child_name)
        if detected_type == normalized_marker_type:
            result.append(child_obj)

    result.sort(key=lambda obj: str(getattr(obj, "name", "") or ""))
    return result


def _marker_world_position(marker_obj) -> Vector:
    return marker_obj.matrix_world.to_translation().copy()


def _choose_best_hole_pair_for_bond(
    peg_positions_world: List[Vector],
    available_hole_markers: List[object],
) -> Optional[Tuple[object, object]]:
    if len(peg_positions_world) != 2:
        return None

    if len(available_hole_markers) < 2:
        return None

    best_pair: Optional[Tuple[object, object]] = None
    best_score: Optional[float] = None

    peg_a = peg_positions_world[0]
    peg_b = peg_positions_world[1]

    for index_a in range(len(available_hole_markers)):
        for index_b in range(index_a + 1, len(available_hole_markers)):
            hole_a = available_hole_markers[index_a]
            hole_b = available_hole_markers[index_b]

            hole_pos_a = _marker_world_position(hole_a)
            hole_pos_b = _marker_world_position(hole_b)

            direct_score = (peg_a - hole_pos_a).length + (peg_b - hole_pos_b).length
            swapped_score = (peg_a - hole_pos_b).length + (peg_b - hole_pos_a).length

            pair_score = min(direct_score, swapped_score)

            if best_score is None or pair_score < best_score:
                best_score = pair_score
                if direct_score <= swapped_score:
                    best_pair = (hole_a, hole_b)
                else:
                    best_pair = (hole_b, hole_a)

    return best_pair


def _compute_rigid_bond_fix_transform(
    bond_root_obj,
    peg_obj_a,
    peg_obj_b,
    target_hole_a,
    target_hole_b,
) -> Optional[Matrix]:
    if bond_root_obj is None:
        return None
    if peg_obj_a is None or peg_obj_b is None:
        return None
    if target_hole_a is None or target_hole_b is None:
        return None

    source_peg_a = _marker_world_position(peg_obj_a)
    source_peg_b = _marker_world_position(peg_obj_b)
    target_pos_a = _marker_world_position(target_hole_a)
    target_pos_b = _marker_world_position(target_hole_b)

    source_axis = source_peg_b - source_peg_a
    target_axis = target_pos_b - target_pos_a

    if source_axis.length <= 1e-8:
        return None
    if target_axis.length <= 1e-8:
        return None

    source_midpoint = (source_peg_a + source_peg_b) * 0.5
    target_midpoint = (target_pos_a + target_pos_b) * 0.5

    rotation = source_axis.normalized().rotation_difference(target_axis.normalized())
    rotation_matrix = rotation.to_matrix().to_4x4()

    translate_to_origin = Matrix.Translation(-source_midpoint)
    translate_to_target = Matrix.Translation(target_midpoint)

    transform = translate_to_target @ rotation_matrix @ translate_to_origin
    return transform


def _apply_world_transform_to_root(root_obj, transform_matrix: Matrix) -> None:
    if root_obj is None:
        return
    root_obj.matrix_world = transform_matrix @ root_obj.matrix_world


def _collect_available_hole_markers_from_contact_graph(
    contact_graph: contact_assembly_graph,
    excluded_root_names: Optional[Set[str]] = None,
) -> List[object]:
    excluded_root_names = set(excluded_root_names or [])
    hole_markers: List[object] = []

    for fragment in contact_graph.fragments.values():
        root_name = str(fragment.blender_obj_name or "").strip()
        if not root_name:
            continue
        if root_name in excluded_root_names:
            continue

        root_obj = bpy.data.objects.get(root_name)
        if root_obj is None:
            continue

        root_holes = _collect_marker_objects_by_type(root_obj, "mk_hole")
        hole_markers.extend(root_holes)

    return hole_markers


def try_to_fix_independent_bonds(
    contact_graph: contact_assembly_graph,
    independent_bond_names: Iterable[str],
) -> List[str]:
    """
    Try to rigidly reposition already-created independent bond roots so their two
    bond pegs lie on the line between the two most plausible hole markers in the
    current scene.

    Rules:
    - no scaling
    - only translation + rotation
    - only bonds with exactly two mk_peg markers are processed
    """
    fixed_bond_names: List[str] = []
    independent_bond_names = [str(name).strip() for name in independent_bond_names if str(name).strip()]

    if not independent_bond_names:
        return fixed_bond_names

    available_hole_markers = _collect_available_hole_markers_from_contact_graph(contact_graph)

    print("MOLMOL_TRY_FIX_BONDS_BEGIN", len(independent_bond_names), "holes", len(available_hole_markers))

    for bond_name in independent_bond_names:
        bond_root_obj = _get_root_object_from_name(bond_name)
        if bond_root_obj is None:
            print("MOLMOL_TRY_FIX_BOND_SKIPPED", bond_name, "missing_bond_root")
            continue

        peg_markers = _collect_marker_objects_by_type(bond_root_obj, "mk_peg")
        if len(peg_markers) != 2:
            print("MOLMOL_TRY_FIX_BOND_SKIPPED", bond_name, "expected_two_pegs", len(peg_markers))
            continue

        peg_positions_world = [
            _marker_world_position(peg_markers[0]),
            _marker_world_position(peg_markers[1]),
        ]

        chosen_holes = _choose_best_hole_pair_for_bond(
            peg_positions_world=peg_positions_world,
            available_hole_markers=available_hole_markers,
        )
        if chosen_holes is None:
            print("MOLMOL_TRY_FIX_BOND_SKIPPED", bond_name, "no_hole_pair")
            continue

        hole_a, hole_b = chosen_holes

        transform = _compute_rigid_bond_fix_transform(
            bond_root_obj=bond_root_obj,
            peg_obj_a=peg_markers[0],
            peg_obj_b=peg_markers[1],
            target_hole_a=hole_a,
            target_hole_b=hole_b,
        )
        if transform is None:
            print("MOLMOL_TRY_FIX_BOND_SKIPPED", bond_name, "no_transform")
            continue

        _apply_world_transform_to_root(bond_root_obj, transform)
        bpy.context.view_layer.update()

        fixed_bond_names.append(bond_root_obj.name)

        print(
            "MOLMOL_TRY_FIX_BOND_APPLIED",
            bond_root_obj.name,
            "peg_a",
            peg_markers[0].name,
            "peg_b",
            peg_markers[1].name,
            "hole_a",
            hole_a.name,
            "hole_b",
            hole_b.name,
        )

    print("MOLMOL_TRY_FIX_BONDS_END", "fixed", len(fixed_bond_names))
    return fixed_bond_names


def _get_top_level_root_object(obj):
    current_obj = obj
    if current_obj is None:
        return None

    while getattr(current_obj, "parent", None) is not None:
        current_obj = current_obj.parent

    return current_obj


def get_selected_root_objects(context) -> List[object]:
    selected_objects = list(getattr(context, "selected_objects", []) or [])
    selected_roots: List[object] = []
    seen_root_names: Set[str] = set()

    for selected_obj in selected_objects:
        root_obj = _get_top_level_root_object(selected_obj)
        if root_obj is None:
            continue

        root_name = str(getattr(root_obj, "name", "") or "").strip()
        if not root_name:
            continue
        if root_name in seen_root_names:
            continue

        seen_root_names.add(root_name)
        selected_roots.append(root_obj)

    selected_roots.sort(key=lambda obj: str(getattr(obj, "name", "") or ""))
    return selected_roots


def _build_contact_graph_from_root_objects(root_objects: Iterable[object]) -> contact_assembly_graph:
    graph = contact_assembly_graph()

    for root_obj in root_objects:
        if root_obj is None:
            continue

        root_name = str(getattr(root_obj, "name", "") or "").strip()
        if not root_name:
            continue

        connection_markers = extract_connection_markers_from_hierarchy(root_obj)

        graph.add_fragment(
            placed_fragment(
                node_id=root_name,
                blender_obj_name=root_name,
                fragment_id=root_name,
                node_kind="fragment",
                connection_markers=connection_markers,
            )
        )

    return graph


def get_independent_bond_root_names_from_root_objects(root_objects: Iterable[object]) -> List[str]:
    bond_names: List[str] = []
    seen_bond_names: Set[str] = set()

    for root_obj in root_objects:
        if root_obj is None:
            continue

        root_name = str(getattr(root_obj, "name", "") or "").strip()
        if not root_name:
            continue
        if not root_name.startswith("bond_independent_"):
            continue
        if root_name in seen_bond_names:
            continue

        seen_bond_names.add(root_name)
        bond_names.append(root_name)

    bond_names.sort()
    return bond_names


def try_to_fix_independent_bonds_for_selected_objects(context) -> List[str]:
    """
    Selection-based entry point for the Blender operator.

    Rules:
    - operate only on selected root pieces
    - selected independent bonds are the only movable targets
    - selected non-bond roots provide the candidate hole markers
    """
    selected_root_objects = get_selected_root_objects(context)
    if not selected_root_objects:
        raise ContactPipelineError("No selected root objects found")

    selected_bond_roots = [
        root_obj
        for root_obj in selected_root_objects
        if str(getattr(root_obj, "name", "") or "").strip().startswith("bond_independent_")
    ]
    if not selected_bond_roots:
        raise ContactPipelineError("No selected independent bonds found")

    selected_fragment_roots = [
        root_obj
        for root_obj in selected_root_objects
        if not str(getattr(root_obj, "name", "") or "").strip().startswith("bond_independent_")
    ]
    if not selected_fragment_roots:
        raise ContactPipelineError("No selected fragment roots found for hole matching")

    contact_graph = _build_contact_graph_from_root_objects(selected_fragment_roots)
    independent_bond_names = get_independent_bond_root_names_from_root_objects(selected_bond_roots)

    if not independent_bond_names:
        print("MOLMOL_TRY_FIX_BONDS_NO_SELECTED_INDEPENDENT_BONDS")
        return []

    fixed_bonds = try_to_fix_independent_bonds(
        contact_graph=contact_graph,
        independent_bond_names=independent_bond_names,
    )
    bpy.context.view_layer.update()
    return fixed_bonds




def _find_latest_child_build_collection(build_root_name: str = "MolMol_Built"):
    build_root = bpy.data.collections.get(str(build_root_name or "").strip())
    if build_root is None:
        return None

    child_collections = list(build_root.children)
    if not child_collections:
        return None

    child_collections.sort(key=lambda item: str(getattr(item, "name", "") or ""))
    return child_collections[-1]


def _iter_collection_objects_recursive(target_collection):
    if target_collection is None:
        return

    for obj in target_collection.objects:
        yield obj

    for child_collection in target_collection.children:
        yield from _iter_collection_objects_recursive(child_collection)


def _build_contact_graph_from_existing_collection(target_collection) -> contact_assembly_graph:
    graph = contact_assembly_graph()

    if target_collection is None:
        return graph

    seen_root_names: Set[str] = set()

    for obj in _iter_collection_objects_recursive(target_collection):
        if getattr(obj, "parent", None) is not None:
            continue

        root_name = str(getattr(obj, "name", "") or "").strip()
        if not root_name:
            continue
        if root_name in seen_root_names:
            continue

        seen_root_names.add(root_name)
        connection_markers = extract_connection_markers_from_hierarchy(obj)

        graph.add_fragment(
            placed_fragment(
                node_id=root_name,
                blender_obj_name=root_name,
                fragment_id=root_name,
                node_kind="fragment",
                connection_markers=connection_markers,
            )
        )

    return graph


def get_independent_bond_root_names_from_collection(target_collection) -> List[str]:
    bond_names: List[str] = []
    seen_bond_names: Set[str] = set()

    if target_collection is None:
        return bond_names

    for obj in _iter_collection_objects_recursive(target_collection):
        if getattr(obj, "parent", None) is not None:
            continue

        root_name = str(getattr(obj, "name", "") or "").strip()
        if not root_name:
            continue
        if not root_name.startswith("bond_independent_"):
            continue
        if root_name in seen_bond_names:
            continue

        seen_bond_names.add(root_name)
        bond_names.append(root_name)

    bond_names.sort()
    return bond_names


def try_to_fix_independent_bonds_in_collection(
    collection_name_prefix: str = "MolMol_Built",
) -> List[str]:
    """
    Convenience entry point for a Blender operator:
    find the latest build collection inside the built root collection,
    reconstruct a lightweight contact graph from scene objects, then try to
    rigidly reposition independent bonds.
    """
    target_collection = _find_latest_child_build_collection(
        build_root_name=collection_name_prefix,
    )
    if target_collection is None:
        raise ContactPipelineError(
            f"No child build collection found inside '{collection_name_prefix}'"
        )

    contact_graph = _build_contact_graph_from_existing_collection(target_collection)
    independent_bond_names = get_independent_bond_root_names_from_collection(target_collection)

    if not independent_bond_names:
        print("MOLMOL_TRY_FIX_BONDS_NO_INDEPENDENT_BONDS", target_collection.name)
        return []

    fixed_bonds = try_to_fix_independent_bonds(
        contact_graph=contact_graph,
        independent_bond_names=independent_bond_names,
    )
    bpy.context.view_layer.update()
    return fixed_bonds




def _find_child_object_by_name(root_obj, target_name: str):
    if root_obj is None:
        return None

    normalized_target_name = str(target_name or "").strip()
    if not normalized_target_name:
        return None

    for child_obj in _iter_object_hierarchy(root_obj):
        if child_obj == root_obj:
            continue

        raw_child_name = getattr(child_obj, "name", "")
        child_name = str(raw_child_name or "").strip()
        if not child_name:
            continue

        if child_name == normalized_target_name:
            return child_obj

        child_name_parts = child_name.split("__")
        if child_name_parts and child_name_parts[-1] == normalized_target_name:
            return child_obj

    return None



def _normalize_marker_token(value: str) -> str:
    value = str(value or "").strip().lower()
    if not value:
        return ""
    return value.replace("__", "_").replace("-", "_")


def _detect_marker_kind_from_name(object_name: str) -> Optional[str]:
    """
    Detect marker kind from object name.

    Canonical names:
        mk_atom, mk_hole, mk_peg, mk_axis

    Legacy names still accepted as fallback:
        mkatom, mkhole, mkpeg
    """
    normalized_name = _normalize_marker_token(object_name)

    if "mk_atom" in normalized_name or "mkatom" in normalized_name:
        return "atom"
    if "mk_hole" in normalized_name or "mkhole" in normalized_name:
        return "hole"
    if "mk_peg" in normalized_name or "mkpeg" in normalized_name:
        return "peg"
    if "mk_axis" in normalized_name or "mkaxis" in normalized_name:
        return "axis"

    return None

def _canonical_slot_key(object_name: str, marker_kind: str) -> str:
    """
    Produce a stable slot key from a Blender empty name.
    Keep the original object name semantics but normalize separators.
    """
    normalized_name = _normalize_marker_token(object_name)

    replacements = {
        "mk_atom": "__atom__",
        "mkatom": "__atom__",
        "mk_hole": "__hole__",
        "mkhole": "__hole__",
        "mk_peg": "__peg__",
        "mkpeg": "__peg__",
        "mk_axis": "__axis__",
        "mkaxis": "__axis__",
    }

    slot_key = normalized_name
    for source, target in replacements.items():
        slot_key = slot_key.replace(source, target)

    slot_key = slot_key.strip("_")
    while "___" in slot_key:
        slot_key = slot_key.replace("___", "__")

    return f"{marker_kind}:{slot_key}"


def apply_initial_node_layout(
    planner_graph,
    bindings: Iterable[planner_node_binding],
    mol,
    extra_offset_factor: float = 0.10,
) -> List[str]:
    nodes = _extract_graph_nodes(planner_graph)
    binding_map = {binding.node_id: binding for binding in bindings}
    moved_objects: List[str] = []
    scale = compute_global_layout_scale(mol, extra_offset_factor=extra_offset_factor)

    for node_id, planner_node in nodes.items():
        binding = binding_map.get(node_id)
        if binding is None:
            continue

        blender_obj = bpy.data.objects.get(binding.blender_obj_name)
        if blender_obj is None:
            raise contact_pipeline_error(
                f"Object '{binding.blender_obj_name}' not found for node '{node_id}'"
            )

        node_kind = _planner_node_kind(planner_node)
        fragment_id = _planner_node_fragment_id(planner_node)

        if node_kind == "fragment":
            if fragment_id == "bond_single":
                neighbors = _neighbor_node_ids(planner_graph, node_id)
                if len(neighbors) != 2:
                    continue

                neighbor_positions = []
                for neighbor_node_id in neighbors:
                    neighbor_node = nodes.get(neighbor_node_id)
                    if neighbor_node is None:
                        continue
                    neighbor_position = _representative_rdkit_position(mol, neighbor_node, scale)
                    if neighbor_position is None:
                        continue
                    neighbor_positions.append(np.asarray(neighbor_position, dtype=float))

                if len(neighbor_positions) != 2:
                    continue

                midpoint = 0.5 * (neighbor_positions[0] + neighbor_positions[1])
                bond_direction = _safe_normalize(neighbor_positions[1] - neighbor_positions[0])
                rotation = _rotation_matrix_between_vectors(
                    np.array([1.0, 0.0, 0.0], dtype=float),
                    bond_direction,
                )
                transform = make_transform_matrix(rotation, midpoint)
                blender_obj.matrix_world = matrix4x4_to_blender(transform)
                moved_objects.append(blender_obj.name)
                continue

            fragment_instance = _get_attr(planner_node, "fragment_instance", None)
            if fragment_instance is None:
                continue

            source_atom_indices = list(_get_attr(fragment_instance, "source_atom_indices", []) or [])

            if fragment_id == "ch3" and len(source_atom_indices) == 1:
                _place_single_atom_fragment_with_axis(
                    blender_obj=blender_obj,
                    planner_graph=planner_graph,
                    planner_node=planner_node,
                    mol=mol,
                    scale=scale,
                    axis_marker_name="ch3_mk_axis",
                )
                moved_objects.append(blender_obj.name)
                continue

            transform = np.asarray(
                fit_fragment_pose(mol, fragment_instance, scale=scale),
                dtype=float,
            )
            blender_obj.matrix_world = matrix4x4_to_blender(transform)
            moved_objects.append(blender_obj.name)
            continue

        if node_kind == "atom":
            linker_atom = _get_attr(planner_node, "linker_atom", None)
            if linker_atom is None:
                continue

            source_atom_index = _get_attr(linker_atom, "source_atom_index", None)
            if source_atom_index is None:
                continue

            atom_position = (
                np.asarray(rdkit_atom_position(mol, int(source_atom_index)), dtype=float) * scale
            )
            blender_obj.matrix_world = matrix4x4_to_blender(
                make_transform_matrix(np.eye(3, dtype=float), atom_position)
            )
            moved_objects.append(blender_obj.name)

    bpy.context.view_layer.update()
    return moved_objects

def extract_marker_slots_from_hierarchy(root_obj) -> Dict[str, peg_hole_slot]:
    """
    Read EMPTY markers from a fragment/template hierarchy and extract
    canonical marker slots from their world-space positions.

    Supported marker kinds:
    - atom
    - hole
    - peg
    - axis

    Legacy names are accepted only as fallback.
    """
    slots: Dict[str, peg_hole_slot] = {}

    if root_obj is None:
        return slots

    parent_root_name = str(getattr(root_obj, "name", "") or "").strip()

    for obj in _iter_object_hierarchy(root_obj):
        if obj is None or obj == root_obj:
            continue

        if getattr(obj, "type", None) != "EMPTY":
            continue

        object_name = str(getattr(obj, "name", "") or "").strip()
        if not object_name:
            continue

        marker_kind = _detect_marker_kind_from_name(object_name)
        if marker_kind is None:
            continue

        world_location = obj.matrix_world.translation
        location = np.array(
            [world_location.x, world_location.y, world_location.z],
            dtype=float,
        )

        slot_key = _canonical_slot_key(object_name, marker_kind)

        slots[slot_key] = peg_hole_slot(
            slot_key=slot_key,
            kind=marker_kind,
            source_object_name=object_name,
            parent_root_name=parent_root_name,
            location=location,
        )

    return slots

def get_atom_live_world_pos(atom_idx, atom_local_positions, planner_graph, contact_graph) -> "Vector":
    """
    Return the current world-space position of an atom in Blender,
    using the live object matrix and the precomputed local atom offset.
    """
    import bpy
    from mathutils import Vector

    local_position = atom_local_positions.get(int(atom_idx))
    if local_position is None:
        raise ContactPipelineError(f"Missing local position for atom index {atom_idx}")

    owner_node_id = None
    for node_id, node in planner_graph.nodes.items():
        if (
            getattr(node, "kind", None) == "linker_atom"
            and getattr(node, "linker_atom", None) is not None
            and int(node.linker_atom.source_atom_index) == int(atom_idx)
        ):
            owner_node_id = node_id
            break

        if (
            getattr(node, "kind", None) == "fragment"
            and getattr(node, "fragment_instance", None) is not None
            and int(atom_idx) in list(node.fragment_instance.source_atom_indices)
        ):
            owner_node_id = node_id
            break

    if owner_node_id is None:
        raise ContactPipelineError(f"No planner node owns atom index {atom_idx}")

    fragment_data = contact_graph.get_fragment(owner_node_id)
    if fragment_data is None:
        raise ContactPipelineError(f"No contact fragment found for node '{owner_node_id}'")

    obj = bpy.data.objects.get(fragment_data.blender_obj_name)
    if obj is None:
        raise ContactPipelineError(
            f"Blender object '{fragment_data.blender_obj_name}' not found for node '{owner_node_id}'"
        )

    return obj.matrix_world @ local_position




def build_from_file_with_contact_pipeline(
    context,
    structure_path: str,
    library_path: str,
    collection_name: str = "MolMol_ContactBuild",
    extra_offset_factor: float = 0.10,
) -> build_execution_result:
    structure_path = str(structure_path or "").strip()
    library_path = str(library_path or "").strip()

    if not structure_path:
        raise ContactPipelineError("Structure file path is empty")
    if not library_path:
        raise ContactPipelineError("Library file path is empty")

    mol = load_molecule(structure_path)

    ctx_scene = getattr(context, "scene", None)
    molmol_settings = getattr(ctx_scene, "molmolsettings", None)

    add_hydrogens = False
    if molmol_settings is not None:
        add_hydrogens = bool(getattr(molmol_settings, "add_hydrogens", False))

    if add_hydrogens:
        mol = Chem.AddHs(mol, addCoords=True)
        print("MOLMOL_ADD_HYDROGENS", mol.GetNumAtoms())

    planner_graph = build_assembly_graph(mol)
    scale = compute_global_layout_scale(mol, extra_offset_factor=extra_offset_factor)

    center_molecule = False
    if ctx_scene and hasattr(ctx_scene, "molmol_settings"):
        center_molecule = bool(ctx_scene.molmol_settings.center_rdkit)

    if center_molecule:
        atom_positions = [
            np.asarray(rdkit_atom_position(mol, atom.GetIdx()), dtype=float)
            for atom in mol.GetAtoms()
        ]
        if atom_positions:
            center_offset = np.mean(np.asarray(atom_positions, dtype=float), axis=0)
        else:
            center_offset = np.array([0.0, 0.0, 0.0], dtype=float)
    else:
        center_offset = np.array([0.0, 0.0, 0.0], dtype=float)

    bindings, instantiated_objects, built_collection_name = instantiate_planner_graph(
        planner_graph=planner_graph,
        collection_name=collection_name,
        library_path=library_path,
    )

    apply_initial_node_layout(
        planner_graph=planner_graph,
        bindings=bindings,
        mol=mol,
        extra_offset_factor=extra_offset_factor,
    )

    if center_molecule:
        offset_value = Vector(center_offset.tolist()) * float(scale)
        for obj_name in instantiated_objects:
            obj = bpy.data.objects.get(obj_name)
            if obj is not None and obj.parent is None:
                obj.location -= offset_value

    contact_graph, unresolved_connections = build_contact_assembly_graph(
        planner_graph,
        bindings,
    )
    bpy.context.view_layer.update()

    binding_map = {binding.node_id: binding for binding in bindings}

    atom_local_positions = {}
    for node_id, node in planner_graph.nodes.items():
        frag_data = contact_graph.get_fragment(node_id)
        if frag_data is None:
            continue

        obj = bpy.data.objects.get(frag_data.blender_obj_name)
        if obj is None:
            continue

        atom_indices = []
        if getattr(node, "kind", None) == "linker_atom" and getattr(node, "linker_atom", None) is not None:
            atom_indices = [node.linker_atom.source_atom_index]
        elif getattr(node, "kind", None) == "fragment" and getattr(node, "fragment_instance", None) is not None:
            atom_indices = list(node.fragment_instance.source_atom_indices)

        obj_inv_matrix = obj.matrix_world.inverted()
        for atom_index in atom_indices:
            raw_pos = np.asarray(rdkit_atom_position(mol, int(atom_index)), dtype=float)
            raw_pos = raw_pos - center_offset
            world_pos_initial = Vector(raw_pos.tolist()) * float(scale)
            atom_local_positions[int(atom_index)] = obj_inv_matrix @ world_pos_initial

    match_result = resolve_all_connections(
        contact_graph,
        unresolved_connections,
    )

    moved_bonds = _generate_independent_single_bonds(
        mol=mol,
        planner_graph=planner_graph,
        contact_graph=contact_graph,
        binding_map=binding_map,
        collection_name=collection_name,
        atom_local_positions=atom_local_positions,
        library_path=library_path,
    )

    bpy.context.view_layer.update()

    all_objects = list(instantiated_objects) + moved_bonds
    return build_execution_result(
        planner_graph=planner_graph,
        contact_graph=contact_graph,
        match_result=match_result,
        built_collection_name=built_collection_name,
        instantiated_objects=list(set(all_objects)),
    )

def append_fragment_collection_from_library(
    library_path: str,
    collection_name: str,
) -> Optional[bpy.types.Collection]:
    """
    Import an entire named collection from a .blend library into bpy.data.

    This preserves the full parent/child hierarchy of all objects inside.
    Returns the imported bpy.data.collections entry, or None if not found.
    The collection is linked to the scene only temporarily if needed for evaluation.
    """
    library_path = _normalize_library_path(library_path)
    collection_name = str(collection_name or "").strip()

    if not library_path or not collection_name:
        return None

    existing = bpy.data.collections.get(collection_name)
    if existing is not None:
        return existing

    try:
        with bpy.data.libraries.load(library_path, link=False) as (data_from, data_to):
            if collection_name not in data_from.collections:
                print(
                    "MOLMOL_LIBRARY_MISSING_COLLECTION",
                    repr(collection_name),
                    f"available={list(data_from.collections)}",
                )
                return None
            data_to.collections = [collection_name]
    except Exception as exc:
        raise ContactPipelineError(
            f"Cannot load collection '{collection_name}' from library '{library_path}': {exc}"
        )

    imported = bpy.data.collections.get(collection_name)
    if imported is None:
        raise ContactPipelineError(
            f"Collection '{collection_name}' not found in bpy.data after import"
        )

    scene_collection = bpy.context.scene.collection
    was_linked_to_scene = any(child == imported for child in scene_collection.children)

    if not was_linked_to_scene:
        scene_collection.children.link(imported)

    bpy.context.view_layer.update()

    if not was_linked_to_scene:
        scene_collection.children.unlink(imported)

    return imported

def _generate_independent_single_bonds(
    mol,
    planner_graph,
    contact_graph,
    binding_map,
    collection_name,
    atom_local_positions,
    library_path: str = "",
) -> List[str]:
    """
    Generate independent bond hierarchies only for single bonds that are not
    already embedded inside fragment templates.

    The bond template is instantiated as a full root-with-children hierarchy so
    all markers (for example mk_peg empties) follow the same world transform as
    the bond mesh.
    """
    bond_order = 1
    primary_template_name = get_bond_template_name(bond_order)
    primary_collection_name = get_bond_template_collection_name(bond_order)

    template_name_candidates = [primary_template_name]
    legacy_name_candidates = [
        "bond_single_medium",
        "bond_single",
    ]

    for candidate_name in legacy_name_candidates:
        if candidate_name not in template_name_candidates:
            template_name_candidates.append(candidate_name)

    template_source_obj = None

    for candidate_name in template_name_candidates:
        template_source_obj = bpy.data.objects.get(candidate_name)
        if template_source_obj is not None:
            print("MOLMOL_BOND_TEMPLATE_FOUND", candidate_name, "FROM_DATA_OBJECTS")
            break

    if template_source_obj is None and library_path:
        imported_collection = append_fragment_collection_from_library(
            library_path=library_path,
            collection_name=primary_collection_name,
        )
        bpy.context.view_layer.update()

        if imported_collection is not None:
            for obj in imported_collection.all_objects:
                if obj.name in template_name_candidates:
                    template_source_obj = obj
                    print("MOLMOL_BOND_TEMPLATE_FOUND", obj.name, "FROM_IMPORTED_COLLECTION")
                    break

    if template_source_obj is None:
        for candidate_name in template_name_candidates:
            try:
                template_source_obj = resolve_template_source_object(
                    source_name=candidate_name,
                    library_path=library_path,
                )
            except Exception:
                template_source_obj = None

            if template_source_obj is not None:
                print("MOLMOL_BOND_TEMPLATE_FOUND", candidate_name, "FROM_TEMPLATE_LOOKUP")
                break

    if template_source_obj is None:
        print(
            "MOLMOL_WARN_MISSING_BOND_TEMPLATE",
            primary_template_name,
            primary_collection_name,
        )
        return []

    target_collection = bpy.data.collections.get(collection_name)
    if target_collection is None:
        raise ContactPipelineError(
            f"Collection '{collection_name}' not found for independent bond generation"
        )

    covered_template_bond_pairs = _build_internal_template_bond_pairs(
        mol=mol,
        planner_graph=planner_graph,
        binding_map=binding_map,
    )

    print("MOLMOL_COVERED_TEMPLATE_BOND_PAIRS", sorted(covered_template_bond_pairs))

    moved_bonds: List[str] = []

    for bond_index, bond in enumerate(mol.GetBonds()):
        current_bond_order = int(round(float(bond.GetBondTypeAsDouble())))
        current_bond_order = max(1, min(3, current_bond_order))
        if current_bond_order != 1:
            continue

        idx_a = int(bond.GetBeginAtomIdx())
        idx_b = int(bond.GetEndAtomIdx())
        pair = tuple(sorted((idx_a, idx_b)))

        if pair in covered_template_bond_pairs:
            print("SKIP_TEMPLATE_INTERNAL_BOND", bond_index, pair, "order", current_bond_order)
            continue

        v_a = get_atom_live_world_pos(
            atom_idx=idx_a,
            atom_local_positions=atom_local_positions,
            planner_graph=planner_graph,
            contact_graph=contact_graph,
        )
        v_b = get_atom_live_world_pos(
            atom_idx=idx_b,
            atom_local_positions=atom_local_positions,
            planner_graph=planner_graph,
            contact_graph=contact_graph,
        )

        direction = v_b - v_a
        if direction.length < 1e-6:
            print("SKIP_ZERO_LENGTH_BOND", bond_index, idx_a, idx_b)
            continue

        new_bond_name = f"bond_independent_{bond_index:03d}_{idx_a}_{idx_b}"

        existing = bpy.data.objects.get(new_bond_name)
        if existing is not None:
            bpy.data.objects.remove(existing, do_unlink=True)

        duplicate_result = duplicate_root_with_children(
            source_obj=template_source_obj,
            target_collection=target_collection,
            new_name=new_bond_name,
        )

        if duplicate_result is None:
            raise ContactPipelineError(
                f"duplicate_root_with_children returned None for independent bond '{new_bond_name}'"
            )

        new_bond_root = duplicate_result[0] if isinstance(duplicate_result, tuple) else duplicate_result
        if new_bond_root is None:
            raise ContactPipelineError(
                f"duplicate_root_with_children returned no root object for '{new_bond_name}'"
            )

        midpoint = (v_a + v_b) / 2.0
        rotation = Vector((0.0, 0.0, 1.0)).rotation_difference(direction.normalized())

        new_bond_root.rotation_mode = "QUATERNION"
        new_bond_root.location = midpoint
        new_bond_root.rotation_quaternion = rotation

        moved_bonds.append(new_bond_root.name)

        print(
            "MOLMOL_INDEPENDENT_BOND_CREATED",
            new_bond_root.name,
            "atoms",
            idx_a,
            idx_b,
        )
        debug_print_connection_markers(new_bond_root)

    bpy.context.view_layer.update()
    return moved_bonds


def _build_template_local_to_rdkit_atom_map(
    mol,
    fragment_template,
    source_atom_indices: List[int],
) -> Dict[int, int]:
    """
    Match template local atom indices to RDKit atom indices by global geometric overlap.
    No ordering assumption is used.
    """
    from scipy.optimize import linear_sum_assignment

    template_positions = np.asarray(fragment_template.atom_local_positions, dtype=float)

    rdkit_positions = np.asarray(
        [
            np.asarray(rdkit_atom_position(mol, int(atom_index)), dtype=float)
            for atom_index in source_atom_indices
        ],
        dtype=float,
    )

    if template_positions.shape[0] != rdkit_positions.shape[0]:
        raise ContactPipelineError(
            f"Template atom count ({template_positions.shape[0]}) does not match "
            f"fragment source atom count ({rdkit_positions.shape[0]}) "
            f"for fragment '{fragment_template.id}'"
        )

    template_centroid = template_positions.mean(axis=0)
    rdkit_centroid = rdkit_positions.mean(axis=0)

    template_centered = template_positions - template_centroid
    rdkit_centered = rdkit_positions - rdkit_centroid

    distance_matrix = np.linalg.norm(
        template_centered[:, None, :] - rdkit_centered[None, :, :],
        axis=2,
    )

    row_ind, col_ind = linear_sum_assignment(distance_matrix)

    mapping: Dict[int, int] = {}
    for template_local_index, rdkit_column in zip(row_ind.tolist(), col_ind.tolist()):
        mapping[int(template_local_index)] = int(source_atom_indices[int(rdkit_column)])

    return mapping