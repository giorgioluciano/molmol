from __future__ import annotations

from typing import Dict, List

import bpy
import mathutils
from mathutils import Matrix, Vector
from rdkit import Chem

from .fragment_recognizer import load_molecule
from .fragment_library import FRAGMENT_LIBRARY
from .geometry_utils import rdkit_atom_position


ATOM_COLORS = {
    "H": (1.0, 1.0, 1.0, 1.0),
    "C": (0.15, 0.15, 0.15, 1.0),
    "N": (0.2, 0.3, 0.9, 1.0),
    "O": (0.9, 0.2, 0.2, 1.0),
    "Cl": (0.2, 0.8, 0.2, 1.0),
    "Br": (0.6, 0.25, 0.1, 1.0),
}

TEMPLATE_COLLECTION_NAME = "MolMol_Templates"
BUILD_ROOT_COLLECTION_NAME = "MolMol_Built"
FRAGMENT_BUILD_COLLECTION_NAME = "MolMol_Built_From_File"

HALF_BOND_TEMPLATE = "half_single_bond"
HALF_BOND_HOLE_ATOM = "half_bond_hole1"
HALF_BOND_HOLE_OUTER = "half_bond_hole2"




# ===========================================================
# BASIC UTILS
# ===========================================================

def update_view() -> None:
    bpy.context.view_layer.update()


def get_or_create_collection(context, name: str, parent=None):
    coll = bpy.data.collections.get(name)
    if coll is None:
        coll = bpy.data.collections.new(name)

    if parent is None:
        parent = context.scene.collection

    if not any(child == coll for child in parent.children):
        parent.children.link(coll)

    return coll


def ensure_child_collection(parent_collection, child_name: str):
    child = bpy.data.collections.get(child_name)
    if child is None:
        child = bpy.data.collections.new(child_name)

    if not any(c == child for c in parent_collection.children):
        parent_collection.children.link(child)

    return child


def ensure_fragment_collection(name: str, parent=None):
    coll = bpy.data.collections.get(name)
    if coll is None:
        coll = bpy.data.collections.new(name)
        if parent is None:
            parent = bpy.context.scene.collection
        parent.children.link(coll)
    return coll


def clear_collection_objects_only(collection) -> None:
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)


def clear_fragment_collection_objects(collection) -> None:
    clear_collection_objects_only(collection)


def get_or_create_material(name: str, color):
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
        mat.use_nodes = True

    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf is not None:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.4

    return mat


def apply_fallback_material(obj, element: str) -> None:
    if obj.type != "MESH":
        return

    color = ATOM_COLORS.get(element, (0.8, 0.8, 0.8, 1.0))
    mat = get_or_create_material(f"MolMol_{element}", color)

    if obj.data.materials:
        obj.data.materials.clear()
    obj.data.materials.append(mat)


def tint_template_hierarchy_if_possible(root_obj, child_map, element: str) -> None:
    if root_obj.type == "MESH":
        apply_fallback_material(root_obj, element)

    for child in child_map.values():
        if child.type == "MESH":
            apply_fallback_material(child, element)


def clone_simple_object(source_obj, collection, new_name: str, location=(0, 0, 0)):
    obj = source_obj.copy()
    if source_obj.data is not None:
        obj.data = source_obj.data.copy()

    obj.animation_data_clear()
    obj.parent = None
    obj.constraints.clear()
    obj.location = Vector(location)
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
    obj.scale = source_obj.scale.copy()

    collection.objects.link(obj)

    for existing_collection in list(obj.users_collection):
        if existing_collection is not collection:
            existing_collection.objects.unlink(obj)

    obj.name = new_name
    update_view()
    return obj


# ===========================================================
# TEMPLATE LOOKUP
# ===========================================================

def _build_fragment_root_to_collection_map() -> Dict[str, str]:
    """
    Build a normalized mapping:
        blender_template_name -> blender_collection_name
    from fragment_library.FRAGMENT_LIBRARY.
    """
    root_to_collection: Dict[str, str] = {}

    for fragment in FRAGMENT_LIBRARY.values():
        root_name = (fragment.blender_template_name or "").strip().lower()
        collection_name = (fragment.blender_collection_name or "").strip()

        if not root_name or not collection_name:
            continue

        root_to_collection[root_name] = collection_name

    return root_to_collection


_FRAGMENT_ROOT_TO_COLLECTION: Dict[str, str] = _build_fragment_root_to_collection_map()
_FRAGMENT_ROOT_NAMES = set(_FRAGMENT_ROOT_TO_COLLECTION.keys())


def is_fragment_root_name(name: str) -> bool:
    normalized_name = (name or "").strip().lower()
    if not normalized_name:
        return False
    return normalized_name in _FRAGMENT_ROOT_NAMES


def map_fragment_root_to_collection(root_name: str) -> str:
    normalized_root_name = (root_name or "").strip().lower()
    if not normalized_root_name:
        return ""
    return _FRAGMENT_ROOT_TO_COLLECTION.get(normalized_root_name, "")


def is_fragment_root_template_name(template_name: str) -> bool:
    if not template_name:
        return False

    normalized_name = str(template_name).strip().lower()
    return normalized_name in _FRAGMENT_ROOT_NAMES


def find_fragment_root_source_object(template_name):
    if not template_name:
        return None

    normalized_name = str(template_name).strip()
    exact_object = bpy.data.objects.get(normalized_name)

    if exact_object is not None:
        children_recursive = getattr(exact_object, "children_recursive", exact_object.children)
        if exact_object.parent is None and len(children_recursive) > 0:
            print("MOLMOL_TEMPLATE_ROOT_RESOLVED", normalized_name, "FROM_DATA_OBJECTS_EXACT")
            return exact_object

    expected_collection_name = map_fragment_root_to_collection(normalized_name)
    if expected_collection_name:
        collection = bpy.data.collections.get(expected_collection_name)
        if collection is not None:
            for obj in collection.objects:
                if obj.name != normalized_name:
                    continue
                children_recursive = getattr(obj, "children_recursive", obj.children)
                if obj.parent is None and len(children_recursive) > 0:
                    print(
                        "MOLMOL_TEMPLATE_ROOT_RESOLVED",
                        normalized_name,
                        "FROM_COLLECTION",
                        expected_collection_name,
                    )
                    return obj

    for collection in bpy.data.collections:
        if not collection.name.startswith("frag_"):
            continue

        for obj in collection.objects:
            if obj.name != normalized_name:
                continue

            children_recursive = getattr(obj, "children_recursive", obj.children)
            if obj.parent is None and len(children_recursive) > 0:
                print(
                    "MOLMOL_TEMPLATE_ROOT_RESOLVED",
                    normalized_name,
                    "FROM_COLLECTION",
                    collection.name,
                )
                return obj

    for obj in bpy.data.objects:
        if obj.name != normalized_name:
            continue

        children_recursive = getattr(obj, "children_recursive", obj.children)
        if obj.parent is None and len(children_recursive) > 0:
            print("MOLMOL_TEMPLATE_ROOT_RESOLVED", normalized_name, "FROM_DATA_OBJECTS_SCAN")
            return obj

    if exact_object is not None:
        print("MOLMOL_TEMPLATE_ROOT_FALLBACK", normalized_name, "NO_CHILDREN")
        return exact_object

    return None


def find_template_object_in_template_collection(template_name):
    template_collection = bpy.data.collections.get(TEMPLATE_COLLECTION_NAME)
    if template_collection is None:
        return None

    for obj in template_collection.objects:
        if obj.name == template_name:
            return obj

    return None


def get_template_source_object(template_name):
    normalized_name = str(template_name).strip() if template_name else ""
    if not normalized_name:
        return None

    if is_fragment_root_template_name(normalized_name):
        source_object = find_fragment_root_source_object(normalized_name)
        if source_object is not None:
            print("MOLMOL_TEMPLATE_LOOKUP", normalized_name, "FRAGMENT_ROOT")
            return source_object

    source_object = find_template_object_in_template_collection(normalized_name)
    if source_object is not None:
        print("MOLMOL_TEMPLATE_LOOKUP", normalized_name, "TEMPLATE_COLLECTION")
        return source_object

    source_object = bpy.data.objects.get(normalized_name)
    if source_object is not None:
        print("MOLMOL_TEMPLATE_LOOKUP", normalized_name, "DATA_OBJECTS")
        return source_object

    print("MOLMOL_TEMPLATE_LOOKUP_FAILED", normalized_name)
    return None


# ===========================================================
# HIERARCHY DUPLICATION
# ===========================================================

def duplicate_root_with_children(source_obj, target_collection, new_name: str):
    """
    Deep-copy root object and all recursive children into target_collection.
    Mesh data is copied to preserve vertex groups and per-instance editability.

    Returns
    -------
    tuple
        (root_copy, child_map) where child_map maps original child names to copied children.
    """
    root_copy = source_obj.copy()
    if source_obj.data is not None:
        root_copy.data = source_obj.data.copy()

    root_copy.animation_data_clear()
    root_copy.parent = None
    root_copy.constraints.clear()
    root_copy.name = new_name
    root_copy.rotation_mode = "QUATERNION"
    root_copy.rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
    target_collection.objects.link(root_copy)

    child_map = {}

    def _recurse(src_parent, dst_parent):
        for child in src_parent.children:
            child_copy = child.copy()
            if child.data is not None:
                child_copy.data = child.data.copy()

            child_copy.animation_data_clear()
            child_copy.constraints.clear()
            child_copy.name = f"{new_name}__{child.name}"
            target_collection.objects.link(child_copy)

            original_world_matrix = child.matrix_world.copy()
            child_copy.parent = dst_parent
            child_copy.matrix_parent_inverse = dst_parent.matrix_world.inverted()
            child_copy.matrix_world = original_world_matrix

            child_map[child.name] = child_copy
            _recurse(child, child_copy)

    _recurse(source_obj, root_copy)
    update_view()
    return root_copy, child_map


# ===========================================================
# ATOM OBJECTS
# ===========================================================

def create_atom_sphere(collection, atom, location, radius: float):
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=float(radius),
        location=location,
        segments=32,
        ring_count=16,
    )
    obj = bpy.context.active_object
    obj.name = f"Atom_{atom['id']}_{atom['element']}"

    if obj.users_collection:
        for existing_collection in list(obj.users_collection):
            existing_collection.objects.unlink(obj)

    collection.objects.link(obj)
    apply_fallback_material(obj, atom["element"])
    return obj


def create_text_label(collection, text: str, location, size: float = 0.35):
    curve = bpy.data.curves.new(type="FONT", name=f"Label_{text}")
    curve.body = text
    curve.size = float(size)

    obj = bpy.data.objects.new(f"Label_{text}", curve)
    obj.location = Vector(location) + Vector((0.0, 0.0, 0.7))
    collection.objects.link(obj)
    return obj


def instantiate_atom_object(
    context,
    collection,
    atom,
    location,
    atom_radius: float = 0.35,
    add_label: bool = True,
):
    del context

    template_name = atom.get("template")
    object_name = f"Atom_{atom['id']}_{atom['element']}"

    source_obj = None
    if template_name:
        source_obj = get_template_source_object(template_name)

    if source_obj is not None:
        obj, child_map = duplicate_root_with_children(
            source_obj=source_obj,
            target_collection=collection,
            new_name=object_name,
        )
        obj.location = Vector(location)
        obj.rotation_mode = "QUATERNION"
        obj.rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
        tint_template_hierarchy_if_possible(obj, child_map, atom["element"])
    else:
        obj = create_atom_sphere(
            collection=collection,
            atom=atom,
            location=location,
            radius=float(atom_radius),
        )
        child_map = {}

    if add_label:
        create_text_label(collection, atom["element"], location)

    return obj, child_map


# ===========================================================
# SLOT / HOLE ACCESS
# ===========================================================

def find_slot_child(atom_obj, slot_index: int):
    token = f"_hole{slot_index + 1}".lower()

    for child in atom_obj.children:
        child_name = child.name.lower()
        if child_name.endswith(token):
            return child

    for child in atom_obj.children:
        child_name = child.name.lower()
        if token in child_name:
            return child

    return None


def get_slot_world_pose(atom_obj, slot_index: int):
    slot_child = find_slot_child(atom_obj, slot_index)
    if slot_child is not None:
        pos = slot_child.matrix_world.translation.copy()
        mx = slot_child.matrix_world.to_3x3()
        xw = (mx @ Vector((1.0, 0.0, 0.0))).normalized()
        yw = (mx @ Vector((0.0, 1.0, 0.0))).normalized()
        zw = (mx @ Vector((0.0, 0.0, 1.0))).normalized()

        print(
            f"[MolMol][HOLE] atom={atom_obj.name} slot={slot_index} "
            f"pos={tuple(round(v, 4) for v in pos)} "
            f"X={tuple(round(v, 4) for v in xw)} "
            f"Y={tuple(round(v, 4) for v in yw)} "
            f"Z={tuple(round(v, 4) for v in zw)}"
        )

        local_pos = slot_child.matrix_local.translation.copy()
        if local_pos.length > 1e-8:
            direction = (atom_obj.matrix_world.to_3x3() @ local_pos.normalized()).normalized()
        else:
            direction = Vector((1.0, 0.0, 0.0))

        return {
            "position": pos,
            "direction": direction,
            "slot_child": slot_child,
        }

    return {
        "position": atom_obj.matrix_world.translation.copy(),
        "direction": Vector((1.0, 0.0, 0.0)),
        "slot_child": None,
    }


# ===========================================================
# SIMPLE ATOM ORIENTATION
# ===========================================================

def build_heavy_neighbor_map(plan):
    neighbor_map = {}

    for connection in plan["connection_instructions"]:
        atom_a = connection["atom_a"]
        atom_b = connection["atom_b"]

        neighbor_map.setdefault(atom_a, []).append(
            {"other_atom": atom_b, "slot": connection["atom_a_slot"]}
        )
        neighbor_map.setdefault(atom_b, []).append(
            {"other_atom": atom_a, "slot": connection["atom_b_slot"]}
        )

    return neighbor_map


# ===========================================================
# RDKIT BALLS AND STICKS
# ===========================================================

def _rdkit_visible_atom_indices(mol, show_hydrogens: bool) -> List[int]:
    indices: List[int] = []
    for atom in mol.GetAtoms():
        if not show_hydrogens and atom.GetSymbol() == "H":
            continue
        indices.append(int(atom.GetIdx()))
    return indices


def _rdkit_center_offset(mol, visible_atom_indices: List[int]) -> Vector:
    if not visible_atom_indices:
        return Vector((0.0, 0.0, 0.0))

    accumulator = Vector((0.0, 0.0, 0.0))
    for atom_idx in visible_atom_indices:
        position = rdkit_atom_position(mol, atom_idx)
        accumulator += Vector((float(position[0]), float(position[1]), float(position[2])))

    return accumulator / len(visible_atom_indices)


def _rdkit_atom_location(mol, atom_idx: int, scale: float, center_offset: Vector) -> Vector:
    position = rdkit_atom_position(mol, atom_idx)
    vector = Vector((float(position[0]), float(position[1]), float(position[2])))
    return (vector - center_offset) * float(scale)


def create_bond_cylinder(collection, p1, p2, radius: float = 0.12, name: str = "Bond"):
    p1 = Vector(p1)
    p2 = Vector(p2)
    delta = p2 - p1
    length = delta.length

    if length <= 1e-8:
        return None

    midpoint = (p1 + p2) * 0.5
    direction = delta.normalized()

    bpy.ops.mesh.primitive_cylinder_add(
        radius=float(radius),
        depth=float(length),
        location=tuple(midpoint),
        vertices=24,
    )
    obj = bpy.context.active_object
    obj.name = name

    up = Vector((0.0, 0.0, 1.0))
    quat = up.rotation_difference(direction)
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = quat

    if obj.users_collection:
        for existing_collection in list(obj.users_collection):
            existing_collection.objects.unlink(obj)
    collection.objects.link(obj)

    mat = get_or_create_material("MolMol_Bond", (0.7, 0.7, 0.7, 1.0))
    if obj.data.materials:
        obj.data.materials.clear()
    obj.data.materials.append(mat)

    return obj


def build_rdkit_balls_and_sticks(
    context,
    structure_path: str,
    molecule_name: str = "MolMol_RDKit",
    atom_radius: float = 0.35,
    bond_radius: float = 0.12,
    coordinate_scale: float = 1.0,
    add_hydrogens: bool = False,
    show_hydrogens: bool = False,
    center_molecule: bool = True,
    add_labels: bool = False,
):
    """
    Build a debug/validation RDKit balls+sticks representation in Blender.
    This is a visual reference tool, not the main fragment assembly pipeline.
    """
    mol = load_molecule(structure_path)
    if add_hydrogens:
        mol = Chem.AddHs(mol, addCoords=True)

    root_collection = get_or_create_collection(context, BUILD_ROOT_COLLECTION_NAME)
    molecule_collection = ensure_child_collection(root_collection, molecule_name)
    clear_collection_objects_only(molecule_collection)

    visible_atom_indices = _rdkit_visible_atom_indices(
        mol,
        show_hydrogens=show_hydrogens,
    )
    center_offset = (
        _rdkit_center_offset(mol, visible_atom_indices)
        if center_molecule
        else Vector((0.0, 0.0, 0.0))
    )

    atom_locations = {}
    atom_objects = {}

    for atom in mol.GetAtoms():
        atom_idx = int(atom.GetIdx())
        element = atom.GetSymbol()

        if atom_idx not in visible_atom_indices:
            continue

        location = _rdkit_atom_location(
            mol=mol,
            atom_idx=atom_idx,
            scale=coordinate_scale,
            center_offset=center_offset,
        )
        atom_locations[atom_idx] = location

        atom_payload = {
            "id": atom_idx,
            "element": element,
            "template": None,
        }

        obj = create_atom_sphere(
            collection=molecule_collection,
            atom=atom_payload,
            location=tuple(location),
            radius=float(atom_radius),
        )
        obj.name = f"RDKitAtom_{atom_idx:03d}_{element}"
        atom_objects[atom_idx] = obj

        if add_labels:
            create_text_label(
                collection=molecule_collection,
                text=f"{atom_idx}:{element}",
                location=tuple(location),
                size=0.25,
            )

    bond_count = 0
    for bond in mol.GetBonds():
        atom_a = int(bond.GetBeginAtomIdx())
        atom_b = int(bond.GetEndAtomIdx())

        if atom_a not in atom_locations or atom_b not in atom_locations:
            continue

        create_bond_cylinder(
            collection=molecule_collection,
            p1=atom_locations[atom_a],
            p2=atom_locations[atom_b],
            radius=float(bond_radius),
            name=f"RDKitBond_{bond_count:03d}_{atom_a}_{atom_b}",
        )
        bond_count += 1

    update_view()

    return {
        "collection": molecule_collection,
        "atom_count": len(atom_objects),
        "bond_count": bond_count,
        "centered": bool(center_molecule),
        "show_hydrogens": bool(show_hydrogens),
    }


def register():
    pass


def unregister():
    pass