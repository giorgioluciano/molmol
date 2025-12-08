import bpy, os, itertools
from mathutils import Vector, Matrix, Quaternion
from math import radians
from ase.io import read as ase_read
from ase.neighborlist import NeighborList
from ase.data import covalent_radii

GEN_COLLECTIONS = ["Atom_sp3", "Atom_sp2", "Atom_sp", "Atom_bent", "Atom_sp3d2"]
HALOGENS = {"F", "Cl", "Br", "I"}

def abspath(path):
    """Get Blender-absolute or OS-absolute path."""
    return bpy.path.abspath(path)

def ensure_hidden_bucket():
    """Get or create a collection for hidden technical objects."""
    name = "_MolymodHiddenObjs"
    coll = bpy.data.collections.get(name)
    if not coll:
        coll = bpy.data.collections.new(name)
    return coll

def unlink_collection_everywhere(coll):
    """Recursively remove a collection from all possible parents in the .blend."""
    target_name = coll.name
    def rec(parent):
        for ch in list(parent.children):
            if ch.name == target_name:
                try: parent.children.unlink(ch)
                except: pass
            else:
                rec(ch)
    rec(bpy.context.scene.collection)
    for parent in list(bpy.data.collections):
        if parent.children.get(target_name) is not None:
            try: parent.children.unlink(parent.children[target_name])
            except: pass

def append_collection(lib_path, coll_name):
    """Append a collection from an external .blend library if not present."""
    lib_path = abspath(lib_path)
    colldir = os.path.join(lib_path, "Collection")
    if not os.path.isfile(lib_path):
        raise FileNotFoundError(f"Library .blend not found: {lib_path}")
    coll = bpy.data.collections.get(coll_name)
    if not coll:
        bpy.ops.wm.append(directory=colldir, filename=coll_name)
        coll = bpy.data.collections[coll_name]
    unlink_collection_everywhere(coll)
    return coll

def append_object(lib_path, obj_name):
    """Append a mesh/object from an external library into the hidden bucket."""
    lib_path = abspath(lib_path)
    objdir = os.path.join(lib_path, "Object")
    if not os.path.isfile(lib_path):
        raise FileNotFoundError(f"Library .blend not found: {lib_path}")
    if obj_name not in bpy.data.objects:
        bpy.ops.wm.append(directory=objdir, filename=obj_name)
    ob = bpy.data.objects[obj_name]
    hidden = ensure_hidden_bucket()
    for c in list(ob.users_collection):
        try: c.objects.unlink(ob)
        except: pass
    if ob.name not in hidden.objects:
        hidden.objects.link(ob)
    ob.hide_render = True
    ob.hide_select = True
    return ob

def load_hole_dirs(lib_path, coll_key):
    """Load direction vectors (as Blender Vector) for all holes in a collection."""
    holes = []
    with bpy.data.libraries.load(abspath(lib_path), link=False) as (src, dst):
        cand = [n for n in src.objects if n.startswith(f"{coll_key}_hole")]
    for h in cand:
        ob = append_object(lib_path, h)
        if ob.location.length > 1e-9:
            holes.append(ob.location.normalized())
    if bpy.context.scene and bpy.context.scene.molymod_settings.debug_mode:
        print(f"[HOLES] {coll_key}: found {len(holes)} hole vectors")
    return holes

def kabsch_rotation(from_vecs, to_vecs):
    """Calculate optimal rotation matrix (Kabsch algorithm) to align from_vecs to to_vecs."""
    import numpy as np
    A = np.array([[v.x, v.y, v.z] for v in from_vecs], dtype=float).T
    B = np.array([[v.x, v.y, v.z] for v in to_vecs], dtype=float).T
    H = A @ B.T
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return Matrix(((R[0,0], R[0,1], R[0,2]),
                   (R[1,0], R[1,1], R[1,2]),
                   (R[2,0], R[2,1], R[2,2]))).to_4x4()

def align_one_vector(src: Vector, dst: Vector):
    """Align a 'src' vector to a 'dst' vector using quaternion rotation."""
    s = src.normalized(); d = dst.normalized()
    q = s.rotation_difference(d)
    return q.to_matrix().to_4x4()

def hungarian_assign(cost):
    """Optimal assignment (Hungarian method) for hole-to-bond vector matching."""
    try:
        import numpy as np
        import scipy.optimize as spopt
        r, c = spopt.linear_sum_assignment(np.array(cost, dtype=float))
        return list(r), list(c)
    except Exception:
        m = len(cost); n = len(cost[0]) if m else 0
        best_perm, best_val = None, 1e18
        for perm in itertools.permutations(range(n), m):
            s = sum(cost[i][perm[i]] for i in range(m))
            if s < best_val:
                best_val, best_perm = s, perm
        return list(range(m)), list(best_perm) if best_perm is not None else ([], [])

def get_bonds(atoms):
    """Get bonds using neighbor list and covalent radii."""
    cutoffs = [covalent_radii[n] * 1.2 for n in atoms.numbers]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    bonds = set()
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        for j in indices:
            if i < j:
                bonds.add((i, j))
    return list(bonds)

def _add_cap_at(point, direction, P, cap_mat, cap_template=None):
    dirn = direction.normalized()
    sR = P.cap_radius * P.cap_scale
    sL = P.cap_length * P.cap_scale
    q = cap_quaternion(dirn, P.cap_forward_axis, P.cap_roll_deg)
    if cap_template is None:
        bpy.ops.mesh.primitive_cone_add(
            vertices=24, radius1=sR, radius2=0.0, depth=sL,
            location=point, rotation=q.to_euler()
        )
        cap = bpy.context.object
        if len(cap.data.materials) == 0:
            cap.data.materials.append(cap_mat)
        else:
            cap.data.materials[0] = cap_mat
        if P.debug_mode:
            print(f"[CAP] Built-in cone at {tuple(point)}")
        return cap
    kind, ref = cap_template
    if kind == "OBJECT":
        cap = ref.copy()
        cap.data = ref.data.copy()
        cap.name = "bond_cap"
        bpy.context.scene.collection.objects.link(cap)
        cap.matrix_world = Matrix.Identity(4)
        cap.location = point
        cap.rotation_euler = q.to_euler()
        cap.scale = (sR, sR, sL)
        if len(cap.data.materials) == 0:
            cap.data.materials.append(cap_mat)
        else:
            cap.data.materials[0] = cap_mat
        if P.debug_mode:
            print(f"[CAP] Duplicated OBJECT '{ref.name}' at {tuple(point)}")
        return cap
    if kind == "COLLECTION":
        bpy.ops.object.collection_instance_add(collection=ref.name, location=(0, 0, 0))
        inst = bpy.context.object
        inst.name = "bond_cap"
        inst.matrix_world = Matrix.Identity(4)
        inst.location = point
        inst.rotation_euler = q.to_euler()
        inst.scale = (sR, sR, sL)
        if P.debug_mode:
            print(f"[CAP] Instanced COLLECTION '{ref.name}' at {tuple(point)}")
        return inst
    print("[CAP] Unknown template kind:", kind)
    return None
    



def parse_atoms_bonds(path, scale):
    """Parse coordinates and bonds from PDB or CIF file using ASE.
    Returns: atoms, bonds, coords, types."""
    atoms, bonds, coords, types = [], [], {}, {}

    ext = os.path.splitext(path)[1].lower()
    try:
        if ext == ".cif":
            molecule = ase_read(path, format="cif")
        else:
            molecule = ase_read(path)
    except StopIteration:
        raise ValueError(f"Il file {path} non contiene strutture leggibili o è vuoto")

    for i, atom in enumerate(molecule):
        idx = i+1  # 1-based indexing
        sym = atom.symbol
        pos = Vector(atom.position) * scale
        atoms.append((idx, sym, pos))
        coords[idx] = pos
        types[idx] = sym

    for (i1, i2) in get_bonds(molecule):
        if i1 < i2:
            bonds.append((i1+1, i2+1))  # Convert to 1-based
        else:
            bonds.append((i2+1, i1+1))

    return atoms, bonds, coords, types

def axis_vec(label: str) -> Vector:
    return {
        'X+': Vector((1,0,0)), 'X-': Vector((-1,0,0)),
        'Y+': Vector((0,1,0)), 'Y-': Vector((0,-1,0)),
        'Z+': Vector((0,0,1)), 'Z-': Vector((0,0,-1)),
    }[label]

def _load_cap_template(P):
    name = (P.cap_template_name or "").strip()
    if not name:
        if P.debug_mode:
            print("[CAP] No template name provided -> using built-in cone")
        return None
    try:
        obj = append_object(P.lib_path, name)
        if P.debug_mode:
            print(f"[CAP] Loaded OBJECT '{name}'")
        return ("OBJECT", obj)
    except Exception as e_obj:
        if P.debug_mode:
            print(f"[CAP] '{name}' not an OBJECT: {e_obj}")
    try:
        coll = append_collection(P.lib_path, name)
        if P.debug_mode:
            print(f"[CAP] Loaded COLLECTION '{name}'")
        return ("COLLECTION", coll)
    except Exception as e_col:
        print(f"[CAP] Template '{name}' not found as Object or Collection: {e_col}")
        return None

def _get_or_make_material(name, rgba):
    mat = bpy.data.materials.get(name)
    if not mat:
        mat = bpy.data.materials.new(name)
        mat.use_nodes = True
    nt = mat.node_tree
    bsdf = next((n for n in nt.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if bsdf:
        bsdf.inputs["Base Color"].default_value = (rgba[0], rgba[1], rgba[2], 1)
        bsdf.inputs["Roughness"].default_value = 0.45
    return mat

def _get_or_make_cap_material(P):
    # Se il nome è uno degli elementi comuni, usa la palette, altrimenti un grigio neutro
    if P.cap_mat_name in ["H","C","N","O","S","P","F","Cl","Br","I"]:
        col = getattr(P, f"col_{P.cap_mat_name}", (0.85, 0.85, 0.85, 1.0))
        return _get_or_make_material(f"Mol_{P.cap_mat_name}", col)
    else:
        return _get_or_make_material(P.cap_mat_name, (0.85, 0.85, 0.85, 1.0))


def cap_quaternion(dirn: Vector, forward_axis: str, roll_deg: float) -> Quaternion:
    forward_local = axis_vec(forward_axis)
    q_pre   = forward_local.rotation_difference(Vector((0,0,1)))
    q_align = Vector((0,0,1)).rotation_difference(dirn.normalized())
    q_roll  = Quaternion(dirn.normalized(), radians(roll_deg))
    return q_roll @ (q_align @ q_pre)

def choose_geometry_key(element: str, nn: int):
    e = element
    if nn <= 1 and (e == "H" or e in HALOGENS):
        return "Atom_sp"
    if nn == 2 and e in {"O", "S", "Se", "Te"}:
        return "Atom_bent"
    if e in {"N", "P", "As", "Sb"} and nn == 3:
        return "Atom_sp3"
    if e == "S" and nn >= 6:
        return "Atom_sp3d2"
    if e == "C":
        if nn >= 4: return "Atom_sp3"
        if nn == 3: return "Atom_sp2"
        if nn <= 2: return "Atom_sp"
    if nn >= 6: return "Atom_sp3d2"
    if nn == 5:  return "Atom_sp3d2"
    if nn == 4:  return "Atom_sp3"
    if nn == 3:  return "Atom_sp2"
    return "Atom_sp"

# Alias per compatibilità col vecchio naming
_axis_vec = axis_vec
