import bpy, os, itertools
from mathutils import Vector, Matrix, Quaternion
from math import radians

# Import ASE for atomic file parsing and bond inference
from ase.io import read as ase_read

GEN_COLLECTIONS = ["Atom_sp3", "Atom_sp2", "Atom_sp", "Atom_bent", "Atom_sp3d2"]
HALOGENS = {"F", "Cl", "Br", "I"}


import bpy

class MOLYMOD_PG_Settings(bpy.types.PropertyGroup):
    # Libreria & file PDB
    lib_path: bpy.props.StringProperty(
        name="Library .blend", subtype='FILE_PATH', default="//molymod_library.blend"
    )
    molecule_path: bpy.props.StringProperty(
        name="Molecule file", subtype='FILE_PATH', default=""
    )
    scale: bpy.props.FloatProperty(
        name="Scale", default=3.0, min=0.001, soft_min=0.01, soft_max=50.0
    )
    compact_factor: bpy.props.FloatProperty(
        name="Compact Factor", default=1.0, min=0.1, max=2.0,
        description="Multiply all coordinates to shrink/expand molecule distances"
    )
    clear_previous: bpy.props.BoolProperty(name="Clear previous build", default=True)
    debug_mode: bpy.props.BoolProperty(name="Debug log & guides", default=False)

    # Bonds
    bond_radius: bpy.props.FloatProperty(name="Bond Radius", default=0.1, min=0.001, soft_max=1.0)
    bond_gap_each_side: bpy.props.FloatProperty(name="Gap Each Side", default=0.2, min=0.0, soft_max=3.0)
    bond_start_offset: bpy.props.FloatProperty(name="Start Offset", default=0.0, soft_min=-2.0, soft_max=2.0)
    bond_end_offset: bpy.props.FloatProperty(name="End Offset", default=0.0, soft_min=-2.0, soft_max=2.0)
    bond_length_factor: bpy.props.FloatProperty(name="Length Factor", default=1.0, min=0.1, max=1.5)
    bond_vertices: bpy.props.IntProperty(name="Vertices", default=30, min=3, soft_max=128)
    bond_mat_name: bpy.props.StringProperty(name="Bond Material", default="MolBond")

    # Caps
    use_caps: bpy.props.BoolProperty(name="Add Caps on Bonds", default=False)
    cap_template_name: bpy.props.StringProperty(name="Cap Template", default="")
    cap_scale: bpy.props.FloatProperty(name="Cap Scale", default=1.0, min=0.01, soft_max=10.0)
    cap_radius: bpy.props.FloatProperty(name="Base Radius", default=0.06, min=0.001, soft_max=1.0)
    cap_length: bpy.props.FloatProperty(name="Base Length", default=0.20, min=0.001, soft_max=2.0)
    cap_roll_deg: bpy.props.FloatProperty(name="Cap Roll (deg)", default=0.0, soft_min=-180.0, soft_max=180.0)
    cap_forward_axis: bpy.props.EnumProperty(
        name="Cap Forward Axis",
        items=[('Z+','Z+',''),('Z-','Z-',''),('X+','X+',''),('X-','X-',''),('Y+','Y+',''),('Y-','Y-','')],
        default='Z+'
    )
    cap_offset: bpy.props.FloatProperty(name="Cap Offset", default=0.10, min=0.0, soft_max=2.0)
    cap_start_offset: bpy.props.FloatProperty(name="Start Cap Offset", default=0.0, soft_min=-2.0, soft_max=2.0)
    cap_end_offset: bpy.props.FloatProperty(name="End Cap Offset", default=0.0, soft_min=-2.0, soft_max=2.0)
    cap_mat_name: bpy.props.StringProperty(name="Cap Material", default="MolCap")

    # Hydrogens / monovalent
    H_forward_axis: bpy.props.EnumProperty(
        name="H Forward Axis",
        items=[('Z+','Z+',''),('Z-','Z-',''),('X+','X+',''),('X-','X-',''),('Y+','Y+',''),('Y-','Y-','')],
        default='X+'
    )
    H_roll_deg: bpy.props.FloatProperty(name="H Roll (deg)", default=0.0, soft_min=-180.0, soft_max=180.0)

    # Palette dei colori
    col_H:  bpy.props.FloatVectorProperty(name="H",  subtype='COLOR', size=4, default=(1,1,1,1))
    col_C:  bpy.props.FloatVectorProperty(name="C",  subtype='COLOR', size=4, default=(0.2,0.2,0.2,1))
    col_N:  bpy.props.FloatVectorProperty(name="N",  subtype='COLOR', size=4, default=(0.1,0.3,0.9,1))
    col_O:  bpy.props.FloatVectorProperty(name="O",  subtype='COLOR', size=4, default=(0.9,0.1,0.1,1))
    col_S:  bpy.props.FloatVectorProperty(name="S",  subtype='COLOR', size=4, default=(1.0,0.85,0.1,1))
    col_P:  bpy.props.FloatVectorProperty(name="P",  subtype='COLOR', size=4, default=(1.0,0.5,0.0,1))
    col_F:  bpy.props.FloatVectorProperty(name="F",  subtype='COLOR', size=4, default=(0.1,0.8,0.1,1))
    col_Cl: bpy.props.FloatVectorProperty(name="Cl", subtype='COLOR', size=4, default=(0.0,0.6,0.0,1))
    col_Br: bpy.props.FloatVectorProperty(name="Br", subtype='COLOR', size=4, default=(0.6,0.2,0.0,1))
    col_I:  bpy.props.FloatVectorProperty(name="I",  subtype='COLOR', size=4, default=(0.5,0.0,0.5,1))



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

