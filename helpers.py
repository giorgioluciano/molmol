import bpy, os, itertools
from mathutils import Vector, Matrix, Quaternion
from math import radians
from collections import defaultdict

# ============ ASE (fallback minimale) ============
try:
    from ase.io import read as ase_read
    ASE_AVAILABLE = True
except ImportError:
    ase_read = None
    ASE_AVAILABLE = False
    print("[Molymod] ASE not found.")

# ============ RDKIT (fonte primaria) ============
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolops
    RDKIT_AVAILABLE = True
    print("[Molymod] RDKit available!")
except ImportError:
    Chem = None
    AllChem = None
    RDKIT_AVAILABLE = False
    print("[Molymod] RDKit not found! Only SMILES/SDF/MOL2 supported.")

# ============ COSTANTI ============
GEN_COLLECTIONS = ["Atom_sp3", "Atom_sp2", "Atom_sp", "Atom_bent", "Atom_sp3d2"]
HALOGENS = {"F", "Cl", "Br", "I"}

MIN_VALENCE = {
    'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 3,
    'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
}

VALID_HYBRIDIZATIONS = {
    'C':  ['sp', 'sp2', 'sp3'],
    'N':  ['sp', 'sp2', 'sp3'],
    'O':  ['sp2', 'sp3', 'bent'],
    'S':  ['sp2', 'sp3', 'sp3d2'],
    'P':  ['sp2', 'sp3', 'sp3d2'],
    'H':  ['sp'],
    'F':  ['sp3'], 'Cl': ['sp3'], 'Br': ['sp3'], 'I': ['sp3'],
}

ATOM_VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'S': 1.80, 'Cl': 1.75, 'P': 1.80,
    'Br': 1.85, 'I': 1.98,
}

# ============ FORMATI SUPPORTATI ============
SUPPORTED_FORMATS = [".sdf", ".mol", ".mol2", ".smi", ".smiles"]

OPENBABEL_MESSAGE = """
[ERROR] Format '{ext}' is not supported by Molymod.

Molymod requires 3D molecular files with explicit bond information.
Supported formats: SDF (3D), MOL2 (3D), SMILES

Convert your file with OpenBabel:
  obabel yourfile{ext}  -O yourfile.sdf --gen3d
  obabel yourfile.pdb   -O yourfile.sdf --gen3d
  obabel yourfile.xyz   -O yourfile.sdf --gen3d
  obabel yourfile.cif   -O yourfile.sdf --gen3d
  obabel yourfile.mol   -O yourfile.sdf --gen3d  (2D to 3D)

Download OpenBabel : https://openbabel.org
3D SDF from PubChem: https://pubchem.ncbi.nlm.nih.gov
3D SDF from ChemSpider: https://www.chemspider.com
"""

MOL_2D_MESSAGE = """
[ERROR] File '{name}' has 2D coordinates (Z=0).

Molymod needs 3D coordinates to build physical models.

Convert to 3D with OpenBabel:
  obabel {name} -O out_3d.sdf --gen3d

Or download 3D SDF directly from:
  PubChem  : https://pubchem.ncbi.nlm.nih.gov
  ChemSpider: https://www.chemspider.com
"""

# ============ BLENDER UTILITIES ============
def abspath(path):
    return bpy.path.abspath(path)

def ensure_hidden_bucket():
    name = "_MolymodHiddenObjs"
    coll = bpy.data.collections.get(name)
    if not coll:
        coll = bpy.data.collections.new(name)
    return coll

def unlink_collection_everywhere(coll):
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
    holes = []
    with bpy.data.libraries.load(abspath(lib_path), link=False) as (src, dst):
        cand = [n for n in src.objects if n.startswith(f"{coll_key}_hole")]
    for h in cand:
        ob = append_object(lib_path, h)
        if ob.location.length > 1e-9:
            holes.append(ob.location.normalized())
    return holes

# ============ MATH UTILITIES ============
def kabsch_rotation(from_vecs, to_vecs):
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
    s = src.normalized(); d = dst.normalized()
    q = s.rotation_difference(d)
    return q.to_matrix().to_4x4()

def hungarian_assign(cost):
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

def axis_vec(label: str) -> Vector:
    return {
        'X+': Vector((1,0,0)), 'X-': Vector((-1,0,0)),
        'Y+': Vector((0,1,0)), 'Y-': Vector((0,-1,0)),
        'Z+': Vector((0,0,1)), 'Z-': Vector((0,0,-1)),
    }[label]

# ============ GEOMETRY ============
def choose_geometry_key(element: str, nn_holes: int):
    if nn_holes <= 1:
        candidate = "sp"
    elif nn_holes == 2:
        candidate = "bent" if element in {"O", "S", "Se", "Te"} else "sp"
    elif nn_holes == 3:
        candidate = "sp2"
    elif nn_holes == 4:
        candidate = "sp3"
    else:
        candidate = "sp3d2"

    valid = VALID_HYBRIDIZATIONS.get(element, ['sp3'])
    if candidate not in valid:
        print(f"[WARNING] {element} cannot be {candidate}! Falling back to sp3.")
        # Trova il più grande valido
        for fallback in ['sp3', 'sp2', 'sp']:
            if fallback in valid:
                return fallback
    return candidate

def assign_double_bond_holes(hole_vecs_world, neighbors_dirs):
    slots = []
    for n_idx, dirn, order in neighbors_dirs:
        n_slots = round(order)
        for k in range(n_slots):
            slots.append((n_idx, dirn))

    if not slots or not hole_vecs_world:
        return {}

    n_holes = len(hole_vecs_world)
    use = min(n_holes, len(slots))

    cost = [
        [1.0 - max(-1.0, min(1.0, hole_vecs_world[h].dot(slots[s][1])))
         for s in range(use)]
        for h in range(n_holes)
    ]

    try:
        import numpy as np
        import scipy.optimize as spopt
        row_ind, col_ind = spopt.linear_sum_assignment(np.array(cost, dtype=float)[:use, :use])
    except Exception:
        row_ind = list(range(use))
        col_ind = list(range(use))

    result = {}
    for r, c in zip(row_ind, col_ind):
        if c >= len(slots): continue
        n_idx = slots[c][0]
        if n_idx not in result:
            result[n_idx] = []
        result[n_idx].append(hole_vecs_world[r])

    return result

# ============ VALIDATION ============
def detect_missing_hydrogens(atoms_list, bonds, bond_orders, types):
    missing_H = {}
    for idx, sym, pos in atoms_list:
        if sym == 'H':
            continue
        neighs = [t for s, t in bonds if s == idx] + \
                 [s for s, t in bonds if t == idx]
        bonds_used = sum(
            round(bond_orders.get((min(idx,n), max(idx,n)), 1))
            for n in neighs
        )
        valence = MIN_VALENCE.get(sym, 4)
        n_H_needed = valence - bonds_used
        if n_H_needed > 0:
            missing_H[idx] = n_H_needed
            print(f"[INFO] Atom {idx} ({sym}): needs {n_H_needed} H")
    return missing_H

# ============ PARSE - UNICA FUNZIONE ============
def parse_atoms_bonds(path, scale):
    """
    Unico entry point per il parsing.
    Accetta SOLO: SDF/MOL 3D, MOL2, SMILES
    Per tutto il resto: errore con istruzioni OpenBabel
    """
    ext = os.path.splitext(path)[1].lower()
    name = os.path.basename(path)

    print(f"\n{'='*50}")
    print(f"PARSING: {name}")
    print(f"{'='*50}")

    # ============ CHECK FORMATO ============
    if ext not in SUPPORTED_FORMATS:
        print(OPENBABEL_MESSAGE.format(ext=ext))
        return None

    # ============ CHECK RDKIT ============
    if not RDKIT_AVAILABLE:
        print("[ERROR] RDKit is required but not installed.")
        print("Install with: pip install rdkit")
        print("(Use Blender's Python pip)")
        return None

    # ============ LEGGI CON RDKIT ============
    mol = None
    try:
        if ext in (".mol", ".sdf"):
            suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
            mol = next((m for m in suppl if m is not None), None)
        elif ext == ".mol2":
            mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)
        elif ext in (".smi", ".smiles"):
            with open(path) as f:
                smi = f.readline().strip().split()[0]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                AllChem.UFFOptimizeMolecule(mol)
                print("[OK] SMILES: 3D coordinates generated by RDKit")
    except Exception as e:
        print(f"[ERROR] RDKit failed to read file: {e}")
        return None

    if mol is None:
        print(f"[ERROR] Could not parse {name}")
        return None

    # ============ CHECK 2D ============
    if ext in (".mol", ".sdf"):
        if mol.GetNumConformers() == 0:
            print(MOL_2D_MESSAGE.format(name=name))
            return None
        
