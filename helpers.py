import bpy, os, itertools
from mathutils import Vector, Matrix, Quaternion
from math import radians

try:
    from ase.io import read as ase_read
    from ase.neighborlist import NeighborList
    from ase.data import covalent_radii
    ASE_AVAILABLE = True
except ImportError:
    ase_read = None
    NeighborList = None
    covalent_radii = None
    ASE_AVAILABLE = False
    print("[Molymod] ASE not found: PDB-only mode active (CONECT records).")

# RDKit for accurate bond orders from SDF/MOL/MOL2/SMILES
try:
    from rdkit import Chem
    from rdkit.Chem import rdmolops
    RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    RDKIT_AVAILABLE = False
    print("[Molymod] RDKit not found: bond orders will use heuristic estimation.")

GEN_COLLECTIONS = ["Atom_sp3", "Atom_sp2", "Atom_sp", "Atom_bent", "Atom_sp3d2"]
HALOGENS = {"F", "Cl", "Br", "I"}

# Thresholds for bond order estimation based on inter-atomic distances (Angstrom)
BOND_ORDER_THRESHOLDS = {
    ('C', 'C'): [(1.60, 1), (1.42, 1.5), (1.34, 2), (1.20, 3)],
    ('C', 'N'): [(1.50, 1), (1.35, 1.5), (1.28, 2), (1.16, 3)],
    ('C', 'O'): [(1.45, 1), (1.35, 1.5), (1.22, 2)],
    ('C', 'S'): [(1.85, 1), (1.65, 2)],
    ('N', 'N'): [(1.45, 1), (1.25, 2), (1.10, 3)],
    ('N', 'O'): [(1.45, 1), (1.22, 2)],
    ('O', 'O'): [(1.48, 1), (1.21, 2)],
}

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

def _guess_bond_order(sym1, sym2, distance):
    """Estimate bond order from distance using empirical thresholds."""
    key = tuple(sorted([sym1, sym2]))
    thresholds = BOND_ORDER_THRESHOLDS.get(key, [])
    for (cutoff, order) in sorted(thresholds, key=lambda x: x[0]):
        if distance <= cutoff:
            return order
    return 1

def _parse_bond_orders_rdkit(path, atom_index_offset=1):
    """
    Read bond orders from SDF/MOL/MOL2/SMILES using RDKit.
    Returns dict {(i1, i2): order} with 1-based indices.
    order values: 1=single, 2=double, 3=triple, 1.5=aromatic
    """
    ext = os.path.splitext(path)[1].lower()
    mol = None
    try:
        if ext in (".sdf", ".mol"):
            suppl = Chem.SDMolSupplier(path, removeHs=False)
            mol = next((m for m in suppl if m is not None), None)
        elif ext == ".mol2":
            mol = Chem.MolFromMol2File(path, removeHs=False)
        elif ext in (".smi", ".smiles"):
            with open(path) as f:
                smi = f.readline().strip().split()[0]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                from rdkit.Chem import AllChem
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    except Exception as e:
        print(f"[Molymod] RDKit failed to read {path}: {e}")
        return {}

    if mol is None:
        print(f"[Molymod] RDKit: could not parse {path}")
        return {}

    bond_orders = {}
    RDKIT_ORDER = {
        Chem.rdchem.BondType.SINGLE:    1,
        Chem.rdchem.BondType.DOUBLE:    2,
        Chem.rdchem.BondType.TRIPLE:    3,
        Chem.rdchem.BondType.AROMATIC:  1.5,
    }
    for bond in mol.GetBonds():
        i1 = bond.GetBeginAtomIdx() + atom_index_offset
        i2 = bond.GetEndAtomIdx()   + atom_index_offset
        order = RDKIT_ORDER.get(bond.GetBondType(), 1)
        bond_orders[(min(i1,i2), max(i1,i2))] = order

    print(f"[Molymod] RDKit: read {len(bond_orders)} bond orders from {os.path.basename(path)}")
    return bond_orders

def assign_double_bond_holes(hole_vecs_world, neighbors_dirs):
    """
    Assegna i fori dell'atomo ai neighbor tenendo conto degli ordini di legame.
    
    Per un legame doppio verso neighbor N, riserva 2 fori (quelli più allineati
    con la direzione verso N). Per legame singolo, 1 foro.
    
    Parametri:
      hole_vecs_world: lista di Vector (fori dell'atomo in world-space, già ruotati)
      neighbors_dirs:  lista di (neighbor_idx, direction_Vector, bond_order)
      
    Ritorna:
      dict {neighbor_idx: [hole_vec, ...]}  (lista di fori assegnati a ciascun neighbor)
    """
    # Espandi i neighbor in "slot" in base all'ordine di legame
    slots = []  # (neighbor_idx, dir)
    for n_idx, dirn, order in neighbors_dirs:
        n_slots = max(1, round(order)) if order != 1.5 else 2
        for k in range(n_slots):
            slots.append((n_idx, dirn))

    if not slots or not hole_vecs_world:
        return {}

    n_holes = len(hole_vecs_world)
    n_slots  = len(slots)
    use = min(n_holes, n_slots)

    # Matrice costo: 1 - dot(hole, dir) per ogni (hole, slot)
    cost = [
        [1.0 - max(-1.0, min(1.0, hole_vecs_world[h].dot(slots[s][1])))
         for s in range(use)]
        for h in range(n_holes)
    ]
    
    # Hungarian
    try:
        import numpy as np
        import scipy.optimize as spopt
        row_ind, col_ind = spopt.linear_sum_assignment(np.array(cost, dtype=float)[:use, :use])
    except Exception:
        row_ind = list(range(use))
        col_ind = list(range(use))

    # Raggruppa i fori assegnati per neighbor
    result = {}
    for r, c in zip(row_ind, col_ind):
        if c >= len(slots): continue
        n_idx = slots[c][0]
        if n_idx not in result:
            result[n_idx] = []
        result[n_idx].append(hole_vecs_world[r])

    return result

def parse_atoms_bonds(path, scale):
    """Parse atoms, bonds, and bond orders from molecular file using ASE + RDKit.
    
    Pipeline:
    1. Open file (PDB/XYZ/CIF/SDF/MOL)
    2. Extract atoms + coordinates
    3. Extract bonds (CONECT if available, else NeighborList)
    4. Extract bond orders (RDKit if available, CONECT multiplicity, else heuristic)
    5. Return: atoms_list, bonds, coords, types, bond_orders
    """
    ext = os.path.splitext(path)[1].lower()

    # STEP 1: Open file
    try:
        if ext == ".cif":
            molecule = ase_read(path, format="cif")
        else:
            molecule = ase_read(path)
    except StopIteration:
        raise ValueError(f"File {path} contains no readable structures or is empty")
    except Exception as e:
        raise ValueError(f"Failed to read {path}: {e}")

    if len(molecule) == 0:
        raise ValueError(f"File {path} contains no atoms")

    # STEP 2: Extract atoms
    atoms_list = []
    coords = {}
    types = {}
    for i, atom in enumerate(molecule):
        idx = i + 1
        sym = atom.symbol
        pos = Vector(atom.position) * scale
        atoms_list.append((idx, sym, pos))
        coords[idx] = pos
        types[idx] = sym
    print(f"[Parse] Found {len(atoms_list)} atoms in {os.path.basename(path)}")

    # STEP 3: Extract bonds
    bonds = []
    if ext in (".pdb", ".ent"):
        # PDB: use CONECT records
        bond_set = set()
        try:
            with open(path, "r") as f:
                for line in f:
                    if line.startswith("CONECT"):
                        fields = line.split()
                        if len(fields) < 3:
                            continue
                        origin = int(fields[1])
                        for target in fields[2:]:
                            t = int(target)
                            pair = (min(origin, t), max(origin, t))
                            bond_set.add(pair)
            bonds = list(bond_set)
            print(f"[Parse] Found {len(bonds)} bonds from CONECT records")
        except Exception as e:
            print(f"[Parse] WARNING: Could not read CONECT records: {e}")

    if not bonds and ASE_AVAILABLE:
        # Fallback: use ASE NeighborList
        for (i1, i2) in get_bonds(molecule):
            if i1 < i2:
                bonds.append((i1 + 1, i2 + 1))
            else:
                bonds.append((i2 + 1, i1 + 1))
        print(f"[Parse] Computed {len(bonds)} bonds from NeighborList")

    # STEP 4: Extract/estimate bond orders
    bond_orders = {}

    # Tentativo 1: RDKit (preciso)
    if RDKIT_AVAILABLE and ext in (".sdf", ".mol", ".mol2", ".smi", ".smiles"):
        bond_orders = _parse_bond_orders_rdkit(path, atom_index_offset=1)

    # Tentativo 2: PDB CONECT multiplicity
    if not bond_orders and ext in (".pdb", ".ent"):
        raw_conect = {}
        try:
            with open(path, "r") as f:
                for line in f:
                    if line.startswith("CONECT"):
                        fields = line.split()
                        if len(fields) < 3:
                            continue
                        origin = int(fields[1])
                        for target in fields[2:]:
                            t = int(target)
                            pair = (min(origin, t), max(origin, t))
                            raw_conect[pair] = raw_conect.get(pair, 0) + 1
            bond_orders = {pair: min(count, 3) for pair, count in raw_conect.items()}
            print(f"[Parse] Read {len(bond_orders)} bond orders from CONECT multiplicity")
        except Exception as e:
            print(f"[Parse] WARNING: Could not read CONECT bond orders: {e}")

    # Tentativo 3: Euristica dalla distanza
    if not bond_orders:
        for (i1, i2) in bonds:
            sym1 = types[i1]
            sym2 = types[i2]
            dist = (coords[i2] - coords[i1]).length / scale
            order = _guess_bond_order(sym1, sym2, dist)
            bond_orders[(min(i1,i2), max(i1,i2))] = order
        print(f"[Parse] Estimated bond orders heuristically for {len(bond_orders)} bonds")
    else:
        # Completa i legami mancanti con ordine 1
        for (i1, i2) in bonds:
            k = (min(i1,i2), max(i1,i2))
            if k not in bond_orders:
                bond_orders[k] = 1

    print(f"[Parse] Final: {len(bond_orders)} bonds with orders")
    return atoms_list, bonds, coords, types, bond_orders

def axis_vec(label: str) -> Vector:
    return {
        'X+': Vector((1,0,0)), 'X-': Vector((-1,0,0)),
        'Y+': Vector((0,1,0)), 'Y-': Vector((0,-1,0)),
        'Z+': Vector((0,0,1)), 'Z-': Vector((0,0,-1)),
    }[label]

def choose_geometry_key(element: str, nn: int):
    """Choose geometry key based on element and neighbor count.
    Returns a key like 'sp2', 'sp3', etc. (without 'Atom_' prefix).
    """
    e = element
    if nn <= 1 and (e == "H" or e in HALOGENS):
        return "sp"
    if nn == 2 and e in {"O", "S", "Se", "Te"}:
        return "bent"
    if e in {"N", "P", "As", "Sb"} and nn == 3:
        return "sp2"
    if e in {"N", "P", "As", "Sb"} and nn == 4:
        return "sp3"
    if e == "S" and nn >= 6:
        return "sp3d2"
    if e == "C":
        if nn >= 4: return "sp3"
        if nn == 3: return "sp2"
        if nn <= 2: return "sp"
    if nn >= 6: return "sp3d2"
    if nn == 5: return "sp3d2"
    if nn == 4: return "sp3"
    if nn == 3: return "sp2"
    return "sp"

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
    if P.cap_mat_name in ["H","C","N","O","S","P","F","Cl","Br","I"]:
        col = getattr(P, f"col_{P.cap_mat_name}", (0.85, 0.85, 0.85, 1.0))
        return _get_or_make_material(f"Mol_{P.cap_mat_name}", col)
    else:
        return _get_or_make_material(P.cap_mat_name, (0.85, 0.85, 0.85, 1.0))

def cap_quaternion(dirn: Vector, forward_axis: str, roll_deg: float) -> Quaternion:
    forward_local = axis_vec(forward_axis)
    q_pre = forward_local.rotation_difference(Vector((0,0,1)))
    q_align = Vector((0,0,1)).rotation_difference(dirn.normalized())
    q_roll = Quaternion(dirn.normalized(), radians(roll_deg))
    return q_roll @ (q_align @ q_pre)

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

# Alias for compatibility
_axis_vec = axis_vec
