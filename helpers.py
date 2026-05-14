import bpy, os, itertools
from mathutils import Vector, Matrix, Quaternion
from math import radians
from collections import defaultdict

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

# Valenze standard per elementi
MIN_VALENCE = {
    'C': 4,  'N': 3,  'O': 2,  'S': 2,  'P': 3,
    'H': 1,  'F': 1,  'Cl': 1, 'Br': 1, 'I': 1,
}

# Ibridazioni VALIDE per elemento (C NON può fare sp3d2!)
VALID_HYBRIDIZATIONS = {
    'C': ['sp', 'sp2', 'sp3'],
    'N': ['sp', 'sp2', 'sp3'],
    'O': ['sp2', 'sp3', 'bent'],
    'S': ['sp2', 'sp3', 'sp3d2'],
    'P': ['sp2', 'sp3', 'sp3d2'],
    'H': ['sp'],
    'F': ['sp3'], 'Cl': ['sp3'], 'Br': ['sp3'], 'I': ['sp3'],
}

# Raggi atomici van der Waals (Å)
ATOM_VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'F': 1.47, 'S': 1.80, 'Cl': 1.75, 'P': 1.80,
    'Br': 1.85, 'I': 1.98,
}

# Thresholds for bond order estimation
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
    return holes

def kabsch_rotation(from_vecs, to_vecs):
    """Calculate optimal rotation matrix (Kabsch algorithm)."""
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
    """Optimal assignment (Hungarian method)."""
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
    """Read bond orders from SDF/MOL/MOL2/SMILES using RDKit."""
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
        print(f"[Molymod] RDKit failed: {e}")
        return {}

    if mol is None:
        return {}

    bond_orders = {}
    RDKIT_ORDER = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.AROMATIC: 1.5,
    }
    for bond in mol.GetBonds():
        i1 = bond.GetBeginAtomIdx() + atom_index_offset
        i2 = bond.GetEndAtomIdx() + atom_index_offset
        order = RDKIT_ORDER.get(bond.GetBondType(), 1)
        bond_orders[(min(i1,i2), max(i1,i2))] = order

    print(f"[OK] RDKit: {len(bond_orders)} bond orders")
    return bond_orders

def find_rings(bonds, max_size=20):
    """Find all simple rings in the molecule using DFS."""
    graph = defaultdict(list)
    for s, t in bonds:
        graph[s].append(t)
        graph[t].append(s)
    
    rings = []
    visited_edges = set()
    
    def dfs_ring(start, current, parent, path):
        if len(path) > max_size:
            return
        for neighbor in graph[current]:
            if neighbor == parent:
                continue
            edge = (min(current, neighbor), max(current, neighbor))
            if edge in visited_edges:
                continue
            if neighbor == start and len(path) >= 3:
                rings.append(list(path))
                return
            if neighbor not in path:
                path.append(neighbor)
                dfs_ring(start, neighbor, current, path)
                path.pop()
    
    all_nodes = set()
    for s, t in bonds:
        all_nodes.add(s)
        all_nodes.add(t)
    
    for node in sorted(all_nodes):
        dfs_ring(node, node, None, [node])
    
    # Remove duplicates
    unique_rings = []
    seen = set()
    for ring in rings:
        min_idx = ring.index(min(ring))
        normalized = tuple(ring[min_idx:] + ring[:min_idx])
        normalized_rev = tuple(reversed(normalized))
        canonical = min(normalized, normalized_rev)
        if canonical not in seen:
            seen.add(canonical)
            unique_rings.append(list(canonical))
    
    return unique_rings

def is_aromatic_ring(ring, bond_orders):
    """Check if all bonds in ring are aromatic (order ~1.5)"""
    for i in range(len(ring)):
        s = ring[i]
        t = ring[(i+1) % len(ring)]
        order = bond_orders.get((min(s,t), max(s,t)), 1)
        if abs(order - 1.5) > 0.2:
            return False
    return True

def find_fusion_bonds(rings):
    """Find bonds shared between rings (fusion bonds)"""
    bond_count = defaultdict(int)
    for ring in rings:
        for i in range(len(ring)):
            s = ring[i]
            t = ring[(i+1) % len(ring)]
            bond_count[(min(s,t), max(s,t))] += 1
    
    fusion = {bond for bond, count in bond_count.items() if count > 1}
    return fusion

def kekulize_rings(rings, bond_orders):
    """Convert aromatic rings (order=1.5) to alternating Kekulé pattern (1,2,1,2...)"""
    aromatic_rings = [r for r in rings if is_aromatic_ring(r, bond_orders)]
    if not aromatic_rings:
        return
    
    print(f"[INFO] Kekulizing {len(aromatic_rings)} aromatic rings")
    fusion_bonds = find_fusion_bonds(aromatic_rings)
    aromatic_rings.sort(key=len, reverse=True)
    processed_bonds = set()
    
    for ring in aromatic_rings:
        start_idx = 0
        fusion_found = False
        
        for i in range(len(ring)):
            s = ring[i]
            t = ring[(i+1) % len(ring)]
            key = (min(s,t), max(s,t))
            
            if key in fusion_bonds and key in processed_bonds:
                start_idx = (i + 1) % len(ring)
                fusion_found = True
                break
        
        for j in range(len(ring)):
            i = (start_idx + j) % len(ring)
            s = ring[i]
            t = ring[(i+1) % len(ring)]
            key = (min(s,t), max(s,t))
            
            if key in processed_bonds:
                continue
            
            if key in fusion_bonds and not fusion_found:
                new_order = 1
            else:
                new_order = 1 if j % 2 == 0 else 2
            
            bond_orders[key] = new_order
            processed_bonds.add(key)
    
    print(f"[OK] Kekulization: {len(processed_bonds)} bonds alternated")

def assign_double_bond_holes(hole_vecs_world, neighbors_dirs):
    """Assign holes to neighbors considering bond orders."""
    slots = []
    for n_idx, dirn, order in neighbors_dirs:
        n_slots = round(order)
        for k in range(n_slots):
            slots.append((n_idx, dirn))

    if not slots or not hole_vecs_world:
        return {}

    n_holes = len(hole_vecs_world)
    n_slots = len(slots)
    use = min(n_holes, n_slots)

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

def detect_2d_coords(coords):
    """Detect if all atoms have Z=0 (2D file)"""
    all_z_zero = all(abs(pos.z) < 0.001 for pos in coords.values())
    if all_z_zero:
        print("[WARNING] 2D coordinates detected (Z=0). Molecule will be flat.")
        print("[SUGGESTION] Use OpenBabel: obabel file.mol -O out.sdf --gen3d")
        return True
    return False

def detect_missing_hydrogens(atoms_list, bonds, bond_orders, types):
    """Calculate missing hydrogens for each atom."""
    missing_H = {}
    for idx, sym, pos in atoms_list:
        if sym == 'H':
            continue
        
        neighs = [t for s, t in bonds if s == idx] + [s for s, t in bonds if t == idx]
        bonds_used = sum(round(bond_orders.get((min(idx, n), max(idx, n)), 1)) for n in neighs)
        
        valence = MIN_VALENCE.get(sym, 4)
        n_H_needed = valence - bonds_used
        
        if n_H_needed > 0:
            missing_H[idx] = n_H_needed
            print(f"[WARNING] Atom {idx} ({sym}) missing {n_H_needed} hydrogens")
    
    return missing_H

def parse_atoms_bonds(path, scale):
    """
    Parse atoms, bonds, and bond orders with full validation.
    Returns: atoms_list, bonds, coords, types, bond_orders
    """
    ext = os.path.splitext(path)[1].lower()
    print(f"\n{'='*50}")
    print(f"PARSING: {os.path.basename(path)}")
    print(f"{'='*50}")

    # STEP 1: Load file
    try:
        if ext == ".cif":
            molecule = ase_read(path, format="cif")
        else:
            molecule = ase_read(path)
    except Exception as e:
        print(f"[ERROR] Cannot read file: {e}")
        return None

    if len(molecule) == 0:
        print(f"[ERROR] File contains no atoms")
        return None

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
    print(f"[OK] Found {len(atoms_list)} atoms")

    # STEP 3: Check 2D vs 3D
    is_2d = detect_2d_coords(coords)

    # STEP 4: Extract bonds
    bonds = []
    if ext in (".pdb", ".ent"):
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
            print(f"[OK] Found {len(bonds)} bonds from CONECT")
        except Exception as e:
            print(f"[WARNING] Could not read CONECT: {e}")

    if not bonds and ASE_AVAILABLE:
        print("[INFO] Computing bonds from distances...")
        for (i1, i2) in get_bonds(molecule):
            bonds.append((i1 + 1, i2 + 1) if i1 < i2 else (i2 + 1, i1 + 1))
        print(f"[OK] Computed {len(bonds)} bonds")
    elif not bonds:
        print("[ERROR] No bonds and ASE unavailable!")
        return None

    # STEP 5: Extract/estimate bond orders
    bond_orders = {}
    
    if RDKIT_AVAILABLE and ext in (".sdf", ".mol", ".mol2", ".smi"):
        bond_orders = _parse_bond_orders_rdkit(path)
    
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
            print(f"[OK] Read {len(bond_orders)} bond orders from CONECT")
        except: pass
    
    if not bond_orders:
        print("[INFO] Using distance heuristic for bond orders")
        for (i1, i2) in bonds:
            dist = (coords[i2] - coords[i1]).length / scale
            order = _guess_bond_order(types[i1], types[i2], dist)
            bond_orders[(min(i1, i2), max(i1, i2))] = order
    
    for (i1, i2) in bonds:
        k = (min(i1, i2), max(i1, i2))
        if k not in bond_orders:
            bond_orders[k] = 1
    
    print(f"[OK] Bond orders: {len(bond_orders)} bonds")

    # STEP 6: Kekulize aromatics
    rings = find_rings(bonds)
    if rings:
        print(f"[INFO] Found {len(rings)} rings")
        kekulize_rings(rings, bond_orders)

    # STEP 7: Summary
    print(f"\n{'='*50}")
    print(f"VALIDATION SUMMARY")
    print(f"{'='*50}")
    print(f"Atoms: {len(atoms_list)}")
    print(f"Bonds: {len(bonds)}")
    print(f"3D coords: {'YES' if not is_2d else 'NO (2D flat)'}")
    print(f"Rings: {len(rings) if rings else 0}")
    print(f"{'='*50}\n")

    return atoms_list, bonds, coords, types, bond_orders

def axis_vec(label: str) -> Vector:
    return {
        'X+': Vector((1,0,0)), 'X-': Vector((-1,0,0)),
        'Y+': Vector((0,1,0)), 'Y-': Vector((0,-1,0)),
        'Z+': Vector((0,0,1)), 'Z-': Vector((0,0,-1)),
    }[label]

def choose_geometry_key(element: str, nn_holes: int):
    """Choose geometry based on holes needed, respecting element constraints."""
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
    
    # Validate for element
    valid = VALID_HYBRIDIZATIONS.get(element, ['sp3'])
    if candidate not in valid:
        print(f"[ERROR] {element} cannot be {candidate}! Using sp3 fallback.")
        return "sp3"
    
    return candidate

def _load_cap_template(P):
    name = (P.cap_template_name or "").strip()
    if not name:
        return None
    try:
        obj = append_object(P.lib_path, name)
        return ("OBJECT", obj)
    except: pass
    try:
        coll = append_collection(P.lib_path, name)
        return ("COLLECTION", coll)
    except: pass
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
    return _get_or_make_material(P.cap_mat_name, (0.85, 0.85, 0.85, 1.0))

def cap_quaternion(dirn: Vector, forward_axis: str, roll_deg: float):
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
        return cap
    
    if kind == "COLLECTION":
        bpy.ops.object.collection_instance_add(collection=ref.name, location=(0, 0, 0))
        inst = bpy.context.object
        inst.name = "bond_cap"
        inst.matrix_world = Matrix.Identity(4)
        inst.location = point
        inst.rotation_euler = q.to_euler()
        inst.scale = (sR, sR, sL)
        return inst
    
    return None

_axis_vec = axis_vec
