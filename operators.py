import bpy, os
from math import radians
from mathutils import Vector, Matrix, Quaternion

from .helpers import (
    GEN_COLLECTIONS, HALOGENS,
    abspath, ensure_hidden_bucket, unlink_collection_everywhere,
    append_collection, append_object, load_hole_dirs,
    kabsch_rotation, align_one_vector, hungarian_assign,
    parse_atoms_bonds, choose_geometry_key, axis_vec,
    _load_cap_template, _get_or_make_cap_material, _add_cap_at
)


def _make_bezier_bond(p0, p1, p2, p3, radius, name, context):
    """Crea un legame come curva Bezier cubica con bevel (tubo 3D).
    p0, p3 = basi dei cap (punti di ancoraggio)
    p1, p2 = maniglie (tangenti ai fori)
    radius  = raggio del tubo
    """
    curve_data = bpy.data.curves.new(name, type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.resolution_u = 12
    curve_data.bevel_depth = radius
    curve_data.bevel_resolution = 4
    curve_data.use_fill_caps = True

    spline = curve_data.splines.new('BEZIER')
    spline.bezier_points.add(1)  # 2 punti totali

    bp0 = spline.bezier_points[0]
    bp0.co             = p0
    bp0.handle_left    = p0  # handle non usato (inizio)
    bp0.handle_right   = p1  # maniglia di uscita
    bp0.handle_left_type  = 'FREE'
    bp0.handle_right_type = 'FREE'

    bp1 = spline.bezier_points[1]
    bp1.co             = p3
    bp1.handle_left    = p2  # maniglia di entrata
    bp1.handle_right   = p3  # handle non usato (fine)
    bp1.handle_left_type  = 'FREE'
    bp1.handle_right_type = 'FREE'

    obj = bpy.data.objects.new(name, curve_data)
    context.scene.collection.objects.link(obj)
    return obj


class MOLYMOD_OT_Build(bpy.types.Operator):
    bl_idname = "molymod.build"
    bl_label = "Build Molecule from File"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        P = context.scene.molymod_settings
        lib     = abspath(P.lib_path)
        molfile = abspath(P.molecule_path)

        if not os.path.isfile(lib):
            self.report({'ERROR'}, f"Library .blend not found: {lib}")
            return {'CANCELLED'}
        if not os.path.isfile(molfile):
            self.report({'ERROR'}, f"Molecule file not found: {molfile}")
            return {'CANCELLED'}

        if P.clear_previous:
            for o in list(bpy.data.objects):
                if o.name.startswith(("mol_", "bond_", "bond_cap", "DEBUG_")):
                    bpy.data.objects.remove(o, do_unlink=True)

        # Assicura che tutte le collection generiche siano appendate
        for key in GEN_COLLECTIONS:
            try:
                append_collection(P.lib_path, key)
            except Exception as e:
                print(f"[Molymod] '{key}' not found (ok if unused). {e}")

        # --- Parsing molecola via ASE ---
        atoms, bonds, coords, types = parse_atoms_bonds(molfile, P.scale)
        # Placeholder: in futuro potrai far restituire anche bond_orders
        bond_orders = None

        # Compattazione opzionale
        if abs(P.compact_factor - 1.0) > 1e-9:
            for k in coords:
                coords[k] = coords[k] * P.compact_factor

        hole_cache    = {}
        placed        = {}
        # {(idx_atomo, idx_vicino): Vector normalizzato dal foro verso il vicino}
        hole_dir_used = {}

        # --- Posizionamento atomi ---
        for idx, sym, _pos_unused in atoms:
            pos    = coords[idx]
            neighs = [t for s, t in bonds if s == idx] + [s for s, t in bonds if t == idx]
            nn     = len(neighs)

            # choose_geometry_key può restituire qualcosa tipo "sp2", "sp3", "H", ecc.
            # Qui lo mappiamo sempre a una collection Atom_*
            base_key = choose_geometry_key(sym, nn)
            if base_key.startswith("Atom_"):
                key = base_key
            else:
                key = f"Atom_{base_key}"

            # Fallback robusto: se la collection non esiste, prova qualche default
            if key not in GEN_COLLECTIONS:
                # Carbonio aromatico/tetraedrico → usa almeno Atom_sp3
                if sym == "C":
                    key = "Atom_sp2" if "Atom_sp2" in GEN_COLLECTIONS else "Atom_sp3"
                elif sym == "H":
                    key = "Atom_H" if "Atom_H" in GEN_COLLECTIONS else "Atom_sp3"
                else:
                    key = "Atom_sp3" if "Atom_sp3" in GEN_COLLECTIONS else GEN_COLLECTIONS[0]

            if key not in hole_cache:
                try:
                    hole_cache[key] = load_hole_dirs(P.lib_path, key)
                except Exception as e:
                    print(f"[Molymod] ERROR loading hole dirs for '{key}':", e)
                    # Nessun hole: l'atomo verrà comunque istanziato, ma senza orientamento avanzato
                    hole_cache[key] = []

            hole_vecs = hole_cache[key]

            # Direzioni dei legami rispetto a questo atomo
            bond_dirs = []
            for n in neighs[:max(1, len(hole_vecs))]:
                v = coords[n] - pos
                if v.length > 1e-9:
                    bond_dirs.append(v.normalized())

            # Istanzia la collection
            try:
                bpy.ops.object.collection_instance_add(collection=key, location=(0, 0, 0))
            except Exception as e:
                print(f"[Molymod] collection_instance_add failed for '{key}':", e)
                continue

            inst = context.object
            inst.name = f"mol_{sym}_{idx}"

            # --- Orientamento istanza ---
            if nn == 1 and (sym == "H" or sym in HALOGENS) and len(bond_dirs) == 1:
                # H / alogeni: orienta asse scelto verso il legame e aggiungi roll opzionale
                b = bond_dirs[0]
                forward_local = axis_vec(P.H_forward_axis)
                q_align = forward_local.rotation_difference(b)
                q_roll  = Quaternion(b, radians(P.H_roll_deg))
                inst.matrix_world = (q_roll @ q_align).to_matrix().to_4x4()

            elif len(hole_vecs) >= 2 and len(bond_dirs) >= 2:
                # Geometrie con almeno due fori: usa Kabsch + Hungarian per la miglior corrispondenza
                cost   = [[1.0 - max(-1.0, min(1.0, h.dot(b))) for b in bond_dirs] for h in hole_vecs]
                rows, cols = hungarian_assign(cost)
                from_v = [hole_vecs[i] for i in rows]
                to_v   = [bond_dirs[j] for j in cols]
                R = kabsch_rotation(from_v, to_v)
                inst.matrix_world = R
                # Memorizza corrispondenza foro→vicino dopo la rotazione
                for r, c in zip(rows, cols):
                    n_idx = neighs[c]
                    rotated_hole = (R.to_3x3() @ hole_vecs[r]).normalized()
                    hole_dir_used[(idx, n_idx)] = rotated_hole

            elif len(bond_dirs) == 1 and len(hole_vecs) >= 1:
                # Un solo legame e almeno un foro: allinea il foro più vicino alla direzione del legame
                b = bond_dirs[0]
                h = max(hole_vecs, key=lambda v: v.dot(b))
                if h.dot(b) < 0.0:
                    h = -h
                R = align_one_vector(h, b)
                inst.matrix_world = R
                n_idx = neighs[0]
                rotated_hole = (R.to_3x3() @ h).normalized()
                hole_dir_used[(idx, n_idx)] = rotated_hole

            elif len(hole_vecs) >= 1 and len(bond_dirs) >= 1:
                # Fallback: usa il primo foro e il primo legame
                R = align_one_vector(hole_vecs[0], bond_dirs[0])
                inst.matrix_world = R
                n_idx = neighs[0]
                rotated_hole = (R.to_3x3() @ hole_vecs[0]).normalized()
                hole_dir_used[(idx, n_idx)] = rotated_hole

            else:
                # Nessuna informazione direzionale → identità
                inst.matrix_world = Matrix.Identity(4)

            inst.location = pos

            # Rendi reali le istanze e tieni un solo oggetto per atomo
            bpy.ops.object.select_all(action='DESELECT')
            inst.select_set(True)
            bpy.context.view_layer.objects.active = inst
            try:
                bpy.ops.object.duplicates_make_real(
                    use_hierarchy=True,
                    use_base_parent=False,
                    use_keep_transform=True,
                )
                new_objs = [o for o in bpy.context.selected_objects if o != inst]
            except Exception as e:
                print("[Instances] duplicates_make_real failed:", e)
                new_objs = []

            if new_objs:
                for o in new_objs:
                    o.select_set(False)
                try:
                    bpy.data.objects.remove(inst, do_unlink=True)
                except:
                    pass
                placed[idx] = new_objs[0]
            else:
                placed[idx] = inst

        # Cap template/material (usato per i cappucci alle estremita')
        cap_template = None
        cap_mat      = None
        if P.use_caps:
            cap_template = _load_cap_template(P)
            cap_mat      = _get_or_make_cap_material(P)

        bond_r = P.bond_radius * P.scale / 3.0
        tf     = P.bond_tangent_factor  # 0..0.49, default 0.35

        # --- Disegno legami con spline Bezier ---
        for s, t in bonds:
            if s not in placed or t not in placed:
                continue
            if not placed[s] or not placed[t]:
                continue

            pos_s = coords[s]
            pos_t = coords[t]
            vec   = pos_t - pos_s
            dist  = vec.length
            if dist <= 1e-9:
                continue

            # Offset di gap dal centro dell'atomo
            offA  = max(0.0, P.bond_gap_each_side) + P.bond_start_offset
            offB  = max(0.0, P.bond_gap_each_side) + P.bond_end_offset
            dirn  = vec.normalized()
            p0    = pos_s + dirn * offA   # base cap lato s
            p3    = pos_t - dirn * offB   # base cap lato t

            # Direzione dei fori (gia' ruotata con l'atomo)
            # Se non disponibile, fallback alla direzione del legame
            dir_s = hole_dir_used.get((s, t), dirn)
            dir_t = hole_dir_used.get((t, s), -dirn)

            # tlen: lunghezza tangente = frazione della distanza effettiva p0-p3
            eff_dist = (p3 - p0).length
            tlen     = eff_dist * tf

            # P1: maniglia di uscita da p0 nella direzione del foro di s
            # P2: maniglia di entrata su p3 nella direzione del foro di t
            p1 = p0 + dir_s * tlen
            p2 = p3 + dir_t * tlen

            bond_name = f"bond_{s}_{t}"
            _make_bezier_bond(p0, p1, p2, p3, bond_r, bond_name, context)

            if P.use_caps and cap_mat:
                _add_cap_at(p0,  dir_s,  P, cap_mat, cap_template)
                _add_cap_at(p3,  dir_t,  P, cap_mat, cap_template)

        self.report({'INFO'}, "Molymod build complete (Bezier bonds) ✅")
        return {'FINISHED'}


class MOLYMOD_OT_ValidateLibrary(bpy.types.Operator):
    bl_idname = "molymod.validate_library"
    bl_label = "Validate Library"
    bl_options = {'REGISTER',}

    def execute(self, context):
        P   = context.scene.molymod_settings
        lib = abspath(P.lib_path)
        if not os.path.isfile(lib):
            self.report({'ERROR'}, f"Library .blend not found: {lib}")
            return {'CANCELLED'}
        missing = []
        with bpy.data.libraries.load(lib, link=False) as (src, dst):
            src_colls = set(src.collections)
            src_objs  = set(src.objects)
        for ck in GEN_COLLECTIONS:
            if ck not in src_colls:
                missing.append(f"[Collection] {ck}")
        for ck in GEN_COLLECTIONS:
            if not any(name.startswith(f"{ck}_hole") for name in src_objs):
                missing.append(f"[Holes] {ck}_hole#")
        msg = ("Missing in library:\n- " + "\n- ".join(missing)) if missing else \
              "Library looks good: all Atom_sp* collections and holes found."
        self.report({'INFO'}, msg)
        print(msg)
        return {'FINISHED'}


class MOLYMOD_OT_ClearAll(bpy.types.Operator):
    bl_idname = "molymod.clear_all"
    bl_label = "Clear All Molecules"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        count = 0
        for o in list(bpy.data.objects):
            if o.name.startswith(("mol_", "bond_", "bond_cap", "DEBUG_")):
                bpy.data.objects.remove(o, do_unlink=True)
                count += 1
        self.report({'INFO'}, f"Removed {count} objects.")
        return {'FINISHED'}
