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
    radius = raggio del tubo
    """
    curve_data = bpy.data.curves.new(name, type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.resolution_u = 12
    curve_data.bevel_depth = radius
    curve_data.bevel_resolution = 4
    curve_data.use_fill_caps = True
    spline = curve_data.splines.new('BEZIER')
    spline.bezier_points.add(1)
    bp0 = spline.bezier_points[0]
    bp0.co           = p0
    bp0.handle_left  = p0
    bp0.handle_right = p1
    bp0.handle_left_type  = 'FREE'
    bp0.handle_right_type = 'FREE'
    bp1 = spline.bezier_points[1]
    bp1.co           = p3
    bp1.handle_left  = p2
    bp1.handle_right = p3
    bp1.handle_left_type  = 'FREE'
    bp1.handle_right_type = 'FREE'
    obj = bpy.data.objects.new(name, curve_data)
    context.scene.collection.objects.link(obj)
    return obj


def _pick_hole(hole_vecs_world, target_dir):
    """Dalla lista di vettori foro in world-space,
    restituisce quello che punta meglio nella direzione target_dir.
    """
    if not hole_vecs_world:
        return target_dir.normalized()
    best = max(hole_vecs_world, key=lambda h: h.dot(target_dir))
    if best.dot(target_dir) < 0.0:
        best = -best
    return best.normalized()


class MOLYMOD_OT_Build(bpy.types.Operator):
    bl_idname  = "molymod.build"
    bl_label   = "Build Molecule from File"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        P       = context.scene.molymod_settings
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

        for key in GEN_COLLECTIONS:
            try:
                append_collection(P.lib_path, key)
            except Exception as e:
                print(f"[Molymod] '{key}' not found (ok if unused). {e}")

        atoms, bonds, coords, types, bond_orders = parse_atoms_bonds(molfile, P.scale)

        if abs(P.compact_factor - 1.0) > 1e-9:
            for k in coords:
                coords[k] = coords[k] * P.compact_factor

        hole_cache    = {}
        placed        = {}
        all_hole_dirs = {}

        for idx, sym, _pos_unused in atoms:
            pos    = coords[idx]
            neighs = [t for s, t in bonds if s == idx] + \
                     [s for s, t in bonds if t == idx]
            nn = len(neighs)

            base_key = choose_geometry_key(sym, nn)
            key = base_key if base_key.startswith("Atom_") else f"Atom_{base_key}"

            if key not in GEN_COLLECTIONS:
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
                    hole_cache[key] = []

            hole_vecs = hole_cache[key]

            bond_dirs = []
            for n in neighs[:max(1, len(hole_vecs))]:
                v = coords[n] - pos
                if v.length > 1e-9:
                    bond_dirs.append(v.normalized())

            try:
                bpy.ops.object.collection_instance_add(collection=key, location=(0, 0, 0))
            except Exception as e:
                print(f"[Molymod] collection_instance_add failed for '{key}':", e)
                continue

            inst = context.object
            inst.name = f"mol_{sym}_{idx}"

            R = Matrix.Identity(4)
            if nn == 1 and (sym == "H" or sym in HALOGENS) and len(bond_dirs) == 1:
                b = bond_dirs[0]
                forward_local = axis_vec(P.H_forward_axis)
                q_align = forward_local.rotation_difference(b)
                q_roll  = Quaternion(b, radians(P.H_roll_deg))
                R = (q_roll @ q_align).to_matrix().to_4x4()
            elif len(hole_vecs) >= 2 and len(bond_dirs) >= 2:
                cost = [[1.0 - max(-1.0, min(1.0, h.dot(b)))
                         for b in bond_dirs] for h in hole_vecs]
                rows, cols = hungarian_assign(cost)
                from_v = [hole_vecs[i] for i in rows]
                to_v   = [bond_dirs[j] for j in cols]
                R = kabsch_rotation(from_v, to_v)
            elif len(bond_dirs) == 1 and len(hole_vecs) >= 1:
                b = bond_dirs[0]
                h = max(hole_vecs, key=lambda v: v.dot(b))
                if h.dot(b) < 0.0:
                    h = -h
                R = align_one_vector(h, b)
            elif len(hole_vecs) >= 1 and len(bond_dirs) >= 1:
                R = align_one_vector(hole_vecs[0], bond_dirs[0])

            inst.matrix_world = R
            inst.location = pos

            R3 = R.to_3x3()
            all_hole_dirs[idx] = [(R3 @ h).normalized() for h in hole_vecs]

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

        cap_template = None
        cap_mat      = None
        if P.use_caps:
            cap_template = _load_cap_template(P)
            cap_mat      = _get_or_make_cap_material(P)

        bond_r = P.bond_radius * P.scale / 3.0
        tf     = P.bond_tangent_factor

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

            offA = max(0.0, P.bond_gap_each_side) + P.bond_start_offset
            offB = max(0.0, P.bond_gap_each_side) + P.bond_end_offset
            dirn = vec.normalized()
            p0   = pos_s + dirn * offA
            p3   = pos_t - dirn * offB

            eff_dist = (p3 - p0).length
            tlen     = eff_dist * tf

            h_s = _pick_hole(all_hole_dirs.get(s, []), dirn)
            h_t = _pick_hole(all_hole_dirs.get(t, []), -dirn)

            p1 = p0 + h_s * tlen
            p2 = p3 + h_t * tlen

            bond_name = f"bond_{s}_{t}"
            _make_bezier_bond(p0, p1, p2, p3, bond_r, bond_name, context)

            if P.use_caps and cap_mat:
                _add_cap_at(p0, h_s, P, cap_mat, cap_template)
                _add_cap_at(p3, h_t, P, cap_mat, cap_template)

        self.report({'INFO'}, "Molymod build complete (hole-guided Bezier bonds)")
        return {'FINISHED'}


class MOLYMOD_OT_ValidateLibrary(bpy.types.Operator):
    bl_idname  = "molymod.validate_library"
    bl_label   = "Validate Library"
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
    bl_idname  = "molymod.clear_all"
    bl_label   = "Clear All Molecules"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        count = 0
        for o in list(bpy.data.objects):
            if o.name.startswith(("mol_", "bond_", "bond_cap", "DEBUG_")):
                bpy.data.objects.remove(o, do_unlink=True)
                count += 1
        self.report({'INFO'}, f"Removed {count} objects.")
        return {'FINISHED'}
