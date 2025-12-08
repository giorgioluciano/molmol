import bpy, os
from math import radians
from mathutils import Vector, Matrix, Quaternion


from .helpers import (
    GEN_COLLECTIONS, HALOGENS,
    abspath, ensure_hidden_bucket, unlink_collection_everywhere,
    append_collection, append_object, load_hole_dirs,
    kabsch_rotation, align_one_vector, hungarian_assign,
    parse_atoms_bonds, choose_geometry_key,axis_vec,  _load_cap_template,
    _get_or_make_cap_material,_add_cap_at )# ASE-driven parser: supports CIF/PDB with bond inference)
# You may add additional helpers as needed.

# Main Operator: Build molecule
class MOLYMOD_OT_Build(bpy.types.Operator):
    bl_idname = "molymod.build"
    bl_label = "Build Molecule from File"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        P = context.scene.molymod_settings
        lib = abspath(P.lib_path)
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

        # Load all atom type template collections
        for key in GEN_COLLECTIONS:
            try:
                append_collection(P.lib_path, key)
            except Exception as e:
                print(f"[Molymod] '{key}' not found (ok if unused). {e}")

        # Use ASE parser to get atoms, bonds, coords, types (works for PDB/CIF and more)
        atoms, bonds, coords, types = parse_atoms_bonds(molfile, P.scale)

        # Optionally apply compact factor to coordinates
        if abs(P.compact_factor - 1.0) > 1e-9:
            for k in coords: coords[k] = coords[k] * P.compact_factor

        # Preload hole vectors for each atom type collection
        hole_cache = {}
        placed = {}

        # Instantiate atoms
        for idx, sym, _pos_unused in atoms:
            pos = coords[idx]
            # List of neighbors (bonded atom indices)
            neighs = [t for s, t in bonds if s == idx] + [s for s, t in bonds if t == idx]
            nn = len(neighs)
            key = choose_geometry_key(sym, nn)
            if key not in hole_cache:
                hole_cache[key] = load_hole_dirs(P.lib_path, key)
            hole_vecs = hole_cache[key]
            bond_dirs = []
            for n in neighs[:max(1, len(hole_vecs))]:
                v = coords[n] - pos
                if v.length > 1e-9:
                    bond_dirs.append(v.normalized())
            bpy.ops.object.collection_instance_add(collection=key, location=(0, 0, 0))
            inst = context.object; inst.name = f"mol_{sym}_{idx}"

            # Orientation logic (unchanged from working script)
            if nn == 1 and (sym == "H" or sym in HALOGENS) and len(bond_dirs) == 1:
                b = bond_dirs[0]
                forward_local = axis_vec(P.H_forward_axis)
                q_align = forward_local.rotation_difference(b)
                q_roll  = Quaternion(b, radians(P.H_roll_deg))
                inst.matrix_world = (q_roll @ q_align).to_matrix().to_4x4()
            elif len(hole_vecs) >= 2 and len(bond_dirs) >= 2:
                cost = [[1.0 - max(-1.0, min(1.0, h.dot(b))) for b in bond_dirs] for h in hole_vecs]
                rows, cols = hungarian_assign(cost)
                from_v = [hole_vecs[i] for i in rows]
                to_v   = [bond_dirs[j] for j in cols]
                inst.matrix_world = kabsch_rotation(from_v, to_v)
            elif len(bond_dirs) == 1 and len(hole_vecs) >= 1:
                b = bond_dirs[0]
                h = max(hole_vecs, key=lambda v: v.dot(b))
                if h.dot(b) < 0.0: h = -h
                inst.matrix_world = align_one_vector(h, b)
            elif len(hole_vecs) >= 1 and len(bond_dirs) >= 1:
                inst.matrix_world = align_one_vector(hole_vecs[0], bond_dirs[0])
            else:
                inst.matrix_world = Matrix.Identity(4)
            inst.location = pos

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
                try: bpy.data.objects.remove(inst, do_unlink=True)
                except: pass
                placed[idx] = new_objs[0]
            else:
                placed[idx] = inst

        # Optional cap template/material logic, as before
        cap_template = None
        cap_mat = None
        if P.use_caps:
            cap_template = _load_cap_template(P)
            cap_mat = _get_or_make_cap_material(P)
        if P.debug_mode and P.use_caps:
            print(f"[CAP] Using template: {cap_template}")

        # Draw bonds
        for s, t in bonds:
            if s in placed and t in placed and placed[s] and placed[t]:
                p1, p2 = coords[s], coords[t]
                vec = p2 - p1
                if vec.length <= 1e-9: continue
                dirn = vec.normalized()
                offA = max(0.0, P.bond_gap_each_side) + P.bond_start_offset
                offB = max(0.0, P.bond_gap_each_side) + P.bond_end_offset
                a = p1 + dirn * offA
                b = p2 - dirn * offB
                seg = b - a
                Leff = max(0.01, seg.length * P.bond_length_factor)
                mid = a + seg * 0.5
                rot = dirn.to_track_quat('Z','Y').to_euler()
                bpy.ops.mesh.primitive_cylinder_add(
                    vertices=P.bond_vertices,
                    radius=P.bond_radius * P.scale / 3.0,
                    depth=Leff,
                    location=mid,
                    rotation=rot
                )
                cyl = context.object; cyl.name = f"bond_{s}_{t}"
                # Material can be handled directly (no palette logic)
                # Add caps if enabled
                if P.use_caps:
                    a_cap = a + dirn * (P.cap_offset + P.cap_start_offset)
                    b_cap = b - dirn * (P.cap_offset + P.cap_end_offset)
                    _add_cap_at(a_cap,  dirn,  P, cap_mat, cap_template)
                    _add_cap_at(b_cap, -dirn,  P, cap_mat, cap_template)

        self.report({'INFO'}, "Molymod build complete ✅")
        return {'FINISHED'}

# Other operators, unchanged
class MOLYMOD_OT_ValidateLibrary(bpy.types.Operator):
    bl_idname = "molymod.validate_library"
    bl_label = "Validate Library"
    bl_options = {'REGISTER',}
    def execute(self, context):
        P = context.scene.molymod_settings
        lib = abspath(P.lib_path)
        if not os.path.isfile(lib):
            self.report({'ERROR'}, f"Library .blend not found: {lib}")
            return {'CANCELLED'}
        missing = []
        with bpy.data.libraries.load(lib, link=False) as (src, dst):
            src_colls = set(src.collections); src_objs = set(src.objects)
        for ck in GEN_COLLECTIONS:
            if ck not in src_colls: missing.append(f"[Collection] {ck}")
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

