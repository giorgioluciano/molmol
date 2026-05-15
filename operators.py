import bpy, os
from math import radians
from mathutils import Vector, Matrix, Quaternion

from .helpers import (
    GEN_COLLECTIONS, HALOGENS, MIN_VALENCE,
    abspath, ensure_hidden_bucket, unlink_collection_everywhere,
    append_collection, append_object, load_hole_dirs,
    kabsch_rotation, align_one_vector, hungarian_assign,
    parse_atoms_bonds, choose_geometry_key, axis_vec,
    _load_cap_template, _get_or_make_cap_material, _add_cap_at,
    assign_double_bond_holes, detect_missing_hydrogens, ATOM_VDW_RADII
)

def _make_bezier_bond(p0, p1, p2, p3, radius, name, context):
    """Crea un legame Bezier cubico con bevel (tubo 3D).
    p0, p3 = punti sulla superficie degli atomi (dai fori)
    p1, p2 = maniglie tangenti nella direzione dei fori
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
    bp0.co = p0
    bp0.handle_left = p0
    bp0.handle_right = p1
    bp0.handle_left_type = 'FREE'
    bp0.handle_right_type = 'FREE'

    bp1 = spline.bezier_points[1]
    bp1.co = p3
    bp1.handle_left = p2
    bp1.handle_right = p3
    bp1.handle_left_type = 'FREE'
    bp1.handle_right_type = 'FREE'

    obj = bpy.data.objects.new(name, curve_data)
    context.scene.collection.objects.link(obj)
    return obj

def _create_straight_cylinder(p0, p3, radius, name, context):
    """Crea un cilindro dritto da p0 a p3."""
    length = (p3 - p0).length
    if length < 1e-9:
        return None
    center = (p0 + p3) / 2
    direction = (p3 - p0).normalized()
    z_axis = Vector((0, 0, 1))
    q = z_axis.rotation_difference(direction)
    bpy.ops.mesh.primitive_cylinder_add(
        radius=radius,
        depth=length,
        location=center,
        rotation=q.to_euler()
    )
    cyl = context.object
    cyl.name = name
    return cyl

def _pick_hole(hole_vecs_world, target_dir):
    """Restituisce il foro più allineato con target_dir."""
    if not hole_vecs_world:
        return target_dir.normalized()
    best = max(hole_vecs_world, key=lambda h: h.dot(target_dir))
    if best.dot(target_dir) < 0.0:
        best = -best
    return best.normalized()

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

        for key in GEN_COLLECTIONS:
            try:
                append_collection(P.lib_path, key)
            except Exception as e:
                print(f"[Molymod] '{key}' not found in library (ok if unused). {e}")

        # PARSE FILE
        result = parse_atoms_bonds(molfile, P.scale)
        if result is None:
            self.report({'ERROR'}, "Failed to parse molecule file. Check console for details.")
            return {'CANCELLED'}

        atoms, bonds, coords, types, bond_orders = result

        if abs(P.compact_factor - 1.0) > 1e-9:
            for k in coords:
                coords[k] = coords[k] * P.compact_factor

        # Rileva H mancanti
        missing_H = detect_missing_hydrogens(atoms, bonds, bond_orders, types)
        if missing_H:
            print(f"[INFO] Total missing H: {sum(missing_H.values())}")

        hole_cache = {}
        placed = {}
        all_hole_dirs = {}
        hole_assignments = {}
        atom_radii = {}
        used_holes = {}

        # ============ LOOP ATOMI ============
        print("\n[ATOMS] Placing and orienting atoms...")
        for idx, sym, _pos_unused in atoms:
            pos = coords[idx]
            neighs = [t for s, t in bonds if s == idx] + \
                     [s for s, t in bonds if t == idx]
            nn = len(neighs)

            # Calcola fori necessari
            nn_effective = 0
            for n in neighs:
                k = (min(idx, n), max(idx, n))
                order = bond_orders.get(k, 1)
                nn_effective += round(order)
            min_valence = MIN_VALENCE.get(sym, nn)
            nn_effective = max(nn_effective, min_valence)

            base_key = choose_geometry_key(sym, nn_effective)
            key = base_key if base_key.startswith("Atom_") else f"Atom_{base_key}"

            # Scegli collection giusta
            if sym == "H":
                if "H_1_90" in GEN_COLLECTIONS:
                    key = "H_1_90"
                elif "Atom_H" in GEN_COLLECTIONS:
                    key = "Atom_H"
                else:
                    key = "Atom_sp"
            elif key not in GEN_COLLECTIONS:
                if sym == "C":
                    key = "Atom_sp2" if "Atom_sp2" in GEN_COLLECTIONS else "Atom_sp3"
                else:
                    key = "Atom_sp3" if "Atom_sp3" in GEN_COLLECTIONS else GEN_COLLECTIONS[0]

            if key not in hole_cache:
                try:
                    hole_cache[key] = load_hole_dirs(P.lib_path, key)
                except Exception as e:
                    print(f"[ERROR] Loading hole dirs for '{key}': {e}")
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
                print(f"[ERROR] collection_instance_add failed for '{key}': {e}")
                continue

            inst = context.object
            inst.name = f"mol_{sym}_{idx}"

            # COORDINATE FISSE
            inst.location = pos

            # Calcola rotazione
            R = Matrix.Identity(4)
            if nn == 1 and (sym == "H" or sym in HALOGENS) and len(bond_dirs) == 1:
                b = bond_dirs[0]
                forward_local = axis_vec(P.H_forward_axis)
                q_align = forward_local.rotation_difference(b)
                q_roll = Quaternion(b, radians(P.H_roll_deg))
                R = (q_roll @ q_align).to_matrix().to_4x4()
            elif len(hole_vecs) >= 2 and len(bond_dirs) >= 2:
                cost = [[1.0 - max(-1.0, min(1.0, h.dot(b))) for b in bond_dirs]
                        for h in hole_vecs]
                rows, cols = hungarian_assign(cost)
                from_v = [hole_vecs[i] for i in rows]
                to_v = [bond_dirs[j] for j in cols]
                R = kabsch_rotation(from_v, to_v)
            elif len(bond_dirs) == 1 and len(hole_vecs) >= 1:
                b = bond_dirs[0]
                h = max(hole_vecs, key=lambda v: v.dot(b))
                if h.dot(b) < 0.0:
                    h = -h
                R = align_one_vector(h, b)
            elif len(hole_vecs) >= 1 and len(bond_dirs) >= 1:
                R = align_one_vector(hole_vecs[0], bond_dirs[0])

            # Applica matrice (rotazione + posizione)
            inst.matrix_world = R
            inst.location = pos

            R3 = R.to_3x3()
            all_hole_dirs[idx] = [(R3 @ h).normalized() for h in hole_vecs]

            # Assegna fori per legami multipli
            has_multiple_bonds = any(
                bond_orders.get((min(idx, n), max(idx, n)), 1) >= 2
                for n in neighs
            )
            if has_multiple_bonds:
                neighbors_with_orders = [
                    (n, (coords[n] - pos).normalized(),
                     bond_orders.get((min(idx,n), max(idx,n)), 1))
                    for n in neighs if (coords[n] - pos).length > 1e-9
                ]
                hole_assignments[idx] = assign_double_bond_holes(
                    all_hole_dirs[idx], neighbors_with_orders
                )

            bpy.ops.object.select_all(action='DESELECT')
            inst.select_set(True)
            bpy.context.view_layer.objects.active = inst

            try:
                bpy.ops.object.duplicates_make_real(
                    use_hierarchy=True,
                    use_base_parent=False,
                )
                new_objs = [o for o in bpy.context.selected_objects if o != inst]
            except Exception as e:
                print(f"[WARNING] duplicates_make_real failed: {e}")
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

            # Forza posizione corretta
            placed[idx].location = pos

            # Calcola radius dalla bounding box
            dims = placed[idx].dimensions
            atom_radii[idx] = max(dims) / 2.0 if max(dims) > 0 else 0.5

            # Inizializza tracking fori
            used_holes[idx] = []

            if P.debug_mode:
                print(f"[DEBUG] Atom {idx} ({sym}): key={key} radius={atom_radii[idx]:.3f}")

        print(f"[OK] Placed {len(placed)} atoms")

        # ============ LOOP BOND ============
        print("\n[BONDS] Drawing bonds...")
        cap_template = None
        cap_mat = None
        if P.use_caps:
            cap_template = _load_cap_template(P)
            cap_mat = _get_or_make_cap_material(P)

        bond_r = P.bond_radius * P.scale / 3.0
        tf = P.bond_tangent_factor
        bonds_drawn = 0

        # Separa singoli e multipli
        single_bonds = [(s,t) for s,t in bonds
                        if bond_orders.get((s,t), bond_orders.get((t,s), 1)) == 1]
        multiple_bonds = [(s,t) for s,t in bonds
                          if bond_orders.get((s,t), bond_orders.get((t,s), 1)) > 1]

        print(f"[INFO] Single bonds: {len(single_bonds)}, Multiple bonds: {len(multiple_bonds)}")

        # ========== PASSATA 1: SINGOLI ==========
        print("[BONDS] Pass 1: Single bonds...")
        for s, t in single_bonds:
            if s not in placed or t not in placed:
                continue

            pos_s = coords[s]
            pos_t = coords[t]
            vec = pos_t - pos_s
            dist = vec.length
            if dist <= 1e-9:
                continue

            dirn = vec.normalized()
            radius_s = atom_radii.get(s, 0.5)
            radius_t = atom_radii.get(t, 0.5)

            holes_s_all = all_hole_dirs.get(s, [])
            holes_t_all = all_hole_dirs.get(t, [])

            free_holes_s = [h for i,h in enumerate(holes_s_all)
                            if i not in used_holes.get(s, [])]
            free_holes_t = [h for i,h in enumerate(holes_t_all)
                            if i not in used_holes.get(t, [])]

            if not free_holes_s or not free_holes_t:
                print(f"[WARNING] No free holes for single bond {s}-{t}!")
                continue

            h_s = max(free_holes_s, key=lambda h: h.dot(dirn))
            h_t_orig = max(free_holes_t, key=lambda h: h.dot(-dirn))

            # Salva indici PRIMA di modificare
            idx_s = holes_s_all.index(h_s)
            idx_t = holes_t_all.index(h_t_orig)
            used_holes[s].append(idx_s)
            used_holes[t].append(idx_t)

            # Orienta h_t nella direzione giusta
            h_t = h_t_orig if h_t_orig.dot(-dirn) >= 0 else -h_t_orig

            # Partenza dalla superficie
            p0 = pos_s + h_s * radius_s
            p3 = pos_t + h_t * radius_t

            #            # Dritto se fori opposti
            if h_s.dot(h_t) < -0.7:
                _create_straight_cylinder(p0, p3, bond_r, f"bond_{s}_{t}", context)
            else:
                eff_dist = (p3 - p0).length
                tlen = eff_dist * tf
                p1 = p0 + h_s * tlen
                p2 = p3 + h_t * tlen
                _make_bezier_bond(p0, p1, p2, p3, bond_r, f"bond_{s}_{t}", context)

            if P.use_caps and cap_mat:
                _add_cap_at(p0, h_s, P, cap_mat, cap_template)
                _add_cap_at(p3, h_t, P, cap_mat, cap_template)

            bonds_drawn += 1

        # ========== PASSATA 2: MULTIPLI ==========
        print("[BONDS] Pass 2: Multiple bonds...")
        for s, t in multiple_bonds:
            if s not in placed or t not in placed:
                continue

            pos_s = coords[s]
            pos_t = coords[t]
            vec = pos_t - pos_s
            dist = vec.length
            if dist <= 1e-9:
                continue

            dirn = vec.normalized()
            radius_s = atom_radii.get(s, 0.5)
            radius_t = atom_radii.get(t, 0.5)

            order = bond_orders.get((s,t), bond_orders.get((t,s), 1))
            n_tubes = round(order)

            holes_s_all = all_hole_dirs.get(s, [])
            holes_t_all = all_hole_dirs.get(t, [])

            free_holes_s = [h for i,h in enumerate(holes_s_all)
                            if i not in used_holes.get(s, [])]
            free_holes_t = [h for i,h in enumerate(holes_t_all)
                            if i not in used_holes.get(t, [])]

            if len(free_holes_s) < n_tubes or len(free_holes_t) < n_tubes:
                print(f"[WARNING] Not enough free holes for bond {s}-{t} "
                      f"(need {n_tubes}, have {len(free_holes_s)}/{len(free_holes_t)})")
                n_tubes = min(len(free_holes_s), len(free_holes_t), n_tubes)

            if n_tubes == 0:
                continue

            # Selezione fori con garanzia di complanarità
            from .helpers import find_coplanar_holes
            
            holes_s_selected, holes_t_selected = find_coplanar_holes(
                free_holes_s, free_holes_t, dirn, n_tubes
            )
            
            # Aggiorna sorted per il resto del codice
            holes_s_sorted = holes_s_selected
            holes_t_sorted = holes_t_selected

                if dot_crossed > dot_direct:
                    holes_t_sorted = [holes_t_sorted[1], holes_t_sorted[0]]
                    print(f"[INFO] Bond {s}-{t}: swapped holes to avoid crossing")

            # Marca fori usati
            for h in holes_s_sorted:
                used_holes[s].append(holes_s_all.index(h))
            for h in holes_t_sorted:
                used_holes[t].append(holes_t_all.index(h))

            eff_dist = dist - radius_s - radius_t
            tlen = eff_dist * tf

            for tube_i in range(n_tubes):
                h_s = holes_s_sorted[tube_i].normalized()
                h_t = holes_t_sorted[tube_i].normalized()

                p0 = pos_s + h_s * radius_s
                p3 = pos_t + h_t * radius_t
                p1 = p0 + h_s * tlen
                p2 = p3 + h_t * tlen

                tube_name = f"bond_{s}_{t}_{tube_i}"
                tube_r = bond_r * 0.85
                _make_bezier_bond(p0, p1, p2, p3, tube_r, tube_name, context)

                if P.use_caps and cap_mat:
                    _add_cap_at(p0, h_s, P, cap_mat, cap_template)
                    _add_cap_at(p3, h_t, P, cap_mat, cap_template)

            bonds_drawn += 1

        print(f"[OK] Drew {bonds_drawn} bonds")

        # ========== PASSATA 3: H MANCANTI ==========
        print("[BONDS] Pass 3: Missing hydrogens...")
        if missing_H:
            max_idx = max(idx for idx, _, _ in atoms)
            H_added = 0

            for parent_idx, n_H in missing_H.items():
                holes_all = all_hole_dirs.get(parent_idx, [])
                free_holes = [h for i,h in enumerate(holes_all)
                             if i not in used_holes.get(parent_idx, [])]

                if not free_holes:
                    print(f"[WARNING] No free holes for H on atom "
                          f"{parent_idx} ({types.get(parent_idx,'?')})!")
                    continue

                parent_pos = coords[parent_idx]
                parent_radius = atom_radii.get(parent_idx, 0.5)

                for i in range(min(n_H, len(free_holes))):
                    max_idx += 1
                    h_dir = free_holes[i].normalized()
                    H_pos = parent_pos + h_dir * 1.09 * P.scale

                    try:
                        h_key = "H_1_90" if "H_1_90" in GEN_COLLECTIONS else \
                                "Atom_H" if "Atom_H" in GEN_COLLECTIONS else "Atom_sp"

                        bpy.ops.object.collection_instance_add(
                            collection=h_key,
                            location=H_pos
                        )
                        h_obj = context.object
                        h_obj.name = f"mol_H_{max_idx}"

                        # Orienta H verso parent
                        h_forward = axis_vec(P.H_forward_axis)
                        h_dir_to_parent = (parent_pos - H_pos).normalized()
                        q_align = h_forward.rotation_difference(h_dir_to_parent)
                        h_obj.rotation_euler = q_align.to_euler()

                        bpy.ops.object.select_all(action='DESELECT')
                        h_obj.select_set(True)
                        bpy.context.view_layer.objects.active = h_obj
                        bpy.ops.object.duplicates_make_real(
                            use_hierarchy=True,
                            use_base_parent=False,
                        )
                        new_h_objs = [o for o in bpy.context.selected_objects
                                      if o != h_obj]
                        if new_h_objs:
                            for o in new_h_objs:
                                o.select_set(False)
                            bpy.data.objects.remove(h_obj, do_unlink=True)

                        # Legame dritto H
                        p0_parent = parent_pos + h_dir * parent_radius
                        _create_straight_cylinder(
                            p0_parent, H_pos, bond_r * 0.7,
                            f"bond_{parent_idx}_{max_idx}", context
                        )

                        # Marca foro usato
                        idx_hole = holes_all.index(free_holes[i])
                        used_holes[parent_idx].append(idx_hole)
                        H_added += 1

                    except Exception as e:
                        print(f"[WARNING] Failed to add H to atom {parent_idx}: {e}")

            print(f"[OK] Added {H_added} hydrogens")

        # ========== DEBUG FORI ==========
        if P.debug_mode:
            print("\n[DEBUG] Hole usage summary:")
            for idx, used in used_holes.items():
                sym = types.get(idx, '?')
                total = len(all_hole_dirs.get(idx, []))
                print(f"  Atom {idx} ({sym}): {len(used)}/{total} holes used")

        self.report({'INFO'}, f"Molymod build complete: {len(placed)} atoms, {bonds_drawn} bonds")
        return {'FINISHED'}


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
            src_colls = set(src.collections)
            src_objs = set(src.objects)

        for ck in GEN_COLLECTIONS:
            if ck not in src_colls:
                missing.append(f"[Collection] {ck}")

        for ck in GEN_COLLECTIONS:
            if not any(name.startswith(f"{ck}_hole") for name in src_objs):
                missing.append(f"[Holes] {ck}_hole#")

        msg = ("Missing in library:\n- " + "\n- ".join(missing)) if missing else \
              "Library looks good!"
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
