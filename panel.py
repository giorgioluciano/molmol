import bpy

# UI Panel definition for the Molymod Builder addon
class MOLYMOD_PT_Panel(bpy.types.Panel):
    """UI panel in the N-panel for molecule builder settings and actions."""
    bl_idname = "MOLYMOD_PT_panel"
    bl_label = "Molymod"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Molymod"

    def draw(self, context):
        # Access settings property group
        P = context.scene.molymod_settings
        layout = self.layout

        # Main parameters
        col = layout.column(align=True)
        col.prop(P, "lib_path")
        col.prop(P, "molecule_path")
        col.prop(P, "scale")
        col.prop(P, "compact_factor")
        col.prop(P, "clear_previous")
        col.prop(P, "debug_mode")
        layout.separator()

        # Bond visualization box
        bond_box = layout.box()
        bond_box.label(text="Bonds")
        bond_box.prop(P, "bond_radius")
        bond_box.prop(P, "bond_gap_each_side")
        rowb = bond_box.row(align=True)
        rowb.prop(P, "bond_start_offset")
        rowb.prop(P, "bond_end_offset")
        bond_box.prop(P, "bond_length_factor")
        bond_box.prop(P, "bond_vertices")
        bond_box.prop(P, "bond_mat_name")

        layout.separator()
        # Optional caps
        cap_box = layout.box(); cap_box.label(text="Caps")
        cap_box.prop(P, "use_caps")
        cap_box.prop(P, "cap_template_name")
        rowcs = cap_box.row(align=True); rowcs.prop(P, "cap_scale")
        rowcb = cap_box.row(align=True); rowcb.prop(P, "cap_radius"); rowcb.prop(P, "cap_length")
        rowca = cap_box.row(align=True); rowca.prop(P, "cap_forward_axis"); rowca.prop(P, "cap_roll_deg")
        cap_box.prop(P, "cap_offset")
        rowc2 = cap_box.row(align=True); rowc2.prop(P, "cap_start_offset"); rowc2.prop(P, "cap_end_offset")
        cap_box.prop(P, "cap_mat_name")
        layout.separator()

        # Hydrogen / monovalent box
        boxh = layout.box(); boxh.label(text="Hydrogen / Monovalent")
        rowh = boxh.row(align=True); rowh.prop(P, "H_forward_axis"); rowh.prop(P, "H_roll_deg")
        layout.separator()

        # Main action buttons
        row = layout.row(align=True)
        row.operator("molymod.build", text="Build Molecule", icon="MESH_CYLINDER")
        row.operator("molymod.validate_library", text="Validate Library", icon="CHECKMARK")
        row.operator("molymod.clear_all", text="Clear All", icon="TRASH")
