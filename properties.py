import bpy

GEN_COLLECTIONS = ["Atom_sp3", "Atom_sp2", "Atom_sp", "Atom_bent", "Atom_sp3d2"]
HALOGENS = {"F", "Cl", "Br", "I"}

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
    bond_gap_each_side: bpy.props.FloatProperty(
        name="Gap Each Side", default=0.2, min=0.0, soft_max=3.0,
        description="DEPRECATED: Not used in new version (atom_radius auto-calculated)"
    )
    bond_start_offset: bpy.props.FloatProperty(
        name="Start Offset", default=0.0, soft_min=-2.0, soft_max=2.0,
        description="DEPRECATED: Not used in new version"
    )
    bond_end_offset: bpy.props.FloatProperty(
        name="End Offset", default=0.0, soft_min=-2.0, soft_max=2.0,
        description="DEPRECATED: Not used in new version"
    )
    bond_length_factor: bpy.props.FloatProperty(
        name="Length Factor", default=1.0, min=0.1, max=1.5,
        description="DEPRECATED: Not used in new version"
    )
    bond_vertices: bpy.props.IntProperty(name="Vertices", default=30, min=3, soft_max=128)
    bond_mat_name: bpy.props.StringProperty(name="Bond Material", default="MolBond")

    bond_tangent_factor: bpy.props.FloatProperty(
        name="Bond Tangent Factor",
        description="Fattore tangente per la spline Bezier dei legami (0=retta, <0.5=curva). tlen = dist * factor",
        default=0.35, min=0.0, max=0.49, soft_min=0.1, soft_max=0.45
    )
    
    # Caps
    use_caps: bpy.props.BoolProperty(name="Add Caps on Bonds", default=False)
    cap_template_name: bpy.props.StringProperty(name="Cap Template", default="")
    cap_scale: bpy.props.FloatProperty(name="Cap Scale", default=1.0, min=0.01, soft_max=10.0)
    cap_radius: bpy.props.FloatProperty(name="Base Radius", default=0.55, min=0.001, soft_max=1.0)
    cap_length: bpy.props.FloatProperty(name="Base Length", default=0.50, min=0.001, soft_max=2.0)
    cap_roll_deg: bpy.props.FloatProperty(name="Cap Roll (deg)", default=0.0, soft_min=-180.0, soft_max=180.0)
    cap_forward_axis: bpy.props.EnumProperty(
        name="Cap Forward Axis",
        items=[('Z+','Z+',''),('Z-','Z-',''),('X+','X+',''),('X-','X-',''),('Y+','Y+',''),('Y-','Y-','')],
        default='Z+'
    )
    cap_offset: bpy.props.FloatProperty(
        name="Cap Offset", default=0.10, min=0.0, soft_max=2.0,
        description="DEPRECATED: Not used in new version"
    )
    cap_start_offset: bpy.props.FloatProperty(
        name="Start Cap Offset", default=0.0, soft_min=-2.0, soft_max=2.0,
        description="DEPRECATED: Not used in new version"
    )
    cap_end_offset: bpy.props.FloatProperty(
        name="End Cap Offset", default=0.0, soft_min=-2.0, soft_max=2.0,
        description="DEPRECATED: Not used in new version"
    )
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

