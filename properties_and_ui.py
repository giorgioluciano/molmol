import bpy
from bpy.props import BoolProperty, FloatProperty, PointerProperty, StringProperty
from bpy.types import Operator, Panel, PropertyGroup

from .contact_pipeline import (
    ContactPipelineError,
    build_from_file_with_contact_pipeline,
    try_to_fix_independent_bonds_for_selected_objects,
)

from .fragment_recognizer import find_fragment_matches, load_molecule
from .smiles_builder import build_rdkit_balls_and_sticks


BUILT_COLLECTION_NAME = "MolMol_Built"


class MOLMOL_PG_Settings(PropertyGroup):
    structure_path: StringProperty(
        name="Structure File",
        description="Path to a .mol or .sdf file",
        default="",
        subtype="FILE_PATH",
    )

    library_path: StringProperty(
        name="Library File",
        description="Path to the MolMol library file",
        default="",
        subtype="FILE_PATH",
    )

    molecule_name: StringProperty(
        name="Name",
        description="Name for the generated molecule / collection",
        default="MolMol_FromFile",
    )

    atom_radius: FloatProperty(
        name="Fallback Radius",
        description="Fallback sphere radius when a template is missing",
        default=0.35,
        min=0.001,
    )

    bond_radius: FloatProperty(
        name="Bond Radius",
        description="Cylinder radius for ball-and-stick bonds",
        default=0.12,
        min=0.001,
    )

    coordinate_scale: FloatProperty(
        name="Coordinate Scale",
        description="Scale factor applied to coordinates read from the input file",
        default=1.0,
        min=0.0001,
    )

    add_hydrogens: BoolProperty(
        name="Add Hydrogens",
        description="Use RDKit to add explicit hydrogens before building the skeleton",
        default=False,
    )

    show_hydrogens: BoolProperty(
        name="Show Hydrogens",
        description="Render hydrogen atoms and hydrogen bonds in the skeleton",
        default=False,
    )

    center_rdkit: BoolProperty(
        name="Center Molecule",
        description="Center visible atoms around the origin before building",
        default=True,
    )

    layout_extra_offset: FloatProperty(
        name="Layout Extra Offset",
        description="Extra spacing factor applied to the initial template layout",
        default=0.10,
        min=0.0,
        max=0.5,
        soft_min=0.0,
        soft_max=0.2,
        precision=3,
    )

    show_labels: BoolProperty(
        name="Show Labels",
        description="Create debug labels for atoms",
        default=False,
    )

    clear_before_build: BoolProperty(
        name="Clear Before Build",
        description=f"Clear {BUILT_COLLECTION_NAME} before generating a new molecule",
        default=True,
    )


def ensure_built_collection():
    collection = bpy.data.collections.get(BUILT_COLLECTION_NAME)
    if collection is None:
        collection = bpy.data.collections.new(BUILT_COLLECTION_NAME)
        bpy.context.scene.collection.children.link(collection)
    return collection


def clear_collection_hierarchy(collection) -> int:
    if collection is None:
        return 0

    removed_object_names = set()

    def _clear_recursive(target_collection):
        for child_collection in list(target_collection.children):
            _clear_recursive(child_collection)
            bpy.data.collections.remove(child_collection)

        for obj in list(target_collection.objects):
            object_name = str(getattr(obj, "name", "") or "").strip()
            if object_name in removed_object_names:
                continue
            removed_object_names.add(object_name)
            bpy.data.objects.remove(obj, do_unlink=True)

    _clear_recursive(collection)
    return len(removed_object_names)


def clear_built_collection() -> int:
    built_collection = bpy.data.collections.get(BUILT_COLLECTION_NAME)
    if built_collection is None:
        return 0
    return clear_collection_hierarchy(built_collection)


def _validate_structure_path(settings):
    structure_path = (settings.structure_path or "").strip()
    if not structure_path:
        raise ValueError("Structure file path is empty")
    return structure_path


def _validate_library_path(settings):
    library_path = (settings.library_path or "").strip()
    if not library_path:
        raise ValueError("Library file path is empty")
    return library_path


class MOLMOL_OT_build_from_file(Operator):
    bl_idname = "molmol.build_from_file"
    bl_label = "Build Templates"
    bl_description = "Read a MOL/SDF file and build the molecule with Blender templates when available"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        settings = context.scene.molmol_settings

        try:
            structure_path = _validate_structure_path(settings)
            library_path = _validate_library_path(settings)

            built_collection = ensure_built_collection()

            if settings.clear_before_build:
                clear_built_collection()

            result = build_from_file_with_contact_pipeline(
                context=context,
                structure_path=structure_path,
                library_path=library_path,
                collection_name=settings.molecule_name or "MolMol_FromFile",
                extra_offset_factor=float(settings.layout_extra_offset),
            )

            built_child_collection = bpy.data.collections.get(result.built_collection_name)
            if built_child_collection is not None:
                already_linked = any(child == built_child_collection for child in built_collection.children)
                if not already_linked:
                    built_collection.children.link(built_child_collection)

        except ValueError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        except ContactPipelineError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        except Exception as exc:
            self.report({"ERROR"}, f"Build failed: {exc}")
            return {"CANCELLED"}

        resolved_count = len(result.match_result.resolved_edges)
        self.report({"INFO"}, f"Template build completed: {resolved_count} connections")
        return {"FINISHED"}


class MOLMOL_OT_build_rdkit_balls_sticks(Operator):
    bl_idname = "molmol.build_rdkit_balls_sticks"
    bl_label = "Build Balls+Sticks"
    bl_description = "Build the molecule directly as spheres and cylinders"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        settings = context.scene.molmol_settings

        try:
            structure_path = _validate_structure_path(settings)

            built_collection = ensure_built_collection()

            if settings.clear_before_build:
                clear_built_collection()

            result = build_rdkit_balls_and_sticks(
                context=context,
                structure_path=structure_path,
                molecule_name=(settings.molecule_name or "MolMol_RDKit") + "_RDKit",
                atom_radius=float(settings.atom_radius),
                bond_radius=float(settings.bond_radius),
                coordinate_scale=float(settings.coordinate_scale),
                add_hydrogens=bool(settings.add_hydrogens),
                show_hydrogens=bool(settings.show_hydrogens),
                center_molecule=bool(settings.center_rdkit),
                add_labels=bool(settings.show_labels),
            )

            created_collection_name = result.get("collection_name")
            if created_collection_name:
                created_collection = bpy.data.collections.get(created_collection_name)
                if created_collection is not None:
                    already_linked = any(child == created_collection for child in built_collection.children)
                    if not already_linked:
                        built_collection.children.link(created_collection)

        except ValueError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        except Exception as exc:
            self.report({"ERROR"}, f"RDKit build failed: {exc}")
            return {"CANCELLED"}

        self.report(
            {"INFO"},
            f"RDKit build completed: {result['atom_count']} atoms, {result['bond_count']} bonds",
        )
        return {"FINISHED"}


class MOLMOL_OT_inspect_fragments(Operator):
    bl_idname = "molmol.inspect_fragments"
    bl_label = "Inspect Fragments"
    bl_description = "Print recognized fragments to the system console and report the count"
    bl_options = {"REGISTER"}

    def execute(self, context):
        settings = context.scene.molmol_settings

        try:
            structure_path = _validate_structure_path(settings)
            mol = load_molecule(structure_path)
            matches = find_fragment_matches(mol)
        except ValueError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        except Exception as exc:
            self.report({"ERROR"}, f"Inspection failed: {exc}")
            return {"CANCELLED"}

        print("MOLMOL_FRAGMENT_MATCHES_BEGIN")
        for match in matches:
            print(match.fragment_id, match.atom_indices)
        print("MOLMOL_FRAGMENT_MATCHES_END")

        self.report({"INFO"}, f"Recognized {len(matches)} fragment matches")
        return {"FINISHED"}


class MOLMOL_OT_clear_built(Operator):
    bl_idname = "molmol.clear_built"
    bl_label = "Clear Built"
    bl_description = f"Clear the {BUILT_COLLECTION_NAME} collection"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        del context
        removed_count = clear_built_collection()
        self.report({"INFO"}, f"Cleared {BUILT_COLLECTION_NAME}: removed {removed_count} objects")
        return {"FINISHED"}


class MOLMOL_OT_try_to_fix_bonds(Operator):
    bl_idname = "molmol.try_to_fix_bonds"
    bl_label = "Try To Fix Selected Bonds"
    bl_description = "Try to rigidly reposition selected independent bonds using selected fragment holes"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        try:
            fixed_bonds = try_to_fix_independent_bonds_for_selected_objects(context)
        except ContactPipelineError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        except Exception as exc:
            self.report({"ERROR"}, f"Bond fix failed: {exc}")
            return {"CANCELLED"}

        if not fixed_bonds:
            self.report({"INFO"}, "No selected independent bonds were fixed")
            return {"FINISHED"}

        self.report({"INFO"}, f"Fixed {len(fixed_bonds)} selected independent bonds")
        return {"FINISHED"}


    

class MOLMOL_PT_panel(Panel):
    bl_label = "MolMol"
    bl_idname = "MOLMOL_PT_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "MolMol"

    def draw(self, context):
        layout = self.layout
        settings = context.scene.molmol_settings

        file_box = layout.box()
        file_box.label(text="Input", icon="FILE_FOLDER")
        file_col = file_box.column(align=True)
        file_col.prop(settings, "structure_path")
        file_col.prop(settings, "library_path")
        file_col.prop(settings, "molecule_name")

        template_box = layout.box()
        template_box.label(text="Template Build", icon="OUTLINER_COLLECTION")
        template_col = template_box.column(align=True)
        template_col.prop(settings, "layout_extra_offset")
        template_col.prop(settings, "clear_before_build")
        template_col.separator()
        template_col.operator("molmol.inspect_fragments", icon="VIEWZOOM")
        template_col.operator("molmol.build_from_file", icon="OUTLINER_COLLECTION")

        layout.operator("molmol.try_to_fix_bonds", icon="CONSTRAINT_BONE")
        layout.label(text="Select fragments and bond_independent roots to fix", icon="INFO")
        layout.label(text="Rigid fix only: rotation + translation, no scaling", icon="INFO")

        rdkit_box = layout.box()
        rdkit_box.label(text="RDKit Debug Build", icon="MESH_UVSPHERE")
        rdkit_col = rdkit_box.column(align=True)
        rdkit_col.prop(settings, "atom_radius")
        rdkit_col.prop(settings, "bond_radius")
        rdkit_col.prop(settings, "coordinate_scale")
        rdkit_col.prop(settings, "add_hydrogens")
        rdkit_col.prop(settings, "show_hydrogens")
        rdkit_col.prop(settings, "center_rdkit")
        rdkit_col.prop(settings, "show_labels")
        rdkit_col.separator()
        rdkit_col.operator("molmol.build_rdkit_balls_sticks", icon="MESH_UVSPHERE")

        tools_box = layout.box()
        tools_box.label(text="Cleanup", icon="TRASH")
        tools_col = tools_box.column(align=True)
        tools_col.operator("molmol.clear_built", icon="TRASH")
        
CLASSES = (
    MOLMOL_PG_Settings,
    MOLMOL_OT_build_from_file,
    MOLMOL_OT_build_rdkit_balls_sticks,
    MOLMOL_OT_inspect_fragments,
    MOLMOL_OT_clear_built,
    MOLMOL_OT_try_to_fix_bonds,
    MOLMOL_PT_panel,
)



def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)

    bpy.types.Scene.molmol_settings = PointerProperty(type=MOLMOL_PG_Settings)


def unregister():
    if hasattr(bpy.types.Scene, "molmol_settings"):
        del bpy.types.Scene.molmol_settings

    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)