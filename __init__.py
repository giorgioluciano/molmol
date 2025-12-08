bl_info = {
    "name": "Molymod Atom Builder (modular, PDB/CIF/ASE)",
    "author": "Giorgio,G.",
    "version": (2, 0, 0),
    "blender": (4, 0, 0),
    "location": "3D View > N-Panel > Molymod",
    "description": "Build Molymod-style molecules from PDB and CIF using external library collections, compact factor, caps, robust orientation, and auto-bond inference.",
    "category": "Import-Export",
}

import bpy

# Import local modules for registration
from .properties import MOLYMOD_PG_Settings
from .panel import MOLYMOD_PT_Panel
from .operators import MOLYMOD_OT_Build, MOLYMOD_OT_ValidateLibrary, MOLYMOD_OT_ClearAll


# Tuple of all custom classes
classes = (
    MOLYMOD_PG_Settings,
    MOLYMOD_OT_Build,
    MOLYMOD_OT_ValidateLibrary,
    MOLYMOD_OT_ClearAll,
    MOLYMOD_PT_Panel,
)

def register():
    """Register all Molymod classes and add property group to Scene."""
    for c in classes:
        bpy.utils.register_class(c)
    bpy.types.Scene.molymod_settings = bpy.props.PointerProperty(type=MOLYMOD_PG_Settings)

def unregister():
    """Unregister all classes and remove property group."""
    del bpy.types.Scene.molymod_settings
    for c in reversed(classes):
        bpy.utils.unregister_class(c)

# Enable running as __main__ for script reload/testing
if __name__ == "__main__":
    register()
