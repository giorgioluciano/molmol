bl_info = {
    "name": "MolMol",
    "author": "MolMol",
    "version": (0, 2, 0),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar > MolMol",
    "description": "Import MOL/SDF, detect fragments, build balls-and-sticks or template-based assemblies",
    "category": "Object",
}

from . import properties_and_ui

register = properties_and_ui.register
unregister = properties_and_ui.unregister