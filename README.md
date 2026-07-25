# MolMol

MolMol is a Blender-oriented molecular assembly pipeline that converts RDKit molecular structures into a Blender scene using fragment templates, planning logic, and marker-based topology resolution.

## What it does

MolMol takes an RDKit molecular structure and builds a Blender scene through:

- RDKit-guided initial fragment layout
- planner / assembly pipeline
- fragment template library
- marker/contact-based connection resolution

Important: marker-based resolution is used to determine valid fragment connections and topology, not to reposition fragments away from the initial RDKit-guided layout.

## Stable in v1.0

- RDKit -> planner -> Blender pipeline
- verified fragment template naming
- initial fragment placement guided by RDKit geometry
- dedicated handling for independent single bonds
- marker/contact resolution for topology only

## Dependencies

### Required

- Blender
- Python
- RDKit

### Install RDKit with pip

```bash
python -m pip install rdkit
```

### Alternative: install RDKit with conda

```bash
conda install -c conda-forge rdkit
```

## Installation

1. Clone this repository.
2. Install RDKit in the Python environment used by the pipeline.
3. Copy or configure the fragment library file separately.
4. Open Blender and install/enable the required addon components if applicable.
5. Run the MolMol pipeline on an RDKit molecule input.

## Fragment library

The fragment library file is distributed separately.
Place it in the expected project path before running the pipeline.

## Architecture

- RDKit defines molecular structure and supports initial geometry decisions.
- The planner determines fragment assembly order and placement logic.
- Blender receives the initial fragment transforms.
- Marker/contact resolution validates and resolves allowed connections between fragments.

## Notes

MolMol is focused on correct initial placement plus topology resolution.
It should not use marker solving as a secondary system to re-layout already correct fragment placements.


