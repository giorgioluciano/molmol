# Atom Builder for Blender

A Blender add-on for **scientific molecular 3D modeling** using the iconic 
teaching model  atomic representation style. Build publication-ready molecular structures 
with **automatic bond detection** via the Atomic Simulation Environment (ASE), 
customizable atom geometries, and advanced rendering options.

## Features

- **real-style atoms**: Import or build molecules with the classic colored-sphere-and-stick representation used in scientific publications.
- **Automatic bond detection**: Uses ASE to infer bonds from atomic coordinates (supports PDB, CIF formats).
- **Flexible atom collections**: Reference external Blender `.blend` libraries with sp³, sp², sp, bent, and other hybridization geometries.
- **Bond customization**: Control bond radius, gap distance, and hydrogen handling.
- **Publication-ready**: Generate high-quality renders suitable for academic papers, theses, and presentations.
- **Batch processing**: Import and process multiple structures in a single workflow.

## Use Cases

- Academic research and scientific visualization
- Materials science and computational chemistry
- Educational content and molecular modeling
- 3D printing of molecular structures
- Publication figures for chemistry and materials papers

## Requirements

- **Blender 5.0+** (or 4.2+, depending on your actual minimum)
- **Python 3.11+** (bundled with Blender)
- **Dependencies** (installed via pip):
  - `ase` (Atomic Simulation Environment) — for bond detection and structure parsing
  - `numpy` — required by ASE
  - `scipy` (optional) — for advanced bond assignment algorithms

## Installation

### Step 1: Locate Blender's Python Interpreter

Blender includes its own Python environment. You need to install packages into it.

#### On Linux:

```

# Example path (adjust version and installation path as needed)

export BLENDER_PYTHON="/snap/blender/current/5.0/python/bin/python3.11"

# or if installed via system package manager:

export BLENDER_PYTHON="\$HOME/blender-5.0-linux-x64/5.0/python/bin/python3.11"

```

#### On macOS:

```

export BLENDER_PYTHON="/Applications/Blender.app/Contents/Resources/5.0/python/bin/python3.11"

```

#### On Windows:

```

set BLENDER_PYTHON="C:\Program Files\Blender Foundation\Blender 5.0\5.0\python\bin\python.exe"

```

### Step 2: Install Required Packages

Open a terminal/command prompt and run:

```

\$BLENDER_PYTHON -m ensurepip
\$BLENDER_PYTHON -m pip install --upgrade pip
\$BLENDER_PYTHON -m pip install ase numpy scipy

```

**Verification** — To confirm the installation worked:

```

\$BLENDER_PYTHON -c "import ase; print(f'ASE version: {ase.__version__}')"

```

### Step 3: Install the Addon

1. Download the addon `.zip` file from the [Releases](../../releases) page.
2. Open **Blender** → **Edit** → **Preferences** → **Add-ons**
3. Click **Install** and select the downloaded `.zip` file
4. Search for "Atomic Structure Builder" and enable it (check the box)
5. The addon will appear in the N-panel (press `N`) under the "Atomic Structure" tab

### Step 4: Configure the Addon

1. In Blender, open the N-panel ("Atomic Structure" tab)
2. Click **"Locate Library File"** and select the example `.blend` library 
   (included in the addon folder)
3. Test with **"Validate Library"** — all required collections should be listed

---

## Troubleshooting

### ImportError: No module named 'ase'

- Confirm you installed ASE in the correct Python interpreter
- Verify the path with `$BLENDER_PYTHON -c "import ase"`
- Try reinstalling: `$BLENDER_PYTHON -m pip install --force-reinstall ase`

### Library file not found

- Make sure the `.blend` library file path is correctly set in the addon panel
- The library must contain collections named: `Atom_sp3`, `Atom_sp2`, `Atom_sp`, `Atom_bent`

### Addon doesn't appear in Preferences

- Restart Blender after installation
- Check that the addon folder structure is correct (should have `__init__.py`)
  
### see also https://www.youtube.com/@jojo75cg
---
```



