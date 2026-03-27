# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**easySpectra** is a local desktop application for predicting and comparing NMR spectra (1H, 13C, 19F, 31P) and circular dichroism (CD) without cloud services. It uses a hybrid C++ GUI frontend + Python computational backend architecture.

## Build Commands

```bash
# Install Python backend dependencies
pip install -r backend/requirements.txt

# Configure and build
cmake -S . -B build
cmake --build build -j

# Run the GUI
./build/easynmr-gui

# Run the CLI
./build/easynmr --input "CCO" --input-format smiles --name ethanol --workflow all

# Validate an experimental NMR file
./build/easynmr-expcheck <path>
```

**Requirements**: CMake 3.20+, C++17, Python3, RDKit, NumPy. FLTK is optional (needed for GUI target).

**Build outputs**: `easynmr` (CLI), `easynmr-gui` (desktop app, if FLTK found), `easynmr-expcheck` (experimental file validator), `easynmr_core` (shared library).

## Testing

```bash
# Quick end-to-end smoke test
./scripts/smoke_products_and_experimental.sh

# Full test matrix (20+ molecules, multiple nuclei/workflows)
./tests/run_spectra_comparison_matrix.sh

# Regenerate bundled examples from current backend code
python scripts/generate_example_pack.py
```

Test cases live in `tests/spectra_comparison_cases.csv` (easy/medium/hard difficulty levels, covering 1H, 13C, 19F, 31P, CD).

## Architecture

### Component Separation

The C++ frontend and Python backend communicate via JSON files on disk. The C++ pipeline:
1. Writes `request.json` to a timestamped output directory
2. Spawns Python as a subprocess
3. Polls `progress.json` for real-time stage updates
4. Reads `response.json` when complete

This JSON contract keeps components independently auditable and swappable.

### Data Flow

```
User input (SMILES/MOL/SDF/XYZ)
  → C++ JobConfig
    → request.json
      → Python backend (RDKit, xTB/MMFF conformers, empirical shifts)
        → spectrum.csv, peaks.csv, structure.svg, structure_atoms.csv, structure_bonds.csv, spectra_manifest.csv, audit.json, response.json
          → C++ GUI renders interactive spectrum + structure
```

### Multi-Product Output

A single job can produce multiple spectra (1H, 13C, 19F, 31P, CD). Products are indexed via `spectra_manifest.csv` — never hardcoded paths. The GUI loads whichever products the backend produced.

### Key C++ Classes

- **`AppWindow`** (`src/gui/app_window.h/.cpp`): Central window managing input panel, job queue (threaded), workflow controls, and coordinates `SpectrumWidget` ↔ `StructureWidget` cross-highlighting.
- **`SpectrumWidget`** (`src/gui/spectrum_widget.h/.cpp`): Interactive plot with zoom/pan, peak picking, experimental overlay (inverted on negative y-axis), and comparison mode.
- **`StructureWidget`** (`src/gui/structure_widget.h/.cpp`): 2D molecular structure with atom click callbacks and stereochemistry rendering.
- **`Pipeline`** (`src/core/pipeline.h/.cpp`): Job orchestration — spawns backend, polls progress, parses outputs.
- **`JobConfig` / `JobOutputs`** (`src/core/job.h`): The canonical data contracts between pipeline and GUI.

### Python Backend (`backend/easynmr_backend.py`)

Key responsibilities:
- Molecule I/O: parses SMILES/MOL/SDF/XYZ via RDKit
- 3D conformer generation: ETKDG embedding → xTB GFN2 or MMFF94 optimization → Boltzmann weighting
- Empirical NMR shift prediction: group-based heuristics for 1H/13C/19F/31P
- J-coupling estimation from dihedral angles (Karplus-like)
- Spectrum simulation: splitting patterns + Lorentzian/Gaussian/Voigt convolution
- CD prediction: heuristic bands from chirality and heteroatom count
- SVG + coordinate CSV export for structure visualization

**xTB integration**: Uses `EASYSPECTRA_XTB` env var to locate the binary; falls back to MMFF94 if unavailable. Set `EASYSPECTRA_XTB=__none__` to force MMFF-only. Timeout controlled by `EASYSPECTRA_XTB_TIMEOUT` (default 25s).

## Solvent and Nucleus Options

- Solvents: `cdcl3`, `dmso`, `h2o`
- Nuclei: `auto`, `1h`, `13c`, `19f`, `31p`
- Workflows: `all`, `nmr`, `cd`

## Experimental File Formats

The backend auto-detects Bruker, MNova, and generic 2-column (ppm, intensity) CSV/text formats. See `docs/EXPERIMENTAL_NMR_CONVERSION_AND_TESTING_CHECKLIST.md` for format details.

## Docs

- `docs/ARCHITECTURE.md`: System design details
- `docs/SPEC.md`: Locked v1 specification (do not break these behaviors)
- `docs/TODO.md`: Active working backlog
- `docs/ROADMAP.md`: Feature planning
