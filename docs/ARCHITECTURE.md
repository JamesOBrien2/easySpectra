# Architecture

## Components

- `src/core/`: C++ job model, workflow-aware pipeline bridge, sampled spectrum loader, and spectral-product manifest model.
- `src/cli/`: CLI front-end (`easynmr`).
- `src/gui/`: FLTK desktop prototype (`easynmr-gui`).
- `backend/`: local Python backend worker (`easynmr_backend.py`).

## Data Flow

1. Frontend builds `JobConfig` with a `workflow_kind` (currently `nmr`).
2. Core pipeline writes `request.json` to job output directory, including `workflow.kind`.
3. Core invokes local Python backend with request/response paths.
4. Backend writes:
   - `spectrum.csv`
   - `peaks.csv`
   - `assignments.json`
   - `assignments.csv`
   - `spectra_manifest.csv` (spectrum label + file paths per product)
   - `structure.svg`
   - `audit.json`
   - `response.json`
5. Frontend reads outputs and updates UI.

## Product Modularity

- Frontend spectrum loading uses product labels from `spectra_manifest.csv` instead of a hard-coded single-spectrum path.
- NMR (`1H`, `13C`, `19F`) is currently the implemented product family.
- The same manifest contract can be extended for additional spectra (for example CD) without changing the core transport layer.

## Why hybrid C++ + Python

- C++/FLTK provides native desktop UX and close alignment with existing xyz tooling.
- Python backend accelerates chemistry workflow iteration and integration with existing QC/ML tools.
- JSON artifact contract keeps components auditable and swappable.
