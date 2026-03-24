# Architecture

## Components

- `src/core/`: C++ job model, workflow-aware pipeline bridge, sampled spectrum loader, and spectral-product manifest model.
- `src/cli/`: easySpectra CLI front-end (`easynmr`).
- `src/gui/`: easySpectra FLTK desktop prototype (`easynmr-gui`).
- `backend/`: easySpectra local Python backend worker (`easynmr_backend.py`).

## Data Flow

1. Frontend builds `JobConfig` with a `workflow_kind` (`all`, `nmr`, `cd`; default is `all`).
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
- NMR (`1H`, `13C`, `19F`, `31P`) and a first CD scaffold are currently implemented products.
- `workflow_kind=all` runs all currently available products in one pipeline execution.
- The same manifest contract can be extended for additional spectra without changing the core transport layer.

## Experimental Overlay

- GUI imports experimental spectra from 2-column text/CSV exports.
- Current parser recognizes Bruker-like and MNova-like text exports, plus generic CSV/text.
- Experimental traces are rendered on the same x-axis with negative y-axis orientation for direct visual comparison against computed traces.
- Multiple experimental overlays can be loaded and switched (`Exp: none` + loaded list) without re-importing files.

## Why hybrid C++ + Python

- C++/FLTK provides native desktop UX and close alignment with existing xyz tooling.
- Python backend accelerates chemistry workflow iteration and integration with existing QC/ML tools.
- JSON artifact contract keeps components auditable and swappable.
