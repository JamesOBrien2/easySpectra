# Architecture

## Components

- `src/core/`: C++ job model, pipeline bridge, spectrum loader.
- `src/cli/`: CLI front-end (`easynmr`).
- `src/gui/`: FLTK desktop prototype (`easynmr-gui`).
- `backend/`: local Python backend worker (`easynmr_backend.py`).

## Data Flow

1. Frontend builds `JobConfig`.
2. Core pipeline writes `request.json` to job output directory.
3. Core invokes local Python backend with request/response paths.
4. Backend writes:
   - `spectrum.csv`
   - `peaks.csv`
   - `assignments.json`
   - `assignments.csv`
   - `structure.svg`
   - `audit.json`
   - `response.json`
5. Frontend reads outputs and updates UI.

## Why hybrid C++ + Python

- C++/FLTK provides native desktop UX and close alignment with existing xyz tooling.
- Python backend accelerates chemistry workflow iteration and integration with existing QC/ML tools.
- JSON artifact contract keeps components auditable and swappable.
