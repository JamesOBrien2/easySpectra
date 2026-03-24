# EasyNMR

EasyNMR is a local-first desktop + CLI tool for fast, explainable NMR prediction from
SMILES or structure files.

Current scaffold includes:
- C++17 core with JSON bridge to local Python backend
- `easynmr` CLI runner
- `easynmr-gui` FLTK prototype (when FLTK is available)
- RDKit + xTB fast-path backend with conformer generation, Boltzmann filtering, and first-order 1H simulation
- Export-ready artifact layout (`spectrum.csv`, `peaks.csv`, `assignments.json`, `audit.json`)
- 2D structure preview (`structure.svg`) and interactive peak-group to atom assignment linking in GUI

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

If FLTK is unavailable, CLI will still build:

```bash
cmake -S . -B build -DEASYNMR_BUILD_GUI=OFF
cmake --build build -j
```

## Quick Start (CLI)

```bash
./build/easynmr --input "CCO" --input-format smiles --name ethanol
```

Outputs are written to `output/<job-id>/` by default.

Backend controls:
- `EASYNMR_XTB=/path/to/xtb` to override xTB binary
- `EASYNMR_XTB_TIMEOUT=25` to limit per-conformer xTB wall-time seconds
- set `EASYNMR_XTB=__none__` to force MMFF fallback mode

## Quick Start (GUI)

```bash
./build/easynmr-gui
```

In the GUI prototype:
- enter a SMILES/string input,
- choose input format,
- click `Queue Job` to enqueue,
- click `Start Queue` to process pending jobs,
- use `Retry Selected` to requeue failed jobs,
- click `Cancel Active` to stop after current running step,
- select a peak group (browser or direct click on spectrum) to see linked hydrogens.

## Docs

- `docs/SPEC.md` - locked v1 scope
- `docs/ARCHITECTURE.md` - implementation architecture and data flow
- `docs/ROADMAP.md` - 1-2 week prototype + 1-2 month usable target plan
