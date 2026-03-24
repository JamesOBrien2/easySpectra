# EasyNMR

EasyNMR is a local-first spectroscopy workbench with both a desktop GUI and a CLI.
It currently supports NMR prediction and a scaffold CD workflow, with a modular workflow design so additional spectra products can be added without reshaping the whole codebase.

Pre-release version: `0.0.1-alpha.0`

## Why this project structure

EasyNMR is organized around **workflow products**, not just one fixed NMR path.

- `workflow.kind` (for example `all`, `nmr`, `cd`) selects the computational route.
- Outputs are emitted through a shared artifact contract (`spectra_manifest.csv`, `audit.json`, per-product CSVs).
- The GUI loads spectra products from the manifest, so new spectra types can plug in cleanly.

This keeps NMR as a first-class product while leaving a clear path for CD and future computational spectroscopy workflows.

## Current capabilities

- C++17 core with JSON bridge to a local Python backend
- CLI runner: `easynmr`
- GUI app (when FLTK is available): `easynmr-gui`
- NMR workflow (`1H`, `13C`, `19F`, `31P`) with conformer ensemble workflow
- CD scaffold workflow (`--workflow cd`) sharing the same conformer pipeline
- Combined workflow mode (`--workflow all`, default) to generate all currently available products in one run
- Product manifest loading in GUI (`spectra_manifest.csv`)
- Experimental overlay import in GUI (Bruker/MNova text exports + generic 2-column text/CSV)
- Overlay rendering on the same x-axis with **experimental trace on negative y-axis**
- Experimental overlay selector in GUI (`Exp: none` + loaded overlays) for fast switching while comparing products
- Interactive linking between peak groups and assigned atoms/hydrogens

## Requirements

- CMake `>= 3.20`
- C++17 compiler
- Python 3
- Python backend deps:

```bash
pip install -r backend/requirements.txt
```

Optional:

- FLTK (for GUI build)
- xTB binary for xTB optimization path

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

CLI-only build (skip GUI):

```bash
cmake -S . -B build -DEASYNMR_BUILD_GUI=OFF
cmake --build build -j
```

## Quick start: CLI

Basic run:

```bash
./build/easynmr --input "CCO" --input-format smiles --name ethanol
```

By default, EasyNMR runs all available workflows (`--workflow all`).
You can still select a specific workflow explicitly:

```bash
./build/easynmr --workflow all --input "CCO" --input-format smiles
./build/easynmr --workflow nmr --input "CCO" --input-format smiles
./build/easynmr --workflow cd  --input "CCO" --input-format smiles
```

Common useful flags:

- `--workflow <all|nmr|cd>`
- `--nucleus <auto|1H|13C|19F|31P>` (NMR workflow)
- `--solvent <cdcl3|dmso|h2o>`
- `--frequency-mhz <number>`
- `--line-shape <lorentzian|gaussian|voigt>`
- `--fwhm-hz <number>`

Outputs are written to `output/<job-id>/` by default.

## Quick start: GUI

```bash
./build/easynmr-gui
```

Typical flow:

1. Enter structure input (SMILES or structure text).
2. Choose `Workflow` (`all`, `nmr`, or `cd`) and `Format`.
3. Click `Preview`.
4. Click `Queue`, then `Run Pending` (or `Run Selected`).
5. Click a peak or peak-group row to inspect assignments.
6. Use `Load Exp` to import experimental spectra overlays; use `Clear Exp` to remove them.
7. Use `Export` in the spectrum panel to save the current plotted view (`.png`/`.ppm`).

## Experimental spectra comparison

Supported import style today:

- Bruker-like text export
- MNova-like text export
- Generic two-column text/CSV (x, y)

Important current limits:

- Bruker raw acquisition folders are not parsed directly yet.
- MNova project files (`.mnova`, `.mnv`) are not parsed directly.
- For both cases above, export to a two-column text/CSV first, then import.

Validation utility:

```bash
./build/easynmr-expcheck <experimental-spectrum-file>
```

## Artifact layout

A completed run writes product outputs into the job directory, including:

- `spectra_manifest.csv` (product labels + file paths)
- `audit.json`
- `structure.svg`
- `structure_atoms.csv`
- `structure_bonds.csv`
- product files such as:
  - `spectrum_1h.csv`, `peaks_1h.csv`, `assignments_1h.json/csv`
  - `spectrum_13c.csv`, `peaks_13c.csv`, ...
  - `spectrum_19f.csv`, `peaks_19f.csv`, ...
  - `spectrum_31p.csv`, `peaks_31p.csv`, ...
  - `spectrum_cd.csv`, `peaks_cd.csv`, `assignments_cd.json/csv` (scaffold)

## Example datasets and smoke checks

Benchmark files are included for quick validation:

- `examples/benchmark_cases.csv` (easy / medium / hard registry)
- `examples/computed/` (computed sample spectra)
- `examples/experimental/`:
  - `easy_mnova_export.txt`
  - `medium_bruker_export.txt`
  - `hard_generic_export.csv`
  - generated easy/medium/hard overlays for CD, `13C`, `19F`, and `31P`

Generate or refresh the extended example pack:

```bash
./scripts/generate_example_pack.py
```

Run the bundled smoke script:

```bash
./scripts/smoke_products_and_experimental.sh
```

This script:

1. builds the project,
2. runs NMR, CD, and `31P` smoke predictions,
3. validates baseline and extended experimental overlay files with `easynmr-expcheck`.

External curated regression pack:

- `examples/external_nmr_pack/README.md` documents canonical easy/medium/hard vendor datasets (Bruker/MNova/JCAMP).
- `./scripts/regression_external_nmr_pack.sh` checks raw-file failure behavior and validates converted text/CSV files (when present).
- Conversion and recurring test procedure: `docs/EXPERIMENTAL_NMR_CONVERSION_AND_TESTING_CHECKLIST.md`.

Extended bundled comparison pack:

- Includes easy/medium/hard computed + experimental overlays for:
  - CD
  - `13C` NMR
  - `19F` NMR
  - `31P` NMR
- Case registry (including SMILES): `examples/benchmark_cases.csv`

Tests folder mapping + runnable matrix:

- `tests/spectra_comparison_cases.csv` gives direct mapping of
  - SMILES
  - workflow/nucleus
  - computed reference file
  - experimental overlay file
- `./tests/run_spectra_comparison_matrix.sh` runs the full comparison matrix checks.

## Backend controls

Environment variables:

- `EASYNMR_XTB=/path/to/xtb` to set xTB binary path
- `EASYNMR_XTB_TIMEOUT=25` to cap per-conformer xTB wall time (seconds)
- `EASYNMR_XTB=__none__` to force MMFF fallback mode

## Repository docs

- `docs/SPEC.md` - locked v1 scope
- `docs/ARCHITECTURE.md` - implementation architecture and data flow
- `docs/ROADMAP.md` - near-term development plan
- `docs/EXPERIMENTAL_NMR_CONVERSION_AND_TESTING_CHECKLIST.md` - recurring external dataset conversion and regression checklist

## Tested NMR References

- Easy (`20-1H`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/11523](https://doi.org/10.14469/hpc/11523)
- Medium (`20-13C`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/11524](https://doi.org/10.14469/hpc/11524)
- Hard (`20-hsqc`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/13944](https://doi.org/10.14469/hpc/13944)
