# easySpectra

`easySpectra` is a local app for predicting and comparing spectra without needing cloud services.

It is built for ease of use:
- paste a structure,
- run predictions,
- load experimental data,
- compare visually,
- export a clean figure.

Current version: `0.0.1`

## What You Can Do Today

- Predict NMR spectra (`1H`, `13C`, `19F`, `31P`)
- Run a CD prediction path (`cd` workflow scaffold)
- Run everything at once (`all` workflow, default)
- Overlay experimental traces against computed traces
- Switch between multiple loaded experimental overlays in the GUI
- Export the plotted spectrum view (`.png` or `.ppm`)

## 5-Minute Setup

### 1) Install Python dependencies

```bash
pip install -r backend/requirements.txt
```

### 2) Build the project

```bash
cmake -S . -B build
cmake --build build -j
```

### 3) Launch the GUI

```bash
./build/easynmr-gui
```

Note on command names:
- The product name is `easySpectra`.
- Current binary names are still `easynmr`, `easynmr-gui`, and `easynmr-expcheck` for compatibility.

## Quick GUI Workflow

1. Paste SMILES or structure text.
2. Leave `Workflow` as `all` (default), or choose `nmr`/`cd`.
3. Click `Preview`.
4. Click `Queue`, then `Run Pending`.
5. Use `Load Exp` to load an experimental file.
6. Switch overlays with the `Exp:` selector (`none` or a loaded file).
7. Click `Export` to save exactly what you see in the spectrum panel.

## Quick CLI Workflow

Basic run:

```bash
./build/easynmr --input "CCO" --input-format smiles --name ethanol
```

Explicit workflow selection:

```bash
./build/easynmr --workflow all --input "CCO" --input-format smiles
./build/easynmr --workflow nmr --input "CCO" --input-format smiles
./build/easynmr --workflow cd  --input "CCO" --input-format smiles
```

Helpful options:
- `--workflow <all|nmr|cd>`
- `--nucleus <auto|1H|13C|19F|31P>`
- `--solvent <cdcl3|dmso|h2o>`
- `--frequency-mhz <number>`
- `--line-shape <lorentzian|gaussian|voigt>`
- `--fwhm-hz <number>`

Outputs are written to `output/<job-id>/`.

## Comparing Computed vs Experimental

### Supported experimental inputs

- Bruker-like exported text
- MNova-like exported text
- Generic 2-column text/CSV (`x, y`)

### Current limits

- Bruker raw acquisition folders are not parsed directly yet.
- MNova project files (`.mnova`, `.mnv`) are not parsed directly yet.
- Export vendor data to 2-column text/CSV first, then load it.

### Validate a file before loading

```bash
./build/easynmr-expcheck <experimental-spectrum-file>
```

In the GUI, experimental overlays are drawn on the same x-axis with negative y-values for easy visual comparison against computed traces.

## Easy/Medium/Hard Test Packs

### Main case maps (includes SMILES)

- `examples/benchmark_cases.csv`
- `tests/spectra_comparison_cases.csv`

### Comparison table (SMILES + computed vs experimental)

Use `tests/spectra_comparison_cases.csv` as the main table when you want direct one-line mapping for each test case.

Key columns:
- `smiles`
- `workflow`
- `target_product`
- `computed_reference_csv`
- `experimental_overlay_file`

### Bundled examples

- `examples/computed/` for computed sample outputs
- `examples/experimental/` for experimental and converted overlays
- `examples/external_nmr_pack/` for curated vendor-style source files and provenance
- `examples/figures/` for optional visual comparison figures

### Visual examples (figures)

Store figures in `examples/figures/` and link/embed them here.

Featured figures:

![easySpectra GUI overview](examples/figures/easy_1H_nmr_comparison_highlighted_full_GUI.png)

*Figure 1. easySpectra GUI overview.*

![easy 1H computed vs experimental comparison](examples/figures/easy_1H_nmr_comparison.png)

*Figure 2. Easy 1H computed vs experimental overlay comparison export.*

Suggested gallery table:

| Case | Workflow/Product | SMILES | Figure |
|---|---|---|---|
| easy_1h_gui | NMR `1H` | `CCO` | ![easy 1H GUI](examples/figures/easy_1H_nmr_comparison_highlighted_full_GUI.png) |
| easy_1h_export | NMR `1H` | `CCO` | ![easy 1H export](examples/figures/easy_1H_nmr_comparison.png) |

### Run the smoke checks

```bash
./scripts/smoke_products_and_experimental.sh
```

### Run the full comparison matrix

```bash
./tests/run_spectra_comparison_matrix.sh
```

### Regenerate the bundled example pack

```bash
./scripts/generate_example_pack.py
```

## Where Important Output Files Live

Each run typically includes:
- `spectra_manifest.csv`
- `audit.json`
- `structure.svg`
- `structure_atoms.csv`
- `structure_bonds.csv`
- spectrum-specific files such as `spectrum_1h.csv`, `peaks_1h.csv`, `assignments_1h.csv` (and equivalent for `13C`, `19F`, `31P`, `cd`)

## Optional Backend Controls

- `EASYSPECTRA_XTB=/path/to/xtb` sets xTB path
- `EASYSPECTRA_XTB_TIMEOUT=25` sets per-conformer timeout (seconds)
- `EASYSPECTRA_XTB=__none__` forces MMFF fallback mode

## Extra Project Docs

- `docs/ARCHITECTURE.md`
- `docs/ROADMAP.md`
- `docs/SPEC.md`
- `docs/EXPERIMENTAL_NMR_CONVERSION_AND_TESTING_CHECKLIST.md`
- `docs/CHANGELOG.md`

## Tested NMR References

- Easy (`20-1H`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/11523](https://doi.org/10.14469/hpc/11523)
- Medium (`20-13C`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/11524](https://doi.org/10.14469/hpc/11524)
- Hard (`20-hsqc`, Bruker + Mnova + JCAMP): DOI [10.14469/hpc/13944](https://doi.org/10.14469/hpc/13944)
