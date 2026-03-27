# TODO

This is a practical working backlog for `easySpectra`.

## Now (Polish and Stability)

- [ ] Verify export-cropping behavior across different zoom ranges and window sizes, and make sure labels are never clipped.
- [ ] Improve axis/tick font sizing for very small exported images.
- [x] Add keyboard shortcuts for common actions: queue, run pending, load experimental, export.
- [x] Add an in-app hint when experimental file format is invalid, including exact expected column format.

## Near-Term (Comparison Quality)

- [ ] Add optional auto-alignment of computed vs experimental traces (shift-only first).
- [ ] Add optional intensity scaling/normalization for overlays.
- [ ] Add basic baseline correction for imported experimental traces.
- [ ] Add an overlay offset slider and style controls (line width, opacity).
- [ ] Add a quick metrics panel (for example MAE/correlation in selected ppm window).

## Near-Term (Workflow/Product Improvements)

- [ ] Move CD from scaffold mode to a fuller prediction path with documented assumptions.
- [ ] Add product-specific presets (for example publication-style NMR and CD export presets).
- [ ] Add optional batch mode from CSV (SMILES + workflow + solvent) in GUI and CLI.
- [ ] Support custom solvent/reference configuration from a user-editable file.

## Property Features (Prioritized by xTB Fit)

Legend:
- Horizon: `Near` = practical with current stack, `Mid` = possible with extra rules/work, `Far` = likely needs additional engines or data.
- Difficulty: `Low`, `Medium`, `High`.

- [ ] IR prediction v1 (frequencies + intensities + broadened plot) from xTB outputs. Horizon: `Near`, Difficulty: `Medium`.
- [ ] IR quality controls (scaling factor presets, line broadening presets, and simple peak labels). Horizon: `Near`, Difficulty: `Low`.
- [ ] Exact mass and isotope pattern calculator from molecular formula (quick MS support without fragmentation). Horizon: `Near`, Difficulty: `Low`.
- [ ] Common adduct support for MS mode (`[M+H]+`, `[M+Na]+`, `[M-H]-`, etc.) with m/z table output. Horizon: `Near`, Difficulty: `Low`.
- [ ] Property panel for per-conformer and ensemble values (dipole, HOMO/LUMO gap, relative energies). Horizon: `Near`, Difficulty: `Medium`.
- [ ] Basic thermochemistry summary from current calculations (relative free energy-style comparisons where available). Horizon: `Mid`, Difficulty: `Medium`.
- [ ] Mass-spec fragmentation prototype (rule-based fragments for simple molecules). Horizon: `Mid`, Difficulty: `High`.
- [ ] Mass-spec fragmentation scoring and ranking against experimental peaks. Horizon: `Far`, Difficulty: `High`.
- [ ] Raman prediction module. Horizon: `Far`, Difficulty: `High`.

## Data and File Support

- [ ] Add a conversion helper for common Bruker/MNova exports into 2-column overlay files.
- [ ] Add direct import support for additional text-like vendor formats where feasible.
- [ ] Expand bundled example coverage with more chemically diverse molecules per nucleus/workflow.
- [ ] Add figure gallery curation script to keep README figures synchronized with current behavior.

## Testing and Reliability

- [ ] Add automatic checks to confirm exported images are cropped correctly and labels are not cut off.
- [ ] Add file-import checks using broken or unusual experimental files, so failures are clear and safe.
- [ ] Add checks to confirm every spectrum listed in `spectra_manifest.csv` can be opened and displayed.
- [ ] Add simple automated checks on each change: build, run a quick example, and run the full example table.
- [ ] Add a few difficult example files (very broad peaks, sparse spectra, empty columns) and keep them in the test set.

## Packaging and Developer Experience

- [ ] Add one-command setup script for local environment bootstrap.
- [ ] Add release checklist for tag + changelog + artifact generation.
- [ ] Add platform packaging notes (macOS app bundle first, Linux next).
- [ ] Standardize user-facing naming to `easySpectra` while preserving CLI compatibility.
