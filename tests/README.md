# Spectra Comparison Test Pack

This folder is the **quick reference + runnable test harness** for comparing computed spectra against bundled experimental overlays.

## Files

- `tests/spectra_comparison_cases.csv`
  - Canonical mapping for each case:
    - difficulty
    - case name
    - workflow (`nmr` or `cd`)
    - target spectrum product (`1H`, `13C`, `19F`, `31P`, `CD`)
    - `--nucleus` argument (when relevant)
    - SMILES
    - computed reference file
    - experimental overlay file

- `tests/run_spectra_comparison_matrix.sh`
  - One-command matrix check across all cases.
  - Validates experimental files with `easynmr-expcheck`.
  - Runs prediction for each row and verifies the expected product label appears in `spectra_manifest.csv`.

## Fast start

Run all matrix cases:

```bash
./tests/run_spectra_comparison_matrix.sh
```

Run with a custom output directory:

```bash
./tests/run_spectra_comparison_matrix.sh output_test_matrix_local
```

Run only matching case names (substring filter):

```bash
./tests/run_spectra_comparison_matrix.sh output_test_matrix_local 31p
./tests/run_spectra_comparison_matrix.sh output_test_matrix_local cd
```

## Manual GUI comparison workflow

For any row in `tests/spectra_comparison_cases.csv`:

1. Copy the row SMILES into GUI input.
2. Set `Workflow` to row `workflow`.
3. If workflow is `nmr`, select row `nucleus_arg` (or `auto` where listed).
4. Run prediction.
5. Use `Load Exp` and choose row `experimental_overlay_file`.
6. Select the matching computed product (`target_product`) in the spectrum selector.

This gives a reproducible computed-vs-experimental comparison path for CD and NMR families.
