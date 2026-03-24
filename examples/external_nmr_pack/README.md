# External NMR Pack (Curated)

This folder contains a curated external dataset pack for recurring regression checks against real Bruker/MNova-style sources.

Canonical tiers:

- `easy`: `20-1H` (1D proton)
- `medium`: `20-13C` (1D carbon)
- `hard`: `20-hsqc` (2D HSQC, more complex archive layout)

Each tier includes raw files under `raw/`:

- Bruker archive (`*.zip`)
- MestreNova project (`*.mnova`)
- JCAMP (`*.jdx`)
- Mnova signature (`*.mnpub`)

Converted files for in-app overlay tests should be exported as two-column text/CSV and stored in `converted/`.

## Provenance

All files in this pack are from Imperial College data records by Henry Rzepa (CC0):

- Easy (`20-1H`): DOI [10.14469/hpc/11523](https://doi.org/10.14469/hpc/11523)
- Medium (`20-13C`): DOI [10.14469/hpc/11524](https://doi.org/10.14469/hpc/11524)
- Hard (`20-hsqc`): DOI [10.14469/hpc/13944](https://doi.org/10.14469/hpc/13944)

Related additional records retained in this folder for optional use:

- DOI [10.14469/hpc/11520](https://doi.org/10.14469/hpc/11520)
- DOI [10.14469/hpc/11525](https://doi.org/10.14469/hpc/11525)

## Regression command

Run:

```bash
./scripts/regression_external_nmr_pack.sh
```

For conversion/export steps, see:

- `docs/EXPERIMENTAL_NMR_CONVERSION_AND_TESTING_CHECKLIST.md`
