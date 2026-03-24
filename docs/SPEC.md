# EasyNMR v1 Spec (Locked)

## MVP

- Input: SMILES, MOL/SDF, XYZ, and embedded 2D structure editor.
- Output: 1H NMR simulation at 400 MHz in ppm and Hz.
- Interaction: atom-to-peak and peak-to-atom linking.
- Batch: up to 50 molecules.
- Runtime target: minutes for small organics (around 40 heavy atoms).
- Warnings preferred over hard failure where possible.
- Local/offline only.
- License: MIT.

## Compute Direction (v1)

- Fast path first: xTB / AIMNet2 + conformer workflow.
- ORCA high-accuracy mode deferred beyond v1.
- Solvent support with ALPB-like workflow assumptions (DMSO, CDCl3, H2O priorities).
- Conformer filtering controls: Boltzmann threshold and energy window.
- Protonation/tautomer handling: auto mode with pH override + manual lock mode.

## Visualization and Exports

- Desktop GUI priority (macOS first, then Linux).
- Integrals included in v1.
- Default line shape as common NMR default; optional additional shapes if feasible.
- Exports: SVG, CSV (sampled spectrum + peak table), PDF.
- Provenance JSON emitted for each run.

## Quality Target

- Chemically sensible assignments.
- Visually plausible 1H spectral pattern.
- Reproducible settings + audit trail.
