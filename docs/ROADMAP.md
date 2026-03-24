# Roadmap

This roadmap prioritizes what can be delivered well with the current stack first (`xTB` + existing C++/Python workflow), then stages features that need heavier heuristics or additional engines.

## Current Baseline

- Local GUI + CLI workflow with queueing and artifact exports.
- Multi-product pipeline with `all`, `nmr`, and `cd`.
- NMR support for `1H`, `13C`, `19F`, `31P`.
- Experimental overlay loading and side-by-side computed vs experimental comparison.

## Near-Term (Best Fit for Current xTB Stack)

- Complete CD from scaffold to a stronger production path.
- IR prediction v1 using xTB-derived frequencies/intensities with practical plotting controls.
- Exact mass + isotope pattern support (formula-based, quick MS utility).
- Adduct table support for common positive/negative ionization modes.
- Property summaries from current runs (dipole, HOMO/LUMO gap, conformer energy spread).
- Continued GUI polish for export quality and comparison UX.

Why near-term:
- These features are either directly available from current calculations, or can be added with lightweight post-processing.

## Mid-Term (xTB-Assisted + Extra Heuristics)

- Basic mass-spec fragmentation prototype (rule-based fragments for common motifs).
- Better computed/experimental alignment tools (auto-shift, normalization options, simple fit metrics).
- Broader import tooling for vendor exports and conversion helpers.
- Batch workflow enhancements for larger curated comparison suites.

Why mid-term:
- Feasible, but requires substantial rules/heuristics and validation to avoid misleading outputs.

## Far-Term (Likely Requires Additional Engines or Data)

- Higher-quality MS/MS fragmentation prediction and confidence scoring.
- Raman module with usable quality targets.
- Optional high-accuracy engine integrations for selected products (for example ORCA-backed workflows).
- Deeper quantitative matching layers for experimental/computed benchmarking.

Why far-term:
- These items typically need methods beyond current xTB-centered capabilities, plus broader benchmarking datasets.
