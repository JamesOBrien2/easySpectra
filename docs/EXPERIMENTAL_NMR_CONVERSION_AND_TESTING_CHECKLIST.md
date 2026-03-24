# Experimental NMR Conversion and Testing Checklist

Use this checklist whenever you modify parsing, plotting, overlay rendering, workflow wiring, or output contracts.

## 1) Keep canonical external datasets in place

Canonical tier files are stored in:

- `examples/external_nmr_pack/easy/raw/20-1H.*`
- `examples/external_nmr_pack/medium/raw/20-13C.*`
- `examples/external_nmr_pack/hard/raw/20-hsqc.*`

Provenance and DOI references are in:

- `examples/external_nmr_pack/README.md`

## 2) Convert raw/vendor files to overlay-ready XY files

EasyNMR currently plots experimental overlays from **2-column text/CSV** only.

For each tier, export one overlay-ready file and place it in:

- `examples/external_nmr_pack/<tier>/converted/`

Accepted extensions:

- `.csv`, `.txt`, `.asc`, `.dat`

Recommended export settings (Mnova/TopSpin):

- X axis: chemical shift (ppm)
- Y axis: intensity
- Decimal point format preferred (`.`)
- No extra columns
- Keep enough points for a continuous trace (typically 500+)

## 3) Run external-pack regression checks

```bash
./scripts/regression_external_nmr_pack.sh
```

This verifies:

- `.mnova` direct import fails with the expected guidance message.
- extracted Bruker directories fail with the expected guidance message.
- every file in `converted/` parses successfully with `easynmr-expcheck`.

## 4) Run product-level smoke tests

```bash
./scripts/smoke_products_and_experimental.sh
```

This verifies:

- build succeeds,
- NMR (including `31P`) and CD workflows still run,
- baseline experimental parser fixtures still pass.

## 4b) Refresh and validate the extended bundled example pack

```bash
./scripts/generate_example_pack.py
./scripts/smoke_products_and_experimental.sh
```

This keeps bundled easy/medium/hard overlays for CD, `13C`, `19F`, and `31P` synchronized with current simulation behavior.

## 5) Manual GUI confirmation

For at least one converted file per tier:

1. Launch `./build/easynmr-gui`
2. Run or load a computed spectrum.
3. Click `Load Exp` and select the converted file.
4. Confirm overlay behavior:
   - same x-axis alignment,
   - experimental trace drawn below zero (negative y),
   - complementary color remains readable,
   - zoom/select interactions still behave correctly.

## 6) Record results

When done, log:

- date,
- branch/commit,
- converted files used,
- pass/fail per tier,
- any visual issues or parser warnings.
