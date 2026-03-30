#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_DIR="${1:-output_smoke_products}"
mkdir -p "$OUT_DIR"

echo "[1/5] Build"
cmake --build build -j >/dev/null

echo "[2/5] NMR workflow smoke run"
EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --input "CCO" \
  --input-format smiles \
  --workflow nmr \
  --name smoke_nmr \
  --output-dir "$OUT_DIR" >/dev/null

echo "[3/5] CD workflow smoke run"
EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --input "CCO" \
  --input-format smiles \
  --workflow cd \
  --name smoke_cd \
  --output-dir "$OUT_DIR" >/dev/null

echo "[4/5] 31P workflow smoke run"
EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --input "CCOP(=O)(OCC)OCC" \
  --input-format smiles \
  --workflow nmr \
  --nucleus 31P \
  --name smoke_31p \
  --output-dir "$OUT_DIR" >/dev/null

echo "[5/5] Experimental parser checks"
./build/easynmr-expcheck examples/experimental/easy_mnova_export.txt
./build/easynmr-expcheck examples/experimental/medium_bruker_export.txt
./build/easynmr-expcheck examples/experimental/hard_generic_export.csv
./build/easynmr-expcheck examples/experimental/easy_cd_overlay.csv
./build/easynmr-expcheck examples/experimental/easy_13c_overlay.csv
./build/easynmr-expcheck examples/experimental/easy_19f_overlay.csv
./build/easynmr-expcheck examples/experimental/easy_31p_overlay.csv

echo "[6/8] Level-of-theory: low (MMFF-only, 5 conformers)"
EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --input "CCCCCCCC" \
  --input-format smiles \
  --workflow nmr \
  --level low \
  --name smoke_lot_low \
  --output-dir "$OUT_DIR" >/dev/null

echo "[7/8] Level-of-theory: high"
EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --input "CCCCCCCC" \
  --input-format smiles \
  --workflow nmr \
  --level high \
  --name smoke_lot_high \
  --output-dir "$OUT_DIR" >/dev/null

echo "[8/8] Batch import: SMI and CSV"
_SMI=$(mktemp /tmp/smoke_batch_XXXXXX.smi)
_CSV=$(mktemp /tmp/smoke_batch_XXXXXX.csv)
printf "CCO ethanol\nCC(=O)O acetic_acid\nc1ccccc1 benzene\n" > "$_SMI"
printf "smiles,name\nCCN,ethylamine\nCCCC,butane\n" > "$_CSV"

EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --batch-smi "$_SMI" \
  --workflow nmr --level low \
  --output-dir "$OUT_DIR" >/dev/null

EASYSPECTRA_XTB=__none__ ./build/easynmr \
  --batch-csv "$_CSV" \
  --workflow nmr --level low \
  --output-dir "$OUT_DIR" >/dev/null

rm -f "$_SMI" "$_CSV"

echo "Smoke checks completed."
