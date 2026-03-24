#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_DIR="${1:-output_smoke_products}"
mkdir -p "$OUT_DIR"

echo "[1/5] Build"
cmake --build build -j >/dev/null

echo "[2/5] NMR workflow smoke run"
EASYNMR_XTB=__none__ ./build/easynmr \
  --input "CCO" \
  --input-format smiles \
  --workflow nmr \
  --name smoke_nmr \
  --output-dir "$OUT_DIR" >/dev/null

echo "[3/5] CD workflow smoke run"
EASYNMR_XTB=__none__ ./build/easynmr \
  --input "CCO" \
  --input-format smiles \
  --workflow cd \
  --name smoke_cd \
  --output-dir "$OUT_DIR" >/dev/null

echo "[4/5] 31P workflow smoke run"
EASYNMR_XTB=__none__ ./build/easynmr \
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

echo "Smoke checks completed."
