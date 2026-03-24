#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CASES_CSV="tests/spectra_comparison_cases.csv"
OUT_ROOT="${1:-output_test_matrix}"
CASE_FILTER="${2:-}"
mkdir -p "$OUT_ROOT"

if [[ ! -f "$CASES_CSV" ]]; then
  echo "error: missing $CASES_CSV" >&2
  exit 1
fi

echo "[1/3] Build"
cmake --build build -j >/dev/null

echo "[2/3] Validate all experimental overlay files with easynmr-expcheck"
expcheck_count=0
while IFS=, read -r difficulty case_name workflow target_product nucleus_arg smiles computed_ref exp_file format_hint notes; do
  if [[ -n "$CASE_FILTER" && "$case_name" != *"$CASE_FILTER"* ]]; then
    continue
  fi
  if [[ ! -f "$exp_file" ]]; then
    echo "error: missing experimental file for $case_name: $exp_file" >&2
    exit 1
  fi
  ./build/easynmr-expcheck "$exp_file" >/dev/null
  echo "  ok expcheck: $case_name -> $exp_file"
  expcheck_count=$((expcheck_count + 1))
done < <(tail -n +2 "$CASES_CSV")

if [[ $expcheck_count -eq 0 ]]; then
  echo "error: no cases matched filter '${CASE_FILTER}'" >&2
  exit 1
fi

echo "[3/3] Run predicted spectra matrix and verify expected product label"
pass_count=0

while IFS=, read -r difficulty case_name workflow target_product nucleus_arg smiles computed_ref exp_file format_hint notes; do
  if [[ -n "$CASE_FILTER" && "$case_name" != *"$CASE_FILTER"* ]]; then
    continue
  fi
  run_cmd=(
    ./build/easynmr
    --input "$smiles"
    --input-format smiles
    --workflow "$workflow"
    --name "$case_name"
    --output-dir "$OUT_ROOT"
  )

  if [[ "$workflow" == "nmr" ]]; then
    if [[ -n "$nucleus_arg" ]]; then
      run_cmd+=(--nucleus "$nucleus_arg")
    fi
  fi

  run_output="$("${run_cmd[@]}")"
  run_dir="$(printf '%s\n' "$run_output" | awk -F': ' '/^Output dir:/{print $2; exit}')"
  if [[ -z "$run_dir" ]]; then
    echo "error: could not parse output dir for $case_name" >&2
    echo "$run_output" >&2
    exit 1
  fi

  manifest="$run_dir/spectra_manifest.csv"
  if [[ ! -f "$manifest" ]]; then
    echo "error: missing manifest for $case_name: $manifest" >&2
    exit 1
  fi

  if ! awk -F, -v want="$target_product" 'NR>1 && $1==want {found=1} END{exit(found?0:1)}' "$manifest"; then
    echo "error: expected product '$target_product' not found for $case_name in $manifest" >&2
    cat "$manifest" >&2
    exit 1
  fi

  if [[ ! -f "$computed_ref" ]]; then
    echo "warning: computed reference file not found on disk for $case_name: $computed_ref"
  fi

  echo "  ok compare: $case_name -> workflow=$workflow product=$target_product"
  pass_count=$((pass_count + 1))
done < <(tail -n +2 "$CASES_CSV")

echo "Comparison matrix checks completed (${pass_count} cases)."

echo "Tip: open tests/spectra_comparison_cases.csv for exact SMILES/file mapping."
