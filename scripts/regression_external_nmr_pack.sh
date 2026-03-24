#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

BUILD_DIR="${BUILD_DIR:-build}"
PACK_DIR="${PACK_DIR:-examples/external_nmr_pack}"
EXPCHECK="./${BUILD_DIR}/easynmr-expcheck"

if [[ ! -x "$EXPCHECK" ]]; then
  echo "[setup] building easynmr-expcheck"
  cmake --build "$BUILD_DIR" -j >/dev/null
fi

if [[ ! -x "$EXPCHECK" ]]; then
  echo "error: missing executable $EXPCHECK"
  exit 1
fi

fail() {
  echo "error: $*" >&2
  exit 1
}

find_bruker_root() {
  local extracted_root="$1"
  while IFS= read -r dir; do
    if [[ -f "$dir/acqus" && ( -f "$dir/fid" || -f "$dir/ser" || -d "$dir/pdata" ) ]]; then
      echo "$dir"
      return 0
    fi
  done < <(find "$extracted_root" -type d)
  return 1
}

total_converted=0

for tier in easy medium hard; do
  case "$tier" in
    easy)
      mnova_name="20-1H.mnova"
      bruker_name="20-1H.zip"
      ;;
    medium)
      mnova_name="20-13C.mnova"
      bruker_name="20-13C.zip"
      ;;
    hard)
      mnova_name="20-hsqc.mnova"
      bruker_name="20-hsqc.zip"
      ;;
    *)
      fail "unknown tier $tier"
      ;;
  esac

  raw_dir="$PACK_DIR/$tier/raw"
  conv_dir="$PACK_DIR/$tier/converted"
  mnova="$raw_dir/$mnova_name"
  bruker_zip="$raw_dir/$bruker_name"

  echo "[tier:$tier] checking raw files"
  [[ -f "$mnova" ]] || fail "missing $mnova"
  [[ -f "$bruker_zip" ]] || fail "missing $bruker_zip"
  mkdir -p "$conv_dir"

  echo "[tier:$tier] expecting MNova direct import failure"
  tmp_out="$(mktemp)"
  if "$EXPCHECK" "$mnova" >"$tmp_out" 2>&1; then
    cat "$tmp_out"
    rm -f "$tmp_out"
    fail "MNova file unexpectedly parsed: $mnova"
  fi
  rg -q "MNova project files are not parsed directly" "$tmp_out" || {
    cat "$tmp_out"
    rm -f "$tmp_out"
    fail "MNova failure message changed for $mnova"
  }
  rm -f "$tmp_out"

  echo "[tier:$tier] expecting Bruker directory import failure"
  tmp_dir="$(mktemp -d)"
  unzip -q "$bruker_zip" -d "$tmp_dir"
  bruker_dir="$(find_bruker_root "$tmp_dir" || true)"
  [[ -n "$bruker_dir" ]] || {
    rm -rf "$tmp_dir"
    fail "could not locate extracted Bruker root inside $bruker_zip"
  }

  tmp_out="$(mktemp)"
  if "$EXPCHECK" "$bruker_dir" >"$tmp_out" 2>&1; then
    cat "$tmp_out"
    rm -f "$tmp_out"
    rm -rf "$tmp_dir"
    fail "Bruker directory unexpectedly parsed: $bruker_dir"
  fi
  rg -q "Bruker raw directories are not parsed directly" "$tmp_out" || {
    cat "$tmp_out"
    rm -f "$tmp_out"
    rm -rf "$tmp_dir"
    fail "Bruker failure message changed for $bruker_dir"
  }
  rm -f "$tmp_out"
  rm -rf "$tmp_dir"

  echo "[tier:$tier] checking converted exports (if present)"
  shopt -s nullglob
  converted=("$conv_dir"/*.csv "$conv_dir"/*.txt "$conv_dir"/*.asc "$conv_dir"/*.dat)
  shopt -u nullglob

  if (( ${#converted[@]} == 0 )); then
    echo "[tier:$tier] no converted files yet (expected until first export)"
    continue
  fi

  for f in "${converted[@]}"; do
    if ! "$EXPCHECK" "$f" >/tmp/easynmr_expcheck_converted.log 2>&1; then
      cat /tmp/easynmr_expcheck_converted.log
      fail "converted file failed parse: $f"
    fi
    ((total_converted+=1))
  done

done

echo "[done] external pack regression checks passed"
echo "[done] converted files validated: $total_converted"
