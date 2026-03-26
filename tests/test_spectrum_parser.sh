#!/usr/bin/env bash
# test_spectrum_parser.sh — Edge-case tests for the experimental spectrum parser.
#
# Tests the C++ load_experimental_spectrum() function through the easynmr-expcheck
# binary.  Creates temporary input files that cover boundary conditions and
# error paths that the integration matrix does not exercise.
#
# Usage:
#   bash tests/test_spectrum_parser.sh [path/to/build]
#
# The build directory defaults to "build" relative to the repo root.
# Exit code 0 = all tests passed.  Non-zero = at least one failure.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

BUILD_DIR="${1:-build}"
EXPCHECK="$BUILD_DIR/easynmr-expcheck"

# ── Guard: binary must be available ────────────────────────────────────────
if [[ ! -x "$EXPCHECK" ]]; then
    echo "SKIP: easynmr-expcheck not found at $EXPCHECK (run cmake --build $BUILD_DIR first)"
    exit 0
fi

TMPDIR_TESTS="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_TESTS"' EXIT

pass=0
fail=0

# ── Helpers ─────────────────────────────────────────────────────────────────

assert_expcheck_succeeds() {
    local label="$1"
    local file="$2"
    local output
    output="$("$EXPCHECK" "$file" 2>&1)"
    local rc=$?
    if [[ $rc -eq 0 ]]; then
        echo "  PASS: $label"
        pass=$((pass + 1))
    else
        echo "  FAIL: $label (expected success, got exit $rc)"
        echo "        output: $output"
        fail=$((fail + 1))
    fi
}

assert_expcheck_fails() {
    local label="$1"
    local file="$2"
    local expected_fragment="${3:-}"
    local output
    output="$("$EXPCHECK" "$file" 2>&1)" || true
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        if [[ -n "$expected_fragment" ]] && ! echo "$output" | grep -qi "$expected_fragment"; then
            echo "  FAIL: $label (failed as expected but missing fragment '$expected_fragment')"
            echo "        output: $output"
            fail=$((fail + 1))
        else
            echo "  PASS: $label"
            pass=$((pass + 1))
        fi
    else
        echo "  FAIL: $label (expected failure, but exited 0)"
        echo "        output: $output"
        fail=$((fail + 1))
    fi
}

assert_expcheck_output_contains() {
    local label="$1"
    local file="$2"
    local fragment="$3"
    local output
    output="$("$EXPCHECK" "$file" 2>&1)"
    local rc=$?
    if [[ $rc -eq 0 ]] && echo "$output" | grep -qi "$fragment"; then
        echo "  PASS: $label"
        pass=$((pass + 1))
    else
        echo "  FAIL: $label (exit=$rc, expected fragment '$fragment' in output)"
        echo "        output: $output"
        fail=$((fail + 1))
    fi
}

# ── Generate a generic CSV with N data rows ──────────────────────────────────

make_csv() {
    local path="$1"
    local n_rows="$2"
    printf "ppm,intensity\n" > "$path"
    for i in $(seq 1 "$n_rows"); do
        printf "%.6f,%.8f\n" "$(echo "scale=6; $i * 0.01" | bc)" "0.00100000" >> "$path"
    done
}

# ════════════════════════════════════════════════════════════════════════════
echo "[1] Boundary: minimum accepted point count (8 rows)"
F="$TMPDIR_TESTS/exactly_8.csv"
make_csv "$F" 8
assert_expcheck_succeeds "exactly 8 data rows accepted" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[2] Boundary: one below minimum (7 rows)"
F="$TMPDIR_TESTS/only_7.csv"
make_csv "$F" 7
assert_expcheck_fails "7 data rows rejected as too few" "$F" "too few"

# ════════════════════════════════════════════════════════════════════════════
echo "[3] Boundary: well above minimum (3000 rows)"
F="$TMPDIR_TESTS/large.csv"
make_csv "$F" 3000
assert_expcheck_succeeds "3000 data rows accepted" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[4] Empty file"
F="$TMPDIR_TESTS/empty.csv"
> "$F"
assert_expcheck_fails "empty file rejected" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[5] File with only header/comments, no data"
F="$TMPDIR_TESTS/header_only.csv"
cat > "$F" <<'EOF'
# This is a comment line
// Another comment
ppm,intensity
## JCAMP header
EOF
assert_expcheck_fails "header-only file rejected as too few points" "$F" "too few"

# ════════════════════════════════════════════════════════════════════════════
echo "[6] MNova semicolon-separated format (European decimal comma)"
F="$TMPDIR_TESTS/mnova_format.txt"
{
    printf "# MestreNova Text Export\n"
    printf "# Source: mnova_test\n"
    for i in $(seq 1 20); do
        # Semicolon separator, comma decimal: e.g. "8,250000;0,001200"
        printf "%d,250000;0,001200\n" "$i"
    done
} > "$F"
assert_expcheck_succeeds "MNova semicolon/comma format accepted" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[7] Bruker-style space-separated format"
F="$TMPDIR_TESTS/bruker_format.txt"
{
    printf "##TITLE= Bruker TopSpin ASCII Export\n"
    printf "##JCAMPDX=5.00\n"
    printf "# X_PPM Y_INTENSITY\n"
    for i in $(seq 1 20); do
        printf "%.6f 0.00100000\n" "$(echo "scale=6; $i * 0.1" | bc)"
    done
} > "$F"
assert_expcheck_succeeds "Bruker space-separated format accepted" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[8] Lines with text mixed in — valid lines should still parse"
F="$TMPDIR_TESTS/mixed.csv"
{
    printf "ppm,intensity\n"
    printf "1.000,0.001\n"
    printf "this line has no numbers\n"          # should be skipped
    printf "2.000,0.002\n"
    printf "bad_line_xyz\n"                      # should be skipped
    for i in $(seq 3 12); do
        printf "%d.000,0.001\n" "$i"
    done
} > "$F"
assert_expcheck_succeeds "mixed valid/invalid lines accepted if enough valid rows" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[9] .mnova extension rejected with informative message"
F="$TMPDIR_TESTS/spectrum.mnova"
printf "dummy content\n" > "$F"
assert_expcheck_fails ".mnova extension rejected" "$F" "MNova project"

# ════════════════════════════════════════════════════════════════════════════
echo "[10] .mnv extension rejected with informative message"
F="$TMPDIR_TESTS/spectrum.mnv"
printf "dummy content\n" > "$F"
assert_expcheck_fails ".mnv extension rejected" "$F" "MNova"

# ════════════════════════════════════════════════════════════════════════════
echo "[11] Non-existent file path"
assert_expcheck_fails "non-existent path rejected" "/tmp/easyspectra_does_not_exist_xyz.csv" "not exist"

# ════════════════════════════════════════════════════════════════════════════
echo "[12] Bruker raw directory layout detected and rejected"
BRUKER_DIR="$TMPDIR_TESTS/bruker_raw"
mkdir -p "$BRUKER_DIR"
touch "$BRUKER_DIR/acqus" "$BRUKER_DIR/fid"
assert_expcheck_fails "Bruker raw directory rejected" "$BRUKER_DIR" "Bruker raw"

# ════════════════════════════════════════════════════════════════════════════
echo "[13] File with all-comment lines plus exactly 8 valid rows"
F="$TMPDIR_TESTS/comments_then_data.csv"
{
    printf "# Header comment\n"
    printf "# Second comment\n"
    printf "ppm,intensity\n"
    for i in $(seq 1 8); do
        printf "%.3f,0.00100000\n" "$(echo "scale=3; $i * 0.5" | bc)"
    done
} > "$F"
assert_expcheck_succeeds "comments plus 8 data rows accepted" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[14] Tab-separated values (generic parser should handle)"
F="$TMPDIR_TESTS/tab_separated.txt"
{
    for i in $(seq 1 20); do
        printf "%.4f\t0.00100000\n" "$(echo "scale=4; $i * 0.1" | bc)"
    done
} > "$F"
assert_expcheck_succeeds "tab-separated values accepted by generic parser" "$F"

# ════════════════════════════════════════════════════════════════════════════
echo "[15] Output fields present on success"
F="$TMPDIR_TESTS/check_output.csv"
make_csv "$F" 50
assert_expcheck_output_contains "success output contains 'ok'" "$F" "ok"
assert_expcheck_output_contains "success output contains 'format'" "$F" "format"
assert_expcheck_output_contains "success output contains 'points'" "$F" "points"

# ════════════════════════════════════════════════════════════════════════════
echo ""
echo "Parser edge-case tests complete: $pass passed, $fail failed."

if [[ $fail -gt 0 ]]; then
    exit 1
fi
exit 0
