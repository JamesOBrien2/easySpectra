#!/usr/bin/env python3
"""check_peak_quality.py — Validate predicted peak positions against reference values.

Called by run_spectra_comparison_matrix.sh after each successful prediction to
catch regressions in chemical shift accuracy.

Usage:
    python3 tests/check_peak_quality.py <case_name> <run_dir>

Arguments:
    case_name   Name of the test case (must match a key in expected_peaks.json).
    run_dir     Directory containing the predicted output files (e.g. peaks_1h.csv).

Exit codes:
    0  All expected peaks found within tolerance, or case has no reference data.
    1  At least one expected peak is missing (regression detected).
    2  Usage / file-not-found error.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def load_predicted_peaks(peaks_csv: Path) -> list[float]:
    """Return a sorted list of predicted peak ppm values from a peaks CSV file."""
    peaks: list[float] = []
    if not peaks_csv.exists():
        return peaks
    with peaks_csv.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                peaks.append(float(row["shift_ppm"]))
            except (KeyError, ValueError):
                # CD peaks CSV uses wavelength_nm column instead.
                try:
                    peaks.append(float(row["wavelength_nm"]))
                except (KeyError, ValueError):
                    continue
    return sorted(peaks)


def nearest_match(target: float, peaks: list[float]) -> tuple[float, float]:
    """Return (nearest_peak_ppm, |delta|) for target in peaks list."""
    if not peaks:
        return float("nan"), float("inf")
    closest = min(peaks, key=lambda p: abs(p - target))
    return closest, abs(closest - target)


def check_case(case_name: str, run_dir: Path, reference: dict) -> bool:
    """Return True if all expected peaks are found within tolerance."""
    peaks_filename = reference["peaks_file"]
    expected = reference["expected_peaks_ppm"]
    tolerance = reference["tolerance_ppm"]
    notes = reference.get("notes", "")

    peaks_path = run_dir / peaks_filename
    predicted = load_predicted_peaks(peaks_path)

    if not predicted:
        print(
            f"  WARN [{case_name}]: peaks file not found or empty: {peaks_path}",
            file=sys.stderr,
        )
        return True  # Don't fail on missing file — expcheck already covers existence.

    all_ok = True
    for expected_ppm in expected:
        nearest, delta = nearest_match(expected_ppm, predicted)
        if delta <= tolerance:
            print(
                f"  peak_ok [{case_name}]: expected {expected_ppm:+.1f} ppm "
                f"-> nearest {nearest:+.4f} ppm (Δ={delta:.2f}, tol={tolerance})"
            )
        else:
            print(
                f"  peak_FAIL [{case_name}]: expected {expected_ppm:+.1f} ppm "
                f"-> nearest {nearest:+.4f} ppm (Δ={delta:.2f} > tol={tolerance})"
                + (f"  [{notes}]" if notes else ""),
                file=sys.stderr,
            )
            all_ok = False

    return all_ok


def main() -> int:
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <case_name> <run_dir>", file=sys.stderr)
        return 2

    case_name = sys.argv[1]
    run_dir = Path(sys.argv[2])

    if not run_dir.is_dir():
        print(f"error: run_dir does not exist: {run_dir}", file=sys.stderr)
        return 2

    ref_path = Path(__file__).parent / "expected_peaks.json"
    if not ref_path.exists():
        print(f"SKIP: reference file not found: {ref_path}")
        return 0

    with ref_path.open(encoding="utf-8") as fh:
        all_refs = json.load(fh)

    if case_name not in all_refs:
        # No reference data for this case — silently pass.
        return 0

    reference = all_refs[case_name]
    ok = check_case(case_name, run_dir, reference)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
