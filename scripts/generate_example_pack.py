#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass
class ExampleCase:
    difficulty: str
    case_name: str
    workflow: str
    smiles: str
    nucleus: str
    label: str
    computed_file: str
    experimental_file: str


CASES: List[ExampleCase] = [
    ExampleCase(
        difficulty="easy",
        case_name="easy_cd_lactic",
        workflow="cd",
        smiles="C[C@H](O)C(=O)O",
        nucleus="",
        label="CD",
        computed_file="examples/computed/easy_cd_computed.csv",
        experimental_file="examples/experimental/easy_cd_overlay.csv",
    ),
    ExampleCase(
        difficulty="medium",
        case_name="medium_cd_menthol",
        workflow="cd",
        smiles="CC[C@H]1CC[C@@H](C)C(C)C1O",
        nucleus="",
        label="CD",
        computed_file="examples/computed/medium_cd_computed.csv",
        experimental_file="examples/experimental/medium_cd_overlay.csv",
    ),
    ExampleCase(
        difficulty="hard",
        case_name="hard_cd_ephedrine",
        workflow="cd",
        smiles="CN[C@@H](C)[C@H](O)c1ccccc1",
        nucleus="",
        label="CD",
        computed_file="examples/computed/hard_cd_computed.csv",
        experimental_file="examples/experimental/hard_cd_overlay.csv",
    ),
    ExampleCase(
        difficulty="easy",
        case_name="easy_13c_ethanol",
        workflow="nmr",
        smiles="CCO",
        nucleus="13C",
        label="13C",
        computed_file="examples/computed/easy_13c_computed.csv",
        experimental_file="examples/experimental/easy_13c_overlay.csv",
    ),
    ExampleCase(
        difficulty="medium",
        case_name="medium_13c_isobutylbenzene",
        workflow="nmr",
        smiles="CC(C)Cc1ccccc1",
        nucleus="13C",
        label="13C",
        computed_file="examples/computed/medium_13c_computed.csv",
        experimental_file="examples/experimental/medium_13c_overlay.csv",
    ),
    ExampleCase(
        difficulty="hard",
        case_name="hard_13c_carbamate_aryl",
        workflow="nmr",
        smiles="CCOC(=O)N1CCC(CC1)C2=CC=CC=C2",
        nucleus="13C",
        label="13C",
        computed_file="examples/computed/hard_13c_computed.csv",
        experimental_file="examples/experimental/hard_13c_overlay.csv",
    ),
    ExampleCase(
        difficulty="easy",
        case_name="easy_19f_fluoroethane",
        workflow="nmr",
        smiles="CCF",
        nucleus="19F",
        label="19F",
        computed_file="examples/computed/easy_19f_computed.csv",
        experimental_file="examples/experimental/easy_19f_overlay.csv",
    ),
    ExampleCase(
        difficulty="medium",
        case_name="medium_19f_fluorobenzene",
        workflow="nmr",
        smiles="Fc1ccccc1",
        nucleus="19F",
        label="19F",
        computed_file="examples/computed/medium_19f_computed.csv",
        experimental_file="examples/experimental/medium_19f_overlay.csv",
    ),
    ExampleCase(
        difficulty="hard",
        case_name="hard_19f_trifluoroanisole",
        workflow="nmr",
        smiles="COc1ccc(C(F)(F)F)cc1",
        nucleus="19F",
        label="19F",
        computed_file="examples/computed/hard_19f_computed.csv",
        experimental_file="examples/experimental/hard_19f_overlay.csv",
    ),
    ExampleCase(
        difficulty="easy",
        case_name="easy_31p_trimethylphosphine_oxide",
        workflow="nmr",
        smiles="CP(=O)(C)C",
        nucleus="31P",
        label="31P",
        computed_file="examples/computed/easy_31p_computed.csv",
        experimental_file="examples/experimental/easy_31p_overlay.csv",
    ),
    ExampleCase(
        difficulty="medium",
        case_name="medium_31p_triethyl_phosphate",
        workflow="nmr",
        smiles="CCOP(=O)(OCC)OCC",
        nucleus="31P",
        label="31P",
        computed_file="examples/computed/medium_31p_computed.csv",
        experimental_file="examples/experimental/medium_31p_overlay.csv",
    ),
    ExampleCase(
        difficulty="hard",
        case_name="hard_31p_triphenylphosphine_oxide",
        workflow="nmr",
        smiles="O=P(c1ccccc1)(c1ccccc1)c1ccccc1",
        nucleus="31P",
        label="31P",
        computed_file="examples/computed/hard_31p_computed.csv",
        experimental_file="examples/experimental/hard_31p_overlay.csv",
    ),
]


LEGACY_BENCHMARK_ROWS = [
    [
        "easy",
        "ethanol_nmr",
        "nmr",
        "CCO",
        "examples/computed/easy_nmr_computed.csv",
        "examples/experimental/easy_mnova_export.txt",
        "mnova_text_export",
    ],
    [
        "medium",
        "isobutylbenzene_nmr",
        "nmr",
        "CC(C)Cc1ccccc1",
        "examples/computed/medium_nmr_computed.csv",
        "examples/experimental/medium_bruker_export.txt",
        "bruker_text_export",
    ],
    [
        "hard",
        "carbamate_aryl_nmr",
        "nmr",
        "CCOC(=O)N1CCC(CC1)C2=CC=CC=C2",
        "examples/computed/hard_nmr_computed.csv",
        "examples/experimental/hard_generic_export.csv",
        "generic_text_export",
    ],
]


def _run(cmd: List[str], cwd: Path, env: dict) -> str:
    proc = subprocess.run(cmd, cwd=str(cwd), env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}\n"
        )
    return proc.stdout


def _extract_output_dir(cli_stdout: str) -> str:
    for line in cli_stdout.splitlines():
        if line.startswith("Output dir:"):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"Could not parse output dir from CLI output:\n{cli_stdout}")


def _read_xy(csv_path: Path) -> List[Tuple[float, float]]:
    points: List[Tuple[float, float]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        _ = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                x = float(row[0].strip())
                y = float(row[1].strip())
            except Exception:
                continue
            points.append((x, y))
    return points


def _write_xy(csv_path: Path, points: List[Tuple[float, float]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "intensity"])
        for x, y in points:
            writer.writerow([f"{x:.6f}", f"{y:.8f}"])


def _difficulty_params(difficulty: str, workflow: str) -> Tuple[float, float, float]:
    if workflow == "cd":
        if difficulty == "easy":
            return (0.6, 0.06, 0.010)
        if difficulty == "medium":
            return (1.2, 0.10, 0.018)
        return (2.0, 0.16, 0.028)

    if difficulty == "easy":
        return (0.010, 0.0012, 0.010)
    if difficulty == "medium":
        return (0.028, 0.0028, 0.020)
    return (0.052, 0.0045, 0.035)


def _make_experimental_like(case: ExampleCase, computed_points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    if not computed_points:
        return []

    stride = max(1, len(computed_points) // 2800)
    base_points = computed_points[::stride]
    if len(base_points) < 400:
        base_points = computed_points

    seed = int(hashlib.sha1(case.case_name.encode("utf-8")).hexdigest()[:16], 16)
    rng = random.Random(seed)
    shift_abs, jitter_x, noise_y = _difficulty_params(case.difficulty, case.workflow)
    direction = -1.0 if (seed % 2) else 1.0
    x_shift = direction * shift_abs

    transformed: List[Tuple[float, float]] = []
    for x, y in base_points:
        x_new = x + x_shift + rng.uniform(-jitter_x, jitter_x)
        y_scale = 0.94 + 0.10 * rng.random()
        y_new = y * y_scale + rng.gauss(0.0, noise_y)
        if case.workflow != "cd":
            y_new = max(0.0, y_new)
        else:
            y_new = max(-1.3, min(1.3, y_new))
        transformed.append((x_new, y_new))

    return transformed


def _copy_case_outputs(root: Path, case: ExampleCase, job_output_dir: Path) -> None:
    manifest = job_output_dir / "spectra_manifest.csv"
    if not manifest.exists():
        raise RuntimeError(f"Missing manifest for {case.case_name}: {manifest}")

    target_spectrum: Path | None = None
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            label = (row.get("spectrum_label") or "").strip()
            if label == case.label:
                raw = (row.get("spectrum_csv") or "").strip()
                p = Path(raw)
                target_spectrum = p if p.is_absolute() else (root / p)
                break
    if target_spectrum is None or not target_spectrum.exists():
        raise RuntimeError(f"Could not locate spectrum for label={case.label} in {manifest}")

    computed_target = root / case.computed_file
    computed_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(target_spectrum, computed_target)

    points = _read_xy(computed_target)
    exp_points = _make_experimental_like(case, points)
    _write_xy(root / case.experimental_file, exp_points)


def _write_benchmark_csv(root: Path) -> None:
    path = root / "examples" / "benchmark_cases.csv"
    rows: List[List[str]] = []
    rows.extend(LEGACY_BENCHMARK_ROWS)
    for case in CASES:
        rows.append(
            [
                case.difficulty,
                case.case_name,
                case.workflow,
                case.smiles,
                case.computed_file,
                case.experimental_file,
                "generic_text_export",
            ]
        )

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "difficulty",
                "case_name",
                "workflow",
                "smiles",
                "computed_reference_csv",
                "experimental_overlay_file",
                "experimental_format_hint",
            ]
        )
        writer.writerows(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    cli = root / "build" / "easynmr"
    if not cli.exists():
        raise RuntimeError("Missing build/easynmr. Build the project first (cmake --build build).")

    generated_out = root / "output_example_generation"
    generated_out.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["EASYNMR_XTB"] = "__none__"

    print(f"[generate] using CLI: {cli}")
    for idx, case in enumerate(CASES, start=1):
        cmd = [
            str(cli),
            "--input",
            case.smiles,
            "--input-format",
            "smiles",
            "--workflow",
            case.workflow,
            "--name",
            case.case_name,
            "--output-dir",
            str(generated_out),
        ]
        if case.workflow == "nmr":
            cmd.extend(["--nucleus", case.nucleus])

        print(f"[{idx:02d}/{len(CASES)}] {case.case_name} ({case.workflow} {case.label})")
        stdout = _run(cmd, cwd=root, env=env)
        out_dir = _extract_output_dir(stdout)
        job_output_dir = Path(out_dir)
        if not job_output_dir.is_absolute():
            job_output_dir = root / job_output_dir
        _copy_case_outputs(root, case, job_output_dir)

    _write_benchmark_csv(root)
    print("[done] wrote computed + experimental examples and updated examples/benchmark_cases.csv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
