#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdDetermineBonds, rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D

KCAL_PER_EH = 627.509474
GAS_R_KCAL = 0.0019872041
TEMP_K = 298.15
MAX_PARALLEL_XTB = 4
# GFN2-xTB harmonic frequency scaling factor (Merz et al., JCTC 2021)
IR_SCALE_FACTOR_GFN2: float = 0.9740
IR_MAX_HESS_CONFS: int = 10
UVVIS_MAX_STDA_CONFS: int = 10
UVVIS_BROADENING_NM: float = 15.0   # Gaussian FWHM in nm

SOLVENT_ALIASES = {
    "cdcl3": "chcl3",
    "chcl3": "chcl3",
    "chloroform": "chcl3",
    "dmso": "dmso",
    "h2o": "h2o",
    "water": "h2o",
}

MULTIPLICITY_LABELS = {
    1: "singlet",
    2: "doublet",
    3: "triplet",
    4: "quartet",
    5: "quintet",
}

# Tolerance used when matching reactant peaks to product peaks per nucleus (ppm).
COMPARE_TOLERANCE_PPM: Dict[str, float] = {
    "1h": 0.30,
    "13c": 2.00,
    "19f": 5.00,
    "31p": 5.00,
}

NUCLEUS_ALIASES = {
    "auto": "auto",
    "all": "auto",
    "1h": "1h",
    "h1": "1h",
    "h": "1h",
    "13c": "13c",
    "c13": "13c",
    "c": "13c",
    "19f": "19f",
    "f19": "19f",
    "f": "19f",
    "31p": "31p",
    "p31": "31p",
    "p": "31p",
}

NUCLEUS_SYMBOL = {
    "1h": "H",
    "13c": "C",
    "19f": "F",
    "31p": "P",
}

NUCLEUS_LABEL = {
    "1h": "1H",
    "13c": "13C",
    "19f": "19F",
    "31p": "31P",
}


@dataclass
class ConformerResult:
    conf_id: int
    energy_eh: float
    method: str


@dataclass
class GroupPrediction:
    group_id: int
    atom_indices: List[int]
    shift_ppm: float
    shift_hz: float
    multiplicity: str
    j_hz: float
    integral: float


def write_progress(progress_path: Optional[Path], stage: str, message: str, fraction: float) -> None:
    if progress_path is None:
        return
    try:
        payload = {
            "stage": stage,
            "message": message,
            "fraction": max(0.0, min(1.0, float(fraction))),
            "timestamp_unix": int(time.time()),
        }
        tmp = progress_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        tmp.replace(progress_path)
    except Exception:
        # Progress is best-effort and should never abort the main workflow.
        return


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="easySpectra backend worker")
    parser.add_argument("--request", required=True, help="Request JSON path")
    parser.add_argument("--response", required=True, help="Response JSON path")
    return parser.parse_args()


def hash_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def parse_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except Exception:
        return None


def parse_bool(value: object, fallback: bool) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"1", "true", "yes", "y", "on"}:
            return True
        if text in {"0", "false", "no", "n", "off"}:
            return False
    return fallback


def normalize_nucleus(text: str) -> str:
    return NUCLEUS_ALIASES.get(text.strip().lower(), "1h")


def maybe_read_text_from_path(raw: str) -> Optional[str]:
    candidate = raw.strip()
    if not candidate:
        return None
    if "\n" in raw or "\r" in raw:
        return None
    if len(candidate) > 1024:
        return None
    try:
        path = Path(candidate)
        if path.exists() and path.is_file():
            return path.read_text(encoding="utf-8")
    except Exception:
        return None
    return None


def normalize_xyz_text(raw: str) -> Tuple[Optional[str], bool]:
    text = raw.strip()
    if not text:
        return None, False

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) < 1:
        return None, False

    # Already canonical XYZ: atom count line + comment + atom rows.
    first_parts = lines[0].split()
    try:
        declared = int(first_parts[0]) if len(first_parts) == 1 else -1
        if declared > 0 and len(lines) >= declared + 2:
            canonical = "\n".join(lines[: declared + 2]) + "\n"
            return canonical, False
    except Exception:
        pass

    parsed_rows: List[Tuple[str, float, float, float]] = []
    symbol_re = re.compile(r"^[A-Za-z]{1,3}$")
    for line in lines:
        parts = line.split()
        if len(parts) >= 4 and symbol_re.match(parts[0]):
            x = parse_float(parts[1])
            y = parse_float(parts[2])
            z = parse_float(parts[3])
            if x is None or y is None or z is None:
                return None, False
            parsed_rows.append((parts[0], x, y, z))
            continue
        if len(parts) >= 5 and parse_float(parts[0]) is not None and symbol_re.match(parts[1]):
            x = parse_float(parts[2])
            y = parse_float(parts[3])
            z = parse_float(parts[4])
            if x is None or y is None or z is None:
                return None, False
            parsed_rows.append((parts[1], x, y, z))
            continue
        return None, False

    if not parsed_rows:
        return None, False

    normalized = [str(len(parsed_rows)), "easySpectra normalized XYZ block"]
    for sym, x, y, z in parsed_rows:
        normalized.append(f"{sym} {x:.10f} {y:.10f} {z:.10f}")
    return "\n".join(normalized) + "\n", True


def nuclei_present(mol: Chem.Mol) -> List[str]:
    targets = [("1h", 1), ("13c", 6), ("19f", 9), ("31p", 15)]
    out: List[str] = []
    for nucleus, atomic_number in targets:
        if any(atom.GetAtomicNum() == atomic_number for atom in mol.GetAtoms()):
            out.append(nucleus)
    return out


def read_xyz(path: Path) -> Tuple[List[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid XYZ file: {path}")

    n_atoms = int(lines[0].strip())
    symbols: List[str] = []
    coords: List[List[float]] = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if len(symbols) != n_atoms:
        raise ValueError(f"Unexpected atom count in XYZ: {path}")

    return symbols, np.array(coords, dtype=float)


def write_conf_xyz(mol: Chem.Mol, conf_id: int, path: Path, title: str) -> None:
    conf = mol.GetConformer(conf_id)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{mol.GetNumAtoms()}\n")
        handle.write(f"{title}\n")
        for atom_idx, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            handle.write(f"{atom.GetSymbol()} {pos.x:.10f} {pos.y:.10f} {pos.z:.10f}\n")


def write_editable_xyz(mol: Chem.Mol, path: Path, warnings: List[str]) -> None:
    editable = Chem.Mol(mol)
    conf_id = 0

    if editable.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xE11A
        conf_id = AllChem.EmbedMolecule(editable, params)

        if conf_id < 0:
            try:
                rdDepictor.SetPreferCoordGen(True)
                rdDepictor.Compute2DCoords(editable)
                conf_id = 0
                warnings.append("3D embedding failed for editor export; wrote 2D XYZ fallback.")
            except Exception:
                raise ValueError("Could not generate coordinates for editable XYZ export.")
    else:
        conf_id = editable.GetConformer().GetId()

    write_conf_xyz(editable, conf_id, path, "easySpectra editable structure")


def load_molecule(input_format: str, value: str, warnings: List[str]) -> Chem.Mol:
    mol: Optional[Chem.Mol] = None
    path_text = maybe_read_text_from_path(value)

    if input_format == "smiles":
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError("Failed to parse SMILES input")
    elif input_format in {"mol", "sdf"}:
        if path_text is not None:
            if input_format == "mol":
                mol = Chem.MolFromMolBlock(path_text, removeHs=False)
            else:
                parts = path_text.split("$$$$")
                mol = Chem.MolFromMolBlock(parts[0], removeHs=False) if parts and parts[0].strip() else None
        else:
            if input_format == "mol":
                mol = Chem.MolFromMolBlock(value, removeHs=False)
            else:
                parts = value.split("$$$$")
                mol = Chem.MolFromMolBlock(parts[0], removeHs=False) if parts and parts[0].strip() else None
        if mol is None:
            raise ValueError(f"Failed to parse {input_format} input")
    elif input_format == "xyz":
        xyz_text = value
        if path_text is not None:
            xyz_text = path_text

        normalized_xyz, was_normalized = normalize_xyz_text(xyz_text)
        if normalized_xyz is None:
            raise ValueError("Failed to parse XYZ input")
        if was_normalized:
            warnings.append("XYZ text was normalized from coordinate rows; inferred atom count/comment.")

        tmp = Path(tempfile.mkdtemp(prefix="easynmr_xyz_")) / "input.xyz"
        tmp.write_text(normalized_xyz, encoding="utf-8")
        mol = Chem.MolFromXYZFile(str(tmp))
        if mol is None:
            raise ValueError("Failed to parse XYZ input")
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except Exception:
            warnings.append("Bond perception from XYZ was uncertain; connectivity may be approximate.")
    else:
        raise ValueError(f"Unsupported input format '{input_format}'")

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        warnings.append("Molecule sanitization raised warnings; proceeding with best-effort interpretation.")

    if any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms()):
        warnings.append("Radical detected: prediction quality may be lower in current fast workflow.")

    if Chem.GetFormalCharge(mol) != 0:
        warnings.append("Charged species detected: prediction quality depends on protonation/state correctness.")

    return Chem.AddHs(mol)


def write_structure_svg(mol: Chem.Mol, path: Path) -> None:
    tagged = Chem.Mol(mol)
    for atom in tagged.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    mol2d = Chem.RemoveHs(tagged)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol2d)
    for atom in mol2d.GetAtoms():
        # Display original 1-based atom labels to match assignment numbering.
        atom.SetProp("atomNote", str(atom.GetAtomMapNum()))

    drawer = rdMolDraw2D.MolDraw2DSVG(380, 260)
    opts = drawer.drawOptions()
    opts.legendFontSize = 14
    opts.addAtomIndices = False
    drawer.DrawMolecule(mol2d)
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()
    path.write_text(svg, encoding="utf-8")


def write_structure_geometry_csv(mol: Chem.Mol, atoms_csv: Path, bonds_csv: Path) -> None:
    tagged = Chem.Mol(mol)
    for atom in tagged.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    mol2d = Chem.RemoveHs(tagged)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol2d)
    conf = mol2d.GetConformer()
    Chem.AssignStereochemistry(mol2d, cleanIt=True, force=True)
    try:
        Chem.WedgeMolBonds(mol2d, conf)
    except Exception:
        # Keep 2D geometry export resilient; fallback is plain bonds.
        pass

    with atoms_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["atom_index", "element", "x", "y", "attached_h"])
        for atom in mol2d.GetAtoms():
            orig_idx1 = atom.GetAtomMapNum()
            orig_idx0 = orig_idx1 - 1
            pos = conf.GetAtomPosition(atom.GetIdx())
            orig_atom = mol.GetAtomWithIdx(orig_idx0)
            attached_h = [nbr.GetIdx() + 1 for nbr in orig_atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
            writer.writerow(
                [
                    orig_idx1,
                    atom.GetSymbol(),
                    f"{pos.x:.6f}",
                    f"{pos.y:.6f}",
                    ";".join(str(h) for h in attached_h),
                ]
            )

    with bonds_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["atom_a", "atom_b", "order", "stereo_style", "stereo_from_atom"])
        for bond in mol2d.GetBonds():
            a = bond.GetBeginAtom().GetAtomMapNum()
            b = bond.GetEndAtom().GetAtomMapNum()
            order = int(round(float(bond.GetBondTypeAsDouble())))
            stereo_style = "none"
            stereo_from_atom = 0
            bond_dir = str(bond.GetBondDir())
            if bond_dir in ("BEGINWEDGE", "BEGINWEDGEBOND", "BEGINUPRIGHT"):
                stereo_style = "wedge"
                stereo_from_atom = a
            elif bond_dir in ("BEGINDASH", "BEGINDASHBOND", "BEGINDOWNRIGHT"):
                stereo_style = "dash"
                stereo_from_atom = a
            elif bond_dir in ("ENDUPRIGHT",):
                stereo_style = "wedge"
                stereo_from_atom = b
            elif bond_dir in ("ENDDOWNRIGHT",):
                stereo_style = "dash"
                stereo_from_atom = b

            writer.writerow([a, b, max(1, order), stereo_style, stereo_from_atom])


def embed_conformers(mol: Chem.Mol, max_conformers: int, seed: int, warnings: List[str]) -> List[int]:
    n_confs = max(1, min(max_conformers, 60))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed & 0x7FFFFFFF
    params.pruneRmsThresh = 0.2

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not conf_ids:
        conf_id = AllChem.EmbedMolecule(mol, randomSeed=params.randomSeed)
        if conf_id < 0:
            raise RuntimeError("Conformer embedding failed")
        conf_ids = [conf_id]

    try:
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=250)
    except Exception:
        try:
            AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=250)
        except Exception:
            warnings.append("Local force-field pre-optimization failed; continuing with raw embedded conformers.")

    return conf_ids


def mmff_energy(mol: Chem.Mol, conf_id: int) -> Optional[float]:
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    if ff is None:
        return None
    return float(ff.CalcEnergy()) / KCAL_PER_EH


def run_xtb_optimize(
    xtb_bin: str,
    mol: Chem.Mol,
    conf_id: int,
    solvent: str,
    output_dir: Path,
    timeout_s: int,
) -> Optional[float]:
    archive_dir = output_dir / f"xtb_conf_{conf_id}"
    archive_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=f"easynmr_xtb_conf_{conf_id}_"))
    log_path = archive_dir / "xtb.log"
    stdout = ""
    stderr = ""

    try:
        xyz_path = run_dir / "input.xyz"
        write_conf_xyz(mol, conf_id, xyz_path, f"conf-{conf_id}")

        cmd = [xtb_bin, str(xyz_path), "--gfn2", "--opt", "loose"]
        if solvent:
            cmd.extend(["--alpb", solvent])

        try:
            proc = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
        except Exception:
            return None

        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                proc.kill()
            stdout, stderr = proc.communicate()
            log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
            return None

        log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
        if proc.returncode != 0:
            return None

        matches = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", stdout)
        energy_eh = float(matches[-1]) if matches else None

        run_opt_xyz = run_dir / "xtbopt.xyz"
        if run_opt_xyz.exists():
            try:
                shutil.copy2(run_opt_xyz, archive_dir / "xtbopt.xyz")
                symbols, coords = read_xyz(run_opt_xyz)
                if len(symbols) == mol.GetNumAtoms():
                    conf = mol.GetConformer(conf_id)
                    for i, (_, xyz) in enumerate(zip(symbols, coords)):
                        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(xyz[0]), float(xyz[1]), float(xyz[2])))
            except Exception:
                pass

        return energy_eh
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def optimize_conformers(
    mol: Chem.Mol,
    conf_ids: Sequence[int],
    solvent: str,
    output_dir: Path,
    warnings: List[str],
    progress_path: Optional[Path] = None,
) -> List[ConformerResult]:
    xtb_bin = os.environ.get("EASYSPECTRA_XTB", "xtb")
    xtb_path = shutil.which(xtb_bin)

    results: List[ConformerResult] = []

    if xtb_path is None:
        warnings.append("xTB not found in PATH; using MMFF energies only.")
        total = max(1, len(conf_ids))
        for i, conf_id in enumerate(conf_ids, start=1):
            e = mmff_energy(mol, conf_id)
            if e is None:
                e = float(conf_id) * 1e-6
            results.append(ConformerResult(conf_id=conf_id, energy_eh=e, method="mmff"))
            write_progress(
                progress_path,
                "conformer_optimization",
                f"Optimizing conformers with MMFF fallback only ({i}/{total})",
                0.38 + 0.30 * (i / total),
            )
        return results

    max_workers = max(1, min(MAX_PARALLEL_XTB, len(conf_ids)))
    timeout_s = int(os.environ.get("EASYSPECTRA_XTB_TIMEOUT", "25"))
    total = max(1, len(conf_ids))
    completed = 0
    xtb_count = 0
    fallback_count = 0
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_map = {
            pool.submit(run_xtb_optimize, xtb_path, mol, conf_id, solvent, output_dir, timeout_s): conf_id for conf_id in conf_ids
        }
        for future in as_completed(future_map):
            conf_id = future_map[future]
            energy_eh = future.result()
            if energy_eh is None:
                fallback = mmff_energy(mol, conf_id)
                if fallback is None:
                    fallback = float(conf_id) * 1e-6
                results.append(ConformerResult(conf_id=conf_id, energy_eh=fallback, method="mmff-fallback"))
                fallback_count += 1
            else:
                results.append(ConformerResult(conf_id=conf_id, energy_eh=energy_eh, method="xtb"))
                xtb_count += 1
            completed += 1
            write_progress(
                progress_path,
                "conformer_optimization",
                (
                    "Optimizing conformers with xTB GFN2-xTB + ALPB("
                    f"{solvent}); completed {completed}/{total}, xTB {xtb_count}, fallback {fallback_count}"
                ),
                0.38 + 0.30 * (completed / total),
            )

    results.sort(key=lambda item: item.energy_eh)
    return results


def boltzmann_weights(energies_kcal: Sequence[float]) -> np.ndarray:
    if not energies_kcal:
        return np.array([], dtype=float)
    rel = np.array(energies_kcal, dtype=float)
    rel -= np.min(rel)
    exponents = -rel / (GAS_R_KCAL * TEMP_K)
    exponents -= np.max(exponents)
    w = np.exp(exponents)
    w /= np.sum(w)
    return w


def select_conformers(
    conformers: Sequence[ConformerResult],
    energy_window_kcal: float,
    boltzmann_cutoff: float,
) -> Tuple[List[ConformerResult], np.ndarray, List[float]]:
    if not conformers:
        return [], np.array([], dtype=float), []

    energies_kcal = [(c.energy_eh - conformers[0].energy_eh) * KCAL_PER_EH for c in conformers]

    filtered: List[ConformerResult] = []
    filtered_rel: List[float] = []
    for conf, rel in zip(conformers, energies_kcal):
        if rel <= energy_window_kcal + 1e-9:
            filtered.append(conf)
            filtered_rel.append(rel)

    if not filtered:
        filtered = [conformers[0]]
        filtered_rel = [0.0]

    weights = boltzmann_weights(filtered_rel)

    idx_sorted = np.argsort(-weights)
    keep: List[int] = []
    cumulative = 0.0
    for idx in idx_sorted:
        keep.append(int(idx))
        cumulative += float(weights[idx])
        if cumulative >= boltzmann_cutoff:
            break

    keep_sorted = sorted(keep)
    selected = [filtered[i] for i in keep_sorted]
    selected_rel = [filtered_rel[i] for i in keep_sorted]
    selected_weights = boltzmann_weights(selected_rel)
    return selected, selected_weights, selected_rel


def _safe_charge(atom: Chem.Atom) -> float:
    val = atom.GetProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else "0.0"
    parsed = parse_float(val)
    return parsed if parsed is not None else 0.0


def estimate_h_shift(mol: Chem.Mol, conf: Chem.Conformer, h_idx: int) -> float:
    h_atom = mol.GetAtomWithIdx(h_idx)
    neighbors = list(h_atom.GetNeighbors())
    if not neighbors:
        return 1.0
    heavy = neighbors[0]

    z = heavy.GetAtomicNum()
    if z == 6:
        if heavy.GetIsAromatic():
            base = 7.15
        elif heavy.GetHybridization() == Chem.HybridizationType.SP2:
            base = 5.35
        else:
            base = 1.15
    elif z == 7:
        base = 3.0
    elif z == 8:
        base = 2.3
    else:
        base = 1.8

    alpha_weights = {7: 0.45, 8: 0.65, 9: 1.1, 15: 0.3, 16: 0.35, 17: 0.9, 35: 0.7}
    beta_weights = {7: 0.14, 8: 0.18, 9: 0.25, 15: 0.10, 16: 0.12, 17: 0.20, 35: 0.15}

    for nbr in heavy.GetNeighbors():
        if nbr.GetIdx() == h_idx:
            continue
        base += alpha_weights.get(nbr.GetAtomicNum(), 0.0)
        for nbr2 in nbr.GetNeighbors():
            if nbr2.GetIdx() in {heavy.GetIdx(), h_idx}:
                continue
            base += beta_weights.get(nbr2.GetAtomicNum(), 0.0)

    base += max(-0.6, min(0.6, -0.65 * _safe_charge(heavy)))

    h_pos = conf.GetAtomPosition(h_idx)
    hetero_bonus = 0.0
    for atom in mol.GetAtoms():
        z2 = atom.GetAtomicNum()
        if z2 in {7, 8, 9, 15, 16, 17, 35}:
            pos = conf.GetAtomPosition(atom.GetIdx())
            dist = math.dist((h_pos.x, h_pos.y, h_pos.z), (pos.x, pos.y, pos.z))
            if dist < 3.0:
                hetero_bonus += (3.0 - dist) * 0.12

    shift = base + hetero_bonus
    return max(0.0, min(12.0, shift))


def group_atoms_by_rank(mol: Chem.Mol, atomic_number: int) -> List[List[int]]:
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    rank_to_atoms: Dict[int, List[int]] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != atomic_number:
            continue
        rank = int(ranks[atom.GetIdx()])
        rank_to_atoms.setdefault(rank, []).append(atom.GetIdx())

    groups = sorted(rank_to_atoms.values(), key=lambda atom_list: min(atom_list))
    return groups


def group_hydrogens(mol: Chem.Mol) -> List[List[int]]:
    return group_atoms_by_rank(mol, 1)


def parent_heavy_atom(mol: Chem.Mol, h_idx: int) -> Optional[int]:
    h_atom = mol.GetAtomWithIdx(h_idx)
    neigh = list(h_atom.GetNeighbors())
    return neigh[0].GetIdx() if neigh else None


def estimate_j_between_groups(mol: Chem.Mol, conf: Chem.Conformer, group_a: List[int], group_b: List[int]) -> float:
    h_a = group_a[0]
    h_b = group_b[0]

    p_a = parent_heavy_atom(mol, h_a)
    p_b = parent_heavy_atom(mol, h_b)
    if p_a is None or p_b is None:
        return 0.0

    if p_a == p_b:
        return 12.0

    atom_a = mol.GetAtomWithIdx(p_a)
    if any(n.GetIdx() == p_b for n in atom_a.GetNeighbors()):
        try:
            phi = rdMolTransforms.GetDihedralDeg(conf, h_a, p_a, p_b, h_b)
            r = math.radians(phi)
            j = 7.5 * (math.cos(r) ** 2) - 1.0 * math.cos(r) + 1.5
            return max(0.5, min(14.0, j))
        except Exception:
            return 7.0

    return 0.0


def estimate_c_shift(mol: Chem.Mol, conf: Chem.Conformer, c_idx: int) -> float:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6:
        return 0.0

    if atom.GetIsAromatic():
        base = 128.0
    elif atom.GetHybridization() == Chem.HybridizationType.SP2:
        base = 120.0
    else:
        base = 25.0

    has_o_double = False
    has_n_double = False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        z = other.GetAtomicNum()
        bo = float(bond.GetBondTypeAsDouble())
        if z == 8 and bo >= 1.9:
            has_o_double = True
        if z == 7 and bo >= 1.9:
            has_n_double = True
    if has_o_double:
        base = 175.0
    elif has_n_double:
        base = 155.0

    alpha_weights = {7: 13.0, 8: 20.0, 9: 36.0, 15: 6.0, 16: 8.0, 17: 12.0, 35: 10.0}
    beta_weights = {7: 2.5, 8: 4.0, 9: 8.0, 15: 1.5, 16: 2.0, 17: 3.0, 35: 2.0}

    for nbr in atom.GetNeighbors():
        base += alpha_weights.get(nbr.GetAtomicNum(), 0.0)
        for nbr2 in nbr.GetNeighbors():
            if nbr2.GetIdx() == atom.GetIdx():
                continue
            base += beta_weights.get(nbr2.GetAtomicNum(), 0.0)

    base += max(-12.0, min(12.0, -22.0 * _safe_charge(atom)))
    return max(0.0, min(230.0, base))


def estimate_f_shift(mol: Chem.Mol, conf: Chem.Conformer, f_idx: int) -> float:
    atom = mol.GetAtomWithIdx(f_idx)
    if atom.GetAtomicNum() != 9:
        return 0.0

    nbrs = list(atom.GetNeighbors())
    if not nbrs:
        return -120.0

    heavy = nbrs[0]
    if heavy.GetIsAromatic():
        base = -112.0
    elif heavy.GetAtomicNum() == 6 and heavy.GetHybridization() == Chem.HybridizationType.SP2:
        base = -98.0
    else:
        base = -128.0

    alpha_weights = {7: 10.0, 8: 14.0, 9: 17.0, 16: 6.0, 17: 5.0, 35: 4.0}
    beta_weights = {7: 2.0, 8: 3.0, 9: 4.0, 16: 1.5, 17: 1.5, 35: 1.0}

    for nbr in heavy.GetNeighbors():
        if nbr.GetIdx() == atom.GetIdx():
            continue
        base += alpha_weights.get(nbr.GetAtomicNum(), 0.0)
        for nbr2 in nbr.GetNeighbors():
            if nbr2.GetIdx() in {heavy.GetIdx(), atom.GetIdx()}:
                continue
            base += beta_weights.get(nbr2.GetAtomicNum(), 0.0)

    base += max(-25.0, min(25.0, -35.0 * _safe_charge(heavy)))
    return max(-260.0, min(120.0, base))


def estimate_p_shift(mol: Chem.Mol, conf: Chem.Conformer, p_idx: int) -> float:
    atom = mol.GetAtomWithIdx(p_idx)
    if atom.GetAtomicNum() != 15:
        return 0.0

    neighbors = list(atom.GetNeighbors())
    if not neighbors:
        return 0.0

    carbon_neighbors = 0
    oxygen_neighbors = 0
    hetero_neighbors = 0
    double_o = 0
    for nbr in neighbors:
        z = nbr.GetAtomicNum()
        if z == 6:
            carbon_neighbors += 1
        if z == 8:
            oxygen_neighbors += 1
        if z not in {1, 6}:
            hetero_neighbors += 1

        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
        order = float(bond.GetBondTypeAsDouble()) if bond is not None else 1.0
        if z == 8 and order >= 1.9:
            double_o += 1

    if double_o > 0:
        base = 8.0
        base += 14.0 * min(2, double_o)
        if oxygen_neighbors >= 2:
            base -= 14.0
        if carbon_neighbors >= 2:
            base += 10.0
    else:
        if carbon_neighbors >= 3:
            base = -24.0
        elif oxygen_neighbors >= 3:
            base = -4.0
        else:
            base = 10.0

    alpha_weights = {7: 8.0, 8: 10.0, 9: 18.0, 15: 6.0, 16: 10.0, 17: 8.0, 35: 6.0}
    beta_weights = {7: 2.0, 8: 2.5, 9: 4.0, 15: 1.5, 16: 2.0, 17: 1.5, 35: 1.0}

    for nbr in neighbors:
        base += alpha_weights.get(nbr.GetAtomicNum(), 0.0)
        for nbr2 in nbr.GetNeighbors():
            if nbr2.GetIdx() == atom.GetIdx():
                continue
            base += beta_weights.get(nbr2.GetAtomicNum(), 0.0)

    base += 2.5 * hetero_neighbors
    base += max(-40.0, min(40.0, -80.0 * _safe_charge(atom)))
    return max(-220.0, min(220.0, base))


def build_h_group_predictions(
    mol: Chem.Mol,
    selected_confs: Sequence[ConformerResult],
    weights: np.ndarray,
    frequency_mhz: float,
) -> List[GroupPrediction]:
    if not selected_confs:
        return []

    work_mol = Chem.Mol(mol)
    try:
        AllChem.ComputeGasteigerCharges(work_mol)
    except Exception:
        pass

    h_groups = group_hydrogens(work_mol)
    conf_map = {c.conf_id: work_mol.GetConformer(c.conf_id) for c in selected_confs}

    weighted_shifts: Dict[int, float] = {}
    for g_id, atoms in enumerate(h_groups, start=1):
        group_shift = 0.0
        for conf_weight, conf_result in zip(weights, selected_confs):
            conf = conf_map[conf_result.conf_id]
            atom_shifts = [estimate_h_shift(work_mol, conf, h_idx) for h_idx in atoms]
            group_shift += float(conf_weight) * float(np.mean(atom_shifts))
        weighted_shifts[g_id] = group_shift

    couplings: Dict[Tuple[int, int], float] = {}
    for i, atoms_i in enumerate(h_groups, start=1):
        for j, atoms_j in enumerate(h_groups, start=1):
            if j <= i:
                continue
            weighted_j = 0.0
            for conf_weight, conf_result in zip(weights, selected_confs):
                conf = conf_map[conf_result.conf_id]
                weighted_j += float(conf_weight) * estimate_j_between_groups(work_mol, conf, atoms_i, atoms_j)
            if weighted_j >= 0.8:
                couplings[(i, j)] = weighted_j

    predictions: List[GroupPrediction] = []
    for g_id, atoms in enumerate(h_groups, start=1):
        partners: List[Tuple[int, float, int]] = []
        for (a, b), j_hz in couplings.items():
            if a == g_id:
                partners.append((b, j_hz, len(h_groups[b - 1])))
            elif b == g_id:
                partners.append((a, j_hz, len(h_groups[a - 1])))

        if not partners:
            multiplicity = "singlet"
            j_hz = 0.0
        elif len(partners) == 1:
            n = partners[0][2]
            multiplicity = MULTIPLICITY_LABELS.get(n + 1, "multiplet")
            j_hz = partners[0][1]
        else:
            multiplicity = "multiplet"
            j_hz = float(np.mean([p[1] for p in partners]))

        shift_ppm = weighted_shifts[g_id]
        predictions.append(
            GroupPrediction(
                group_id=g_id,
                atom_indices=[idx + 1 for idx in atoms],
                shift_ppm=shift_ppm,
                shift_hz=shift_ppm * frequency_mhz,
                multiplicity=multiplicity,
                j_hz=j_hz,
                integral=float(len(atoms)),
            )
        )

    predictions.sort(key=lambda item: item.shift_ppm, reverse=True)
    return predictions


def build_nonproton_group_predictions(
    mol: Chem.Mol,
    selected_confs: Sequence[ConformerResult],
    weights: np.ndarray,
    frequency_mhz: float,
    nucleus: str,
) -> List[GroupPrediction]:
    if not selected_confs:
        return []

    if nucleus == "13c":
        atomic_number = 6
        estimator = estimate_c_shift
    elif nucleus == "19f":
        atomic_number = 9
        estimator = estimate_f_shift
    elif nucleus == "31p":
        atomic_number = 15
        estimator = estimate_p_shift
    else:
        return []

    work_mol = Chem.Mol(mol)
    try:
        AllChem.ComputeGasteigerCharges(work_mol)
    except Exception:
        pass

    groups = group_atoms_by_rank(work_mol, atomic_number)
    if not groups:
        return []

    conf_map = {c.conf_id: work_mol.GetConformer(c.conf_id) for c in selected_confs}
    predictions: List[GroupPrediction] = []
    for g_id, atoms in enumerate(groups, start=1):
        weighted_shift = 0.0
        for conf_weight, conf_result in zip(weights, selected_confs):
            conf = conf_map[conf_result.conf_id]
            atom_shifts = [estimator(work_mol, conf, atom_idx) for atom_idx in atoms]
            weighted_shift += float(conf_weight) * float(np.mean(atom_shifts))
        predictions.append(
            GroupPrediction(
                group_id=g_id,
                atom_indices=[idx + 1 for idx in atoms],
                shift_ppm=weighted_shift,
                shift_hz=weighted_shift * frequency_mhz,
                multiplicity="singlet",
                j_hz=0.0,
                integral=float(len(atoms)),
            )
        )

    predictions.sort(key=lambda item: item.shift_ppm, reverse=True)
    return predictions


def build_group_predictions(
    mol: Chem.Mol,
    selected_confs: Sequence[ConformerResult],
    weights: np.ndarray,
    frequency_mhz: float,
    nucleus: str,
) -> List[GroupPrediction]:
    if nucleus == "1h":
        return build_h_group_predictions(mol, selected_confs, weights, frequency_mhz)
    return build_nonproton_group_predictions(mol, selected_confs, weights, frequency_mhz, nucleus)


def splitting_pattern(multiplicity: str) -> Tuple[np.ndarray, np.ndarray]:
    if multiplicity == "singlet":
        return np.array([0.0]), np.array([1.0])
    if multiplicity == "doublet":
        return np.array([-0.5, 0.5]), np.array([1.0, 1.0])
    if multiplicity == "triplet":
        return np.array([-1.0, 0.0, 1.0]), np.array([1.0, 2.0, 1.0])
    if multiplicity == "quartet":
        return np.array([-1.5, -0.5, 0.5, 1.5]), np.array([1.0, 3.0, 3.0, 1.0])
    if multiplicity == "quintet":
        return np.array([-2.0, -1.0, 0.0, 1.0, 2.0]), np.array([1.0, 4.0, 6.0, 4.0, 1.0])
    return np.array([-1.2, -0.6, 0.0, 0.6, 1.2]), np.array([0.8, 1.0, 1.2, 1.0, 0.8])


def simulate_spectrum(
    groups: Sequence[GroupPrediction],
    frequency_mhz: float,
    line_shape: str,
    fwhm_hz: float,
    nucleus: str,
) -> List[Tuple[float, float]]:
    if not groups:
        return []

    if nucleus == "13c":
        margin = 20.0
    elif nucleus == "19f":
        margin = 35.0
    elif nucleus == "31p":
        margin = 45.0
    else:
        margin = 1.2

    min_ppm = min(g.shift_ppm for g in groups) - margin
    max_ppm = max(g.shift_ppm for g in groups) + margin
    n_points = 6000

    ppm_axis = np.linspace(max_ppm, min_ppm, n_points)
    intensity = np.zeros_like(ppm_axis)

    gamma = max(1e-5, (fwhm_hz / frequency_mhz) / 2.0)
    sigma = max(1e-5, (fwhm_hz / frequency_mhz) / 2.355)

    for group in groups:
        offsets, amps = splitting_pattern(group.multiplicity)
        amps = amps / np.sum(amps)
        for offset, amp in zip(offsets, amps):
            center = group.shift_ppm + (offset * group.j_hz / frequency_mhz)
            scale = amp * group.integral
            dx = ppm_axis - center
            if line_shape == "gaussian":
                line = np.exp(-((dx**2) / (2.0 * sigma * sigma)))
            elif line_shape == "voigt":
                g = np.exp(-((dx**2) / (2.0 * sigma * sigma)))
                l = (gamma * gamma) / (dx * dx + gamma * gamma)
                line = 0.5 * (g + l)
            else:
                line = (gamma * gamma) / (dx * dx + gamma * gamma)
            intensity += scale * line

    max_intensity = float(np.max(intensity)) if len(intensity) else 1.0
    if max_intensity <= 0:
        max_intensity = 1.0
    intensity /= max_intensity

    return list(zip(ppm_axis.tolist(), intensity.tolist()))


def _chirality_sign(mol: Chem.Mol) -> float:
    try:
        centers = Chem.FindMolChiralCenters(Chem.RemoveHs(mol), includeUnassigned=True)
    except Exception:
        centers = []
    if not centers:
        return 1.0
    r_count = sum(1 for _, tag in centers if tag == "R")
    s_count = sum(1 for _, tag in centers if tag == "S")
    if r_count == s_count:
        return 1.0
    return 1.0 if r_count > s_count else -1.0


def _estimate_cd_bands(mol: Chem.Mol, conf_count: int) -> List[Tuple[float, float, float]]:
    heavy = Chem.RemoveHs(mol)
    aromatic_rings = int(Chem.rdMolDescriptors.CalcNumAromaticRings(heavy))
    hetero_count = sum(1 for atom in heavy.GetAtoms() if atom.GetAtomicNum() not in {1, 6})
    atom_count = heavy.GetNumAtoms()
    chirality = _chirality_sign(heavy)

    bands: List[Tuple[float, float, float]] = []
    bands.append((190.0 + min(18.0, 0.9 * hetero_count), 9.0 + 0.25 * atom_count, 1.0 * chirality))
    bands.append((215.0 + min(24.0, 3.0 * aromatic_rings), 12.0 + 0.2 * atom_count, -0.75 * chirality))
    bands.append((245.0 + min(38.0, 2.0 * hetero_count + 4.0 * aromatic_rings), 15.0 + 0.2 * atom_count, 0.58 * chirality))
    if aromatic_rings > 0:
        bands.append((290.0 + min(80.0, 12.0 * aromatic_rings), 18.0 + 0.1 * atom_count, -0.35 * chirality))
    if conf_count > 12:
        bands.append((330.0, 24.0, 0.20 * chirality))
    return bands


def simulate_cd_spectrum(
    mol: Chem.Mol,
    selected_confs: Sequence[ConformerResult],
    weights: np.ndarray,
) -> List[Tuple[float, float]]:
    if not selected_confs:
        return []
    wavelength_axis = np.linspace(450.0, 180.0, 4200)
    intensity = np.zeros_like(wavelength_axis)
    bands = _estimate_cd_bands(mol, len(selected_confs))

    weight_scale = float(np.sum(weights)) if len(weights) > 0 else 1.0
    if weight_scale <= 0.0:
        weight_scale = 1.0

    for center, width, amplitude in bands:
        dx = wavelength_axis - center
        profile = np.exp(-0.5 * ((dx / max(1e-3, width)) ** 2))
        intensity += amplitude * profile * weight_scale

    max_abs = float(np.max(np.abs(intensity))) if len(intensity) else 1.0
    if max_abs <= 0:
        max_abs = 1.0
    intensity /= max_abs
    return list(zip(wavelength_axis.tolist(), intensity.tolist()))


def _find_cd_band_extrema(data: Sequence[Tuple[float, float]], limit: int = 8) -> List[Tuple[float, float]]:
    if len(data) < 5:
        return []
    x = np.array([p[0] for p in data], dtype=float)
    y = np.array([p[1] for p in data], dtype=float)
    extrema: List[Tuple[float, float]] = []
    for i in range(1, len(y) - 1):
        left = y[i - 1]
        mid = y[i]
        right = y[i + 1]
        is_max = mid > left and mid >= right
        is_min = mid < left and mid <= right
        if is_max or is_min:
            extrema.append((x[i], mid))
    extrema.sort(key=lambda item: abs(item[1]), reverse=True)
    return extrema[:limit]


# ---------------------------------------------------------------------------
# IR prediction
# ---------------------------------------------------------------------------

# Heuristic band table: (SMARTS, freq_cm1, rel_intensity_km_mol, label)
# Used as fallback when xTB Hessian is unavailable.
_IR_HEURISTIC_BANDS: List[Tuple[str, float, float, str]] = [
    # sp3 C-H stretches (present in nearly all organics)
    ("[CX4;H]",                                         2930.0, 60.0,  "C-H str (sp3)"),
    ("[CX4;H]",                                         2860.0, 40.0,  "C-H str (sp3 sym)"),
    # O-H stretches
    ("[OX2H;!$([OX2H][CX3]=O)]",                       3380.0, 100.0, "O-H str"),
    ("[CX3](=O)[OX2H]",                                 3000.0, 70.0,  "O-H str (COOH)"),
    # N-H stretches
    ("[NX3H2]",                                         3420.0, 55.0,  "N-H str (NH2)"),
    ("[NX3H2]",                                         3330.0, 55.0,  "N-H str (NH2 2nd)"),
    ("[NX3H1]",                                         3340.0, 50.0,  "N-H str"),
    # sp C-H
    ("[CX2H]#[CX2]",                                    3300.0, 35.0,  "C-H str (terminal alkyne)"),
    # sp2 C-H
    ("[cH]",                                            3020.0, 25.0,  "Ar C-H str"),
    ("[CX3H]=[CX3]",                                    3080.0, 25.0,  "=C-H str"),
    # Aldehyde C-H (Fermi doublet)
    ("[CX3H1](=O)",                                     2820.0, 35.0,  "C-H str (ald)"),
    ("[CX3H1](=O)",                                     2720.0, 30.0,  "C-H str (ald 2nd)"),
    # Carbonyl C=O stretches
    ("[CX3H1](=O)",                                     1725.0, 100.0, "C=O str (ald)"),
    ("[CX3](=O)[CX4]",                                  1715.0, 100.0, "C=O str (ket)"),
    ("[CX3](=O)[OX2H0][#6]",                            1735.0, 100.0, "C=O str (ester)"),
    ("[CX3](=O)[OX2H]",                                 1710.0, 100.0, "C=O str (acid)"),
    ("[CX3](=O)[NX3]",                                  1660.0, 90.0,  "C=O str (amide)"),
    # Triple bonds
    ("[C]#[N]",                                         2230.0, 90.0,  "C≡N str"),
    ("[CX2]#[CX2]",                                     2150.0, 35.0,  "C≡C str"),
    # Aromatic C=C
    ("a:a",                                             1600.0, 40.0,  "C=C str (arom)"),
    ("a:a",                                             1500.0, 30.0,  "C=C str (arom 2nd)"),
    # Nitro
    ("[$([N+](=O)[O-]),$([N](=O)=O)]",                 1540.0, 90.0,  "NO2 asym str"),
    ("[$([N+](=O)[O-]),$([N](=O)=O)]",                 1350.0, 70.0,  "NO2 sym str"),
    # Ethers and esters C-O
    ("[OD2;X2]([CX4])[CX4]",                           1100.0, 70.0,  "C-O-C str (ether)"),
    ("[CX3](=O)[OX2H0][CX4]",                          1200.0, 60.0,  "C-O str (ester)"),
    # Halogens
    ("[F]",                                             1150.0, 80.0,  "C-F str"),
    ("[Cl]",                                             730.0, 60.0,  "C-Cl str"),
    ("[Br]",                                             590.0, 50.0,  "C-Br str"),
    # Sulfur
    ("[SX4](=O)(=O)",                                   1330.0, 90.0,  "SO2 asym str"),
    ("[SX4](=O)(=O)",                                   1150.0, 80.0,  "SO2 sym str"),
    ("[SX3](=O)",                                       1060.0, 80.0,  "S=O str"),
    # Phosphorus
    ("[PX4]=O",                                         1200.0, 90.0,  "P=O str"),
]


def _ir_heuristic_peaks(mol: Chem.Mol) -> List[Tuple[float, float, str]]:
    """Group-based heuristic IR peaks. Returns (freq_cm1, rel_intensity, label)."""
    peaks: List[Tuple[float, float, str]] = []
    seen: set = set()
    for smarts, freq, intensity, label in _IR_HEURISTIC_BANDS:
        key = (smarts, freq)
        if key in seen:
            continue
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if pattern is None:
            continue
        try:
            matches = mol.GetSubstructMatches(pattern)
        except Exception:
            continue
        if not matches:
            continue
        seen.add(key)
        # Scale slightly with match count but cap contribution
        scaled = intensity * min(1.5, 0.7 + 0.3 * math.sqrt(len(matches)))
        peaks.append((freq, scaled, label))
    return peaks


def _label_ir_band(freq_cm1: float) -> str:
    """Assign a rough region label to an xTB-derived vibrational frequency."""
    if freq_cm1 >= 3000:
        return "X-H str"
    if freq_cm1 >= 2500:
        return "C-H str"
    if freq_cm1 >= 2000:
        return "triple bond"
    if freq_cm1 >= 1700:
        return "C=O str"
    if freq_cm1 >= 1500:
        return "C=C/N str"
    if freq_cm1 >= 1200:
        return "C-O/N str"
    return "fingerprint"


def _parse_vibspectrum(path: Path) -> List[Tuple[float, float]]:
    """Parse xTB vibspectrum file. Returns (freq_cm1, I_IR_km_mol) for real modes (>50 cm⁻¹)."""
    peaks: List[Tuple[float, float]] = []
    if not path.exists():
        return peaks
    try:
        in_block = False
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("$vibrational spectrum"):
                    in_block = True
                    continue
                if line.startswith("$end"):
                    break
                if not in_block or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    freq = float(parts[2])
                    ir_int = float(parts[3])
                except (ValueError, IndexError):
                    continue
                if freq < 50.0:
                    continue
                peaks.append((freq, max(0.0, ir_int)))
    except Exception:
        pass
    return peaks


def _parse_xtb_ir_stdout(stdout: str) -> List[Tuple[float, float]]:
    """Extract frequencies + IR intensities from xTB stdout (fallback to vibspectrum file)."""
    peaks: List[Tuple[float, float]] = []
    freq_values: List[float] = []
    intens_values: List[float] = []
    freq_re = re.compile(r"eigval\s*:([\s\d.+\-eE]+)")
    intens_re = re.compile(r"intens\s*:([\s\d.+\-eE]+)")
    for line in stdout.splitlines():
        m = freq_re.search(line)
        if m:
            for tok in m.group(1).split():
                try:
                    freq_values.append(float(tok))
                except ValueError:
                    pass
        m = intens_re.search(line)
        if m:
            for tok in m.group(1).split():
                try:
                    intens_values.append(float(tok))
                except ValueError:
                    pass
    for freq, intens in zip(freq_values, intens_values):
        if freq >= 50.0:
            peaks.append((freq, abs(intens)))
    return peaks


def run_xtb_hessian(
    xtb_bin: str,
    mol: Chem.Mol,
    conf_id: int,
    solvent: str,
    output_dir: Path,
    timeout_s: int,
) -> List[Tuple[float, float]]:
    """Run xTB GFN2 Hessian on one conformer. Returns (freq_cm1, I_IR_km_mol) for real modes."""
    archive_dir = output_dir / f"xtb_hess_conf_{conf_id}"
    archive_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=f"easynmr_xtb_hess_{conf_id}_"))
    log_path = archive_dir / "xtb_hess.log"
    stdout = ""
    stderr = ""
    try:
        xyz_path = run_dir / "input.xyz"
        write_conf_xyz(mol, conf_id, xyz_path, f"conf-{conf_id}")

        cmd = [xtb_bin, str(xyz_path), "--gfn2", "--hess"]
        if solvent:
            cmd.extend(["--alpb", solvent])

        try:
            proc = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
        except Exception:
            return []

        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                proc.kill()
            stdout, stderr = proc.communicate()
            log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
            return []

        log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
        if proc.returncode != 0:
            return []

        vibspec_path = run_dir / "vibspectrum"
        peaks = _parse_vibspectrum(vibspec_path)
        if peaks:
            try:
                shutil.copy2(vibspec_path, archive_dir / "vibspectrum")
            except Exception:
                pass
            return peaks
        return _parse_xtb_ir_stdout(stdout)
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def simulate_ir_spectrum(
    sticks: Sequence[Tuple[float, float]],
    scale_factor: float = IR_SCALE_FACTOR_GFN2,
    fwhm_cm1: float = 8.0,
    n_points: int = 2000,
) -> List[Tuple[float, float]]:
    """Broaden an IR stick spectrum into a continuous trace (Lorentzian).

    Returns list of (wavenumber_cm1, intensity) covering 4000→400 cm⁻¹,
    normalised to a max intensity of 1.0.
    """
    wn_axis = np.linspace(4000.0, 400.0, n_points)
    intensity = np.zeros(n_points)
    half_gamma = fwhm_cm1 / 2.0
    for freq, intens in sticks:
        scaled = freq * scale_factor
        denom = (wn_axis - scaled) ** 2 + half_gamma ** 2
        intensity += intens * half_gamma / (np.pi * np.maximum(denom, 1e-12))
    max_int = float(np.max(intensity))
    if max_int > 0.0:
        intensity /= max_int
    return list(zip(wn_axis.tolist(), intensity.tolist()))


def predict_ir_spectrum(
    mol: Chem.Mol,
    selected_confs: Sequence[ConformerResult],
    weights: np.ndarray,
    solvent: str,
    output_dir: Path,
    warnings: List[str],
    progress_path: Optional[Path] = None,
) -> Tuple[List[Tuple[float, float]], List[Tuple[float, float, str]], bool]:
    """Predict an IR spectrum.

    Runs xTB Hessian on the top Boltzmann-weighted conformers (max IR_MAX_HESS_CONFS),
    averages the stick spectra by Boltzmann weight, then broadens with Lorentzian.
    Falls back to group-based heuristics if xTB is unavailable or fails.

    Returns:
        spectrum: list of (wavenumber_cm1, intensity) – broadened trace
        labeled_peaks: list of (wavenumber_cm1, intensity, label) – key peaks for display
        used_heuristic: True if the heuristic fallback was used
    """
    xtb_bin_env = os.environ.get("EASYSPECTRA_XTB", "xtb")
    xtb_bin = None if xtb_bin_env == "__none__" else shutil.which(xtb_bin_env)
    timeout_s = int(os.environ.get("EASYSPECTRA_XTB_TIMEOUT", "25"))

    used_heuristic = False
    raw_sticks: List[Tuple[float, float]] = []

    if xtb_bin and selected_confs:
        n = min(len(selected_confs), IR_MAX_HESS_CONFS)
        hess_confs = selected_confs[:n]
        hess_weights = weights[:n].copy()
        w_sum = float(np.sum(hess_weights))
        if w_sum <= 0.0:
            w_sum = 1.0
        hess_weights /= w_sum

        valid = 0
        for i, conf in enumerate(hess_confs):
            if progress_path is not None:
                write_progress(
                    progress_path, "ir_hessian",
                    f"Running xTB Hessian on conformer {i + 1}/{n}",
                    0.84 + 0.06 * (i / max(1, n)),
                )
            conf_peaks = run_xtb_hessian(xtb_bin, mol, conf.conf_id, solvent, output_dir, timeout_s)
            if conf_peaks:
                w = float(hess_weights[i])
                for freq, intens in conf_peaks:
                    raw_sticks.append((freq, intens * w))
                valid += 1

        if valid == 0:
            warnings.append(
                "xTB Hessian failed for all selected conformers; "
                "falling back to group-based IR estimate."
            )
            used_heuristic = True
    else:
        warnings.append(
            "xTB not available for IR Hessian prediction; "
            "using group-based IR estimate (qualitative only)."
        )
        used_heuristic = True

    # Build labeled peaks list
    labeled_peaks: List[Tuple[float, float, str]] = []
    if used_heuristic:
        heuristic = _ir_heuristic_peaks(mol)
        raw_sticks = [(f, i) for f, i, _ in heuristic]
        labeled_peaks = heuristic
        fwhm = 20.0
    else:
        # Label xTB peaks by frequency region
        labeled_peaks = [(f, i, _label_ir_band(f)) for f, i in raw_sticks]
        fwhm = 8.0

    spectrum = simulate_ir_spectrum(raw_sticks, IR_SCALE_FACTOR_GFN2, fwhm)
    return spectrum, labeled_peaks, used_heuristic


def write_ir_spectrum_csv(path: Path, data: Sequence[Tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber", "intensity"])
        for wn, intensity in data:
            writer.writerow([f"{wn:.2f}", f"{intensity:.8f}"])


def write_ir_peaks_csv(
    path: Path,
    labeled_peaks: Sequence[Tuple[float, float, str]],
    scale_factor: float = IR_SCALE_FACTOR_GFN2,
    top_n: int = 24,
) -> None:
    """Write IR peaks in a format the C++ peak browser can display.

    Columns mirror the NMR peaks CSV so the C++ read_peak_rows parser can load them:
    peak_id, wavenumber_cm1, wavenumber_cm1, label, j_hz, rel_intensity, assignment, atom_indices
    """
    sorted_peaks = sorted(labeled_peaks, key=lambda x: x[1], reverse=True)[:top_n]
    sorted_peaks.sort(key=lambda x: x[0], reverse=True)  # display high→low
    max_int = max((i for _, i, _ in sorted_peaks), default=1.0) or 1.0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak_id", "wavenumber_cm1", "wavenumber_cm1", "label", "j_hz", "rel_intensity", "assignment", "atom_indices"])
        for idx, (freq, intens, label) in enumerate(sorted_peaks, start=1):
            scaled = freq * scale_factor
            rel = intens / max_int
            writer.writerow([idx, f"{scaled:.1f}", f"{scaled:.1f}", label, "0.00", f"{rel:.4f}", "", ""])


# ── UV-Vis Prediction ─────────────────────────────────────────────────────────

def _parse_xtb_stda_stdout(stdout: str) -> List[Tuple[float, float]]:
    """Extract excited-state energies and oscillator strengths from xTB --stda stdout.

    Returns list of (wavelength_nm, oscillator_strength) for allowed transitions
    (f > 1e-5) with wavelength in 150–900 nm range.
    """
    sticks: List[Tuple[float, float]] = []

    # Primary pattern: tabular block with state, Erel/eV, lambda/nm, f columns
    # xTB sTDA prints:
    #   state    Erel/eV    lambda/nm      f        ...
    #       1       2.4543     505.49     0.1234   ...
    table_re = re.compile(
        r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+(?:[eE][+-]?\d+)?)",
        re.MULTILINE,
    )
    in_block = False
    for line in stdout.splitlines():
        if re.search(r"lambda.*nm|wavelength.*nm|excited\s+state|EXCITED STATE", line, re.IGNORECASE):
            in_block = True
            continue
        if in_block:
            m = table_re.match(line)
            if m:
                try:
                    # columns: state, Erel_eV, lambda_nm, f
                    lambda_nm = float(m.group(3))
                    f_osc = float(m.group(4))
                    if 150.0 <= lambda_nm <= 900.0 and f_osc > 1e-5:
                        sticks.append((lambda_nm, f_osc))
                except ValueError:
                    pass
            elif sticks:
                # Stop at blank/separator line after at least one state found
                stripped = line.strip()
                if not stripped or stripped.startswith("---") or stripped.startswith("==="):
                    in_block = False

    if sticks:
        return sticks

    # Fallback: scan for pattern "N  eV_value eV  (nm_value nm)"
    ev_nm_re = re.compile(r"([\d.]+)\s+eV\s+\(?\s*([\d.]+)\s+nm", re.IGNORECASE)
    for line in stdout.splitlines():
        m = ev_nm_re.search(line)
        if m:
            try:
                lambda_nm = float(m.group(2))
                # Try to find oscillator strength on same line
                f_m = re.search(r"f\s*=\s*([\d.eE+\-]+)", line, re.IGNORECASE)
                f_osc = float(f_m.group(1)) if f_m else 0.05
                if 150.0 <= lambda_nm <= 900.0 and f_osc > 1e-5:
                    sticks.append((lambda_nm, f_osc))
            except ValueError:
                pass

    if sticks:
        return sticks

    # Last resort: look for any line containing an energy in eV and convert
    eV_re = re.compile(r"^\s*\d+\s+([\d.]+)\s+eV", re.MULTILINE)
    for m in eV_re.finditer(stdout):
        try:
            ev = float(m.group(1))
            if 1.0 <= ev <= 10.0:
                lambda_nm = 1240.0 / ev
                if 150.0 <= lambda_nm <= 900.0:
                    sticks.append((lambda_nm, 0.05))
        except ValueError:
            pass

    return sticks


def run_xtb_stda(
    xtb_bin: str,
    mol: Chem.Mol,
    conf_id: int,
    solvent: str,
    output_dir: Path,
    timeout_s: int,
) -> List[Tuple[float, float]]:
    """Run xTB GFN2 --stda on one conformer.

    Returns list of (wavelength_nm, oscillator_strength) for allowed transitions.
    Returns [] on failure or timeout.
    """
    archive_dir = output_dir / f"xtb_stda_conf_{conf_id}"
    archive_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=f"easynmr_xtb_stda_{conf_id}_"))
    log_path = archive_dir / "xtb_stda.log"
    stdout = ""
    stderr = ""
    try:
        xyz_path = run_dir / "input.xyz"
        # Prefer previously optimised geometry if available
        archived_opt = output_dir.parent / f"xtb_conf_{conf_id}" / "xtbopt.xyz"
        if archived_opt.exists():
            shutil.copy2(archived_opt, xyz_path)
        else:
            write_conf_xyz(mol, conf_id, xyz_path, f"conf-{conf_id}")

        cmd = [xtb_bin, str(xyz_path), "--gfn2", "--stda"]
        if solvent:
            cmd.extend(["--alpb", solvent])

        try:
            proc = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
        except Exception:
            return []

        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                proc.kill()
            stdout, stderr = proc.communicate()
            log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
            return []

        log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
        if proc.returncode != 0:
            return []

        return _parse_xtb_stda_stdout(stdout)
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def _uvvis_heuristic_bands(mol: Chem.Mol) -> List[Tuple[float, float, str]]:
    """Return heuristic UV-Vis bands from SMARTS chromophore matching.

    Returns list of (wavelength_nm, relative_intensity, label).
    More specific patterns are checked first; the most specific match wins for
    overlapping chromophores (de-duplicated by match atom sets).
    """
    # (smarts, wavelength_nm, rel_intensity, label)
    # More specific patterns listed first so they shadow generic ones
    table: List[Tuple[str, float, float, str]] = [
        # Extended / fused aromatic
        ("c1ccc2c(c1)ccc3ccccc23",      340.0, 0.90, "\u03c0\u2192\u03c0* (anthracene)"),
        ("c1ccc2ccccc2c1",              315.0, 0.85, "\u03c0\u2192\u03c0* (naph)"),
        # Enone / dienone
        ("[CX3](=O)[CX3]=[CX3][CX3]=[CX3]", 350.0, 0.80, "\u03c0\u2192\u03c0* (cross-conj)"),
        ("[CX3](=O)[CX3]=[CX3]",        328.0, 0.75, "\u03c0\u2192\u03c0* (enone)"),
        # Nitro
        ("[NX3](=O)=O",                 330.0, 0.60, "n\u2192\u03c0* (NO\u2082)"),
        # Thione
        ("[CX3](=S)",                   320.0, 0.40, "n\u2192\u03c0* (C=S)"),
        # Diene
        ("[CX3]=[CX3]-[CX3]=[CX3]",    217.0, 0.80, "\u03c0\u2192\u03c0* (diene)"),
        # Aldehyde / ketone n→π*
        ("[CX3H1](=O)",                 293.0, 0.30, "n\u2192\u03c0* (CHO)"),
        ("[CX3](=O)[#6]",               280.0, 0.20, "n\u2192\u03c0* (C=O)"),
        # Substituted benzene (check before unsubstituted)
        ("c1ccc(-[!#1])cc1",            272.0, 0.65, "\u03c0\u2192\u03c0* (Ar-sub)"),
        # Unsubstituted benzene
        ("c1ccccc1",                    254.0, 0.50, "\u03c0\u2192\u03c0* (ArH)"),
        # Carboxylic acid / ester / amide
        ("[CX3](=O)[OX2H1]",            210.0, 0.40, "n\u2192\u03c0* (COOH)"),
        ("[CX3](=O)[OX2][#6]",          210.0, 0.35, "n\u2192\u03c0* (ester)"),
        ("[CX3](=O)[NX3]",              220.0, 0.45, "n\u2192\u03c0* (amide)"),
        # Isolated alkene
        ("[CX3]=[CX3]",                 185.0, 0.30, "\u03c0\u2192\u03c0* (C=C)"),
    ]

    bands: List[Tuple[float, float, str]] = []
    seen_atoms: set = set()

    for smarts, wl, rel_int, label in table:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern):
            match_set = frozenset(match)
            # Skip if the majority of these atoms were already covered
            if len(match_set & seen_atoms) >= len(match_set) // 2:
                continue
            seen_atoms.update(match_set)
            bands.append((wl, rel_int, label))

    return bands


def simulate_uvvis_spectrum(
    sticks: Sequence[Tuple[float, float]],
    fwhm_nm: float = UVVIS_BROADENING_NM,
    n_points: int = 1500,
) -> List[Tuple[float, float]]:
    """Broaden a UV-Vis stick spectrum into a continuous trace (Gaussian).

    Returns list of (wavelength_nm, intensity) covering 200→800 nm,
    normalised to a max intensity of 1.0.
    """
    wl_axis = np.linspace(200.0, 800.0, n_points)
    intensity = np.zeros(n_points)
    sigma = fwhm_nm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    for wl, f_osc in sticks:
        diff = wl_axis - wl
        intensity += f_osc * np.exp(-0.5 * (diff / sigma) ** 2)
    max_int = float(np.max(intensity))
    if max_int > 0.0:
        intensity /= max_int
    return list(zip(wl_axis.tolist(), intensity.tolist()))


def predict_uvvis_spectrum(
    mol: Chem.Mol,
    selected_confs: Sequence,
    weights: Sequence[float],
    solvent: str,
    output_dir: Path,
    warnings: List[str],
    progress_path: Optional[Path],
) -> Tuple[List[Tuple[float, float]], List[Tuple[float, float, str]], bool]:
    """Predict UV-Vis absorption spectrum via xTB --stda with heuristic fallback.

    Returns:
        (broadened_spectrum, labeled_peaks, used_heuristic)
        labeled_peaks: list of (wavelength_nm, oscillator_strength, label)
    """
    xtb_bin_env = os.environ.get("EASYSPECTRA_XTB", "xtb")
    xtb_bin = None if xtb_bin_env == "__none__" else shutil.which(xtb_bin_env)
    timeout_s = int(os.environ.get("EASYSPECTRA_XTB_TIMEOUT", "25"))

    used_heuristic = False
    combined_sticks: List[Tuple[float, float]] = []

    if xtb_bin:
        n_confs = min(len(selected_confs), UVVIS_MAX_STDA_CONFS)
        confs_subset = selected_confs[:n_confs]
        weights_subset = list(weights[:n_confs])
        weight_sum = sum(weights_subset)
        if weight_sum > 0:
            weights_subset = [w / weight_sum for w in weights_subset]

        for conf_result, w in zip(confs_subset, weights_subset):
            if progress_path:
                write_progress(progress_path, "uvvis_stda",
                               f"xTB sTDA: conformer {conf_result.conf_id}", 0.88)
            sticks = run_xtb_stda(xtb_bin, mol, conf_result.conf_id, solvent, output_dir, timeout_s)
            for wl, f_osc in sticks:
                combined_sticks.append((wl, f_osc * w))

    if not combined_sticks:
        if xtb_bin:
            warnings.append(
                "UV-Vis: xTB sTDA produced no excitations — using chromophore heuristic."
            )
        else:
            warnings.append(
                "UV-Vis: xTB not available — using chromophore heuristic (approximate)."
            )
        raw_bands = _uvvis_heuristic_bands(mol)
        if not raw_bands:
            # Featureless molecule — add a generic far-UV band
            raw_bands = [(190.0, 0.3, "\u03c3\u2192\u03c3* (generic)")]
        combined_sticks = [(wl, rel_int) for wl, rel_int, _ in raw_bands]
        used_heuristic = True

    if progress_path:
        write_progress(progress_path, "uvvis_spectrum_simulation",
                       "Simulating UV-Vis spectrum", 0.93)

    fwhm = UVVIS_BROADENING_NM if not used_heuristic else 25.0
    spectrum = simulate_uvvis_spectrum(combined_sticks, fwhm_nm=fwhm)

    # Build labeled peaks: top N unique bands by intensity
    # For xTB sticks, group nearby transitions (within 20 nm)
    peaks_by_wl: dict = {}
    for wl, f in combined_sticks:
        bucket = round(wl / 10.0) * 10
        if bucket not in peaks_by_wl or f > peaks_by_wl[bucket][1]:
            peaks_by_wl[bucket] = (wl, f, "")

    labeled_peaks: List[Tuple[float, float, str]] = []
    if used_heuristic:
        raw_bands_sorted = sorted(_uvvis_heuristic_bands(mol), key=lambda x: x[1], reverse=True)
        for wl, rel_int, label in raw_bands_sorted[:20]:
            labeled_peaks.append((wl, rel_int, label))
    else:
        for wl, f, label in sorted(peaks_by_wl.values(), key=lambda x: x[1], reverse=True)[:20]:
            labeled_peaks.append((wl, f, "S\u2099"))

    return spectrum, labeled_peaks, used_heuristic


def write_uvvis_spectrum_csv(path: Path, data: Sequence[Tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavelength", "intensity"])
        for wl, intensity in data:
            writer.writerow([f"{wl:.2f}", f"{intensity:.8f}"])


def write_uvvis_peaks_csv(
    path: Path,
    labeled_peaks: Sequence[Tuple[float, float, str]],
    top_n: int = 20,
) -> None:
    """Write UV-Vis peaks in the 8-column format the C++ peak browser expects."""
    sorted_peaks = sorted(labeled_peaks, key=lambda x: x[1], reverse=True)[:top_n]
    sorted_peaks.sort(key=lambda x: x[0])  # display short→long wavelength
    max_f = max((f for _, f, _ in sorted_peaks), default=1.0) or 1.0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak_id", "wavelength_nm", "wavelength_nm", "label",
                         "f_osc", "rel_intensity", "assignment", "atom_indices"])
        for idx, (wl, f_osc, label) in enumerate(sorted_peaks, start=1):
            rel = f_osc / max_f
            writer.writerow([idx, f"{wl:.1f}", f"{wl:.1f}", label,
                             f"{f_osc:.4f}", f"{rel:.4f}", "", ""])


# ── Molecular Properties ──────────────────────────────────────────────────────

# Vacuum-to-NHE shift (Trasatti convention, eV)
_VACUUM_TO_NHE_EV: float = 4.44
# Fc/Fc+ vs NHE (V), used for conversion to ferrocene reference
_FC_VS_NHE_V: float = 0.40


def compute_rdkit_descriptors(mol: Chem.Mol) -> dict:
    """Return a dict of common RDKit molecular descriptors."""
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    mol_no_h = Chem.RemoveHs(mol)
    mw = Descriptors.MolWt(mol_no_h)
    exact_mw = Descriptors.ExactMolWt(mol_no_h)
    formula = rdMolDescriptors.CalcMolFormula(mol_no_h)
    logp = Crippen.MolLogP(mol_no_h)
    tpsa = rdMolDescriptors.CalcTPSA(mol_no_h)
    hbd = rdMolDescriptors.CalcNumHBD(mol_no_h)
    hba = rdMolDescriptors.CalcNumHBA(mol_no_h)
    rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol_no_h)
    ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol_no_h)
    heavy_atom_count = mol_no_h.GetNumAtoms()
    formal_charge = Chem.GetFormalCharge(mol_no_h)
    return {
        "formula": formula,
        "mw": round(mw, 3),
        "exact_mw": round(exact_mw, 5),
        "logp": round(logp, 3),
        "tpsa": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "rotbonds": rotbonds,
        "ar_rings": ar_rings,
        "heavy_atom_count": heavy_atom_count,
        "formal_charge": formal_charge,
    }


def _parse_xtb_cm5_charges(stdout: str, n_atoms: int) -> List[float]:
    """Extract per-atom CM5 charges from xTB stdout."""
    charges: List[float] = []
    # xTB prints a block like:
    #  #   Z          covCN         q      C6AA      α(0)
    #  1   6 C        3.996    -0.069   ...
    # OR a "Mulliken/CM5 charges" block. Try the tabular block first.
    # Pattern: lines with   <idx>   <Z>   <sym>  ...  <charge>  ...
    block_re = re.compile(
        r"^\s*(\d+)\s+\d+\s+\w+\s+[\d.]+\s+([-+]?\d+\.\d+)", re.MULTILINE
    )
    # Try the "#   Z  covCN  q" table (charge is 4th numeric column)
    header_found = False
    in_table = False
    table_charges: List[float] = []
    for line in stdout.splitlines():
        if re.search(r"#\s+Z\s+covCN\s+q", line):
            header_found = True
            in_table = True
            table_charges = []
            continue
        if in_table:
            stripped = line.strip()
            if not stripped or stripped.startswith("---"):
                in_table = False
                continue
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    table_charges.append(float(parts[3]))
                except ValueError:
                    in_table = False
    if header_found and len(table_charges) == n_atoms:
        return table_charges

    # Fallback: "Mulliken/CM5 charges" block
    in_cm5 = False
    cm5_charges: List[float] = []
    for line in stdout.splitlines():
        if re.search(r"Mulliken/CM5 charges", line, re.IGNORECASE):
            in_cm5 = True
            cm5_charges = []
            continue
        if in_cm5:
            stripped = line.strip()
            if not stripped or stripped.startswith("---"):
                in_cm5 = False
                continue
            parts = stripped.split()
            # format:  idx  sym  mulliken  cm5
            if len(parts) >= 4:
                try:
                    cm5_charges.append(float(parts[3]))
                except ValueError:
                    in_cm5 = False
    if len(cm5_charges) == n_atoms:
        return cm5_charges

    # Last resort: generic regex scan
    matches = block_re.findall(stdout)
    if len(matches) == n_atoms:
        return [float(q) for _, q in matches]
    return []


def _parse_xtb_homo_lumo(stdout: str) -> Tuple[float, float, float]:
    """Return (homo_ev, lumo_ev, gap_ev) from xTB stdout. Returns (0, 0, 0) on failure."""
    gap_m = re.search(r"HOMO-LUMO gap\s*[:\s]+([\d.]+)\s*eV", stdout, re.IGNORECASE)
    # xTB orbital table format:  <idx>  <occ>  <energy_Eh>  <energy_eV> (HOMO)
    homo_m = re.search(r"([-\d.]+)\s+\(HOMO\)", stdout)
    lumo_m = re.search(r"([-\d.]+)\s+\(LUMO\)", stdout)
    # Fallback patterns
    if not homo_m:
        homo_m = re.search(r"HOMO\s+([-\d.]+)\s+eV", stdout)
    if not lumo_m:
        lumo_m = re.search(r"LUMO\s+([-\d.]+)\s+eV", stdout)
    try:
        homo_ev = float(homo_m.group(1)) if homo_m else 0.0
        lumo_ev = float(lumo_m.group(1)) if lumo_m else 0.0
        gap_ev = float(gap_m.group(1)) if gap_m else (lumo_ev - homo_ev if homo_ev and lumo_ev else 0.0)
        return homo_ev, lumo_ev, gap_ev
    except (AttributeError, ValueError):
        return 0.0, 0.0, 0.0


def _parse_xtb_dipole(stdout: str) -> float:
    """Return dipole magnitude in Debye from xTB stdout."""
    # xTB prints: "molecular dipole:" block with "full:" line
    m = re.search(r"full:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", stdout)
    if m:
        try:
            return float(m.group(4))  # last field is magnitude in Debye
        except ValueError:
            pass
    # Fallback: first "total" dipole line
    m = re.search(r"dipole moment\s*[:\s]+([\d.]+)\s+Debye", stdout, re.IGNORECASE)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return 0.0


def run_xtb_sp(
    xtb_bin: str,
    mol: Chem.Mol,
    conf_id: int,
    solvent: str,
    output_dir: Path,
    charge: int,
    timeout_s: int,
    label: str = "sp",
) -> Optional[dict]:
    """Run xTB GFN2 single-point with --pop to get CM5 charges + electronic properties.

    Returns dict with keys: charges, homo_ev, lumo_ev, gap_ev, dipole_debye.
    Returns None on failure.
    """
    archive_dir = output_dir / f"xtb_{label}_conf_{conf_id}"
    archive_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=f"easynmr_xtb_{label}_{conf_id}_"))
    log_path = archive_dir / "xtb_sp.log"
    stdout = ""
    stderr = ""

    try:
        # Prefer previously optimised geometry if available
        xyz_path = run_dir / "input.xyz"
        archived_opt = output_dir.parent / f"xtb_conf_{conf_id}" / "xtbopt.xyz"
        if archived_opt.exists():
            shutil.copy2(archived_opt, xyz_path)
        else:
            write_conf_xyz(mol, conf_id, xyz_path, f"conf-{conf_id}")

        # For open-shell calculations (odd total electrons), xTB needs --uhf 1.
        # This applies whenever the N+1 or N-1 charge makes an even-electron molecule
        # into an odd-electron (doublet) system.
        total_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - charge
        n_unpaired = total_electrons % 2  # 1 for doublet, 0 for closed-shell singlet
        cmd = [xtb_bin, str(xyz_path), "--gfn2", "--pop",
               "--chrg", str(charge), "--uhf", str(n_unpaired)]
        if solvent:
            cmd.extend(["--alpb", solvent])

        try:
            proc = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
        except Exception:
            return None

        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                proc.kill()
            stdout, stderr = proc.communicate()
            log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
            return None

        log_path.write_text((stdout or "") + "\n" + (stderr or ""), encoding="utf-8")
        if proc.returncode != 0:
            return None

        n_atoms = mol.GetNumAtoms()
        charges = _parse_xtb_cm5_charges(stdout, n_atoms)
        homo_ev, lumo_ev, gap_ev = _parse_xtb_homo_lumo(stdout)
        dipole_debye = _parse_xtb_dipole(stdout)

        if not charges:
            return None

        return {
            "charges": charges,
            "homo_ev": homo_ev,
            "lumo_ev": lumo_ev,
            "gap_ev": gap_ev,
            "dipole_debye": dipole_debye,
        }
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def compute_fukui_indices(
    charges_n: List[float],
    charges_np1: List[float],
    charges_nm1: List[float],
) -> List[dict]:
    """Compute condensed Fukui indices from CM5 atomic charges.

    f+(i)  = q_N(i) - q_{N+1}(i)   electrophilic attack site
    f-(i)  = q_{N-1}(i) - q_N(i)   nucleophilic attack site
    f0(i)  = (f+(i) + f-(i)) / 2    radical attack site
    Each index is normalized so the positive values sum to 1.0.
    """
    n = len(charges_n)
    f_plus = [charges_n[i] - charges_np1[i] for i in range(n)]
    f_minus = [charges_nm1[i] - charges_n[i] for i in range(n)]
    f_zero = [(f_plus[i] + f_minus[i]) / 2.0 for i in range(n)]

    def _normalize(values: List[float]) -> List[float]:
        clamped = [max(0.0, v) for v in values]
        total = sum(clamped)
        if total > 0.0:
            return [v / total for v in clamped]
        return clamped

    f_plus = _normalize(f_plus)
    f_minus = _normalize(f_minus)
    f_zero = _normalize(f_zero)
    return [
        {"f_plus": round(f_plus[i], 5), "f_minus": round(f_minus[i], 5), "f_zero": round(f_zero[i], 5)}
        for i in range(n)
    ]


def _pka_group_table() -> List[Tuple[str, float, float, str, str]]:
    """Return list of (smarts, pka_low, pka_high, group_name, site_type) for heuristic pKa."""
    return [
        # acidic groups
        ("[CX3](=O)[OX2H1]",           3.5,  5.0, "Carboxylic acid",   "acidic"),
        ("c[OX2H]",                     9.0, 10.5, "Phenol",            "acidic"),
        ("[S](=O)(=O)[NX3H]",           9.5, 10.5, "Sulfonamide NH",    "acidic"),
        ("[SX2H]",                      9.5, 11.0, "Thiol",             "acidic"),
        ("[CX3](=O)[NX3H]",            15.0, 18.0, "Amide NH",          "acidic"),
        ("[CX4][OX2H]",                15.0, 18.0, "Alkyl alcohol",     "acidic"),
        # basic groups
        ("[NX3][CX3](=[NX2])[NX3]",   12.0, 13.5, "Guanidine",         "basic"),
        ("[CX3](=[NX2])[NX3]",        11.0, 12.5, "Amidine",           "basic"),
        ("c[NX3H2]",                   4.0,  5.0, "Aromatic amine",    "basic"),
        ("[CX4][NX3H2]",               9.5, 10.5, "Aliphatic 1° amine","basic"),
        ("[CX4][NX3H1][CX4]",         10.0, 11.0, "Aliphatic 2° amine","basic"),
        ("[NX3]([CX4])[CX4][CX4]",    9.5, 10.5, "Aliphatic 3° amine","basic"),
        ("[nX2;!H]",                   5.0,  6.0, "Pyridine-like N",   "basic"),
        ("[nH]",                        6.5,  7.5, "Imidazole NH",      "basic"),
    ]


def estimate_pka_groups(mol: Chem.Mol) -> List[dict]:
    """Return list of pKa estimates for ionisable groups via SMARTS matching.

    atom_idx is 1-based (matches structure display convention).
    """
    groups: List[dict] = []
    seen_atoms: set = set()
    for smarts, pka_low, pka_high, group_name, site_type in _pka_group_table():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern):
            # Ionisable atom: O for acids, N for amines, S for thiol
            ionisable_idx = None
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                anum = atom.GetAtomicNum()
                if site_type == "acidic" and anum in (7, 8, 16):
                    ionisable_idx = idx
                    break
                if site_type == "basic" and anum == 7:
                    ionisable_idx = idx
                    break
            if ionisable_idx is None:
                ionisable_idx = match[0]
            if ionisable_idx in seen_atoms:
                continue
            seen_atoms.add(ionisable_idx)
            pka_est = round((pka_low + pka_high) / 2.0, 1)
            groups.append({
                "atom_idx": ionisable_idx + 1,  # 1-based
                "pka_low": pka_low,
                "pka_high": pka_high,
                "pka_est": pka_est,
                "group_name": group_name,
                "site_type": site_type,
            })
    return groups


def estimate_redox_potentials(homo_ev: float, lumo_ev: float) -> dict:
    """Estimate oxidation/reduction potentials via Trasatti vacuum-to-NHE shift.

    E_ox_NHE  ≈ −HOMO_eV − 4.44   (V vs NHE)
    E_red_NHE ≈ −LUMO_eV − 4.44   (V vs NHE)
    Also reported vs Fc/Fc+ (Fc ≈ +0.40 V vs NHE).
    Warning: gas-phase estimates, errors typically ±0.5 V.
    """
    e_ox_nhe = round(-homo_ev - _VACUUM_TO_NHE_EV, 3)
    e_red_nhe = round(-lumo_ev - _VACUUM_TO_NHE_EV, 3)
    e_ox_fc = round(e_ox_nhe - _FC_VS_NHE_V, 3)
    e_red_fc = round(e_red_nhe - _FC_VS_NHE_V, 3)
    return {
        "e_ox_nhe": e_ox_nhe,
        "e_red_nhe": e_red_nhe,
        "e_ox_fc": e_ox_fc,
        "e_red_fc": e_red_fc,
        "warning": "Estimates based on gas-phase GFN2-xTB HOMO/LUMO; errors \u00b10.5 V typical.",
    }


def predict_molecular_properties(
    mol: Chem.Mol,
    selected_confs: Sequence,
    weights: Sequence[float],
    solvent: str,
    output_dir: Path,
    warnings: List[str],
    progress_path: Optional[Path],
) -> dict:
    """Orchestrate all property calculations and return a unified properties dict."""
    props: dict = {
        "rdkit": {},
        "electronic": {},
        "fukui": [],
        "pka_groups": [],
        "redox": {},
        "warnings": [],
        "used_xtb": False,
    }

    # (a) RDKit descriptors — always available
    try:
        props["rdkit"] = compute_rdkit_descriptors(mol)
    except Exception as exc:
        props["warnings"].append(f"RDKit descriptors failed: {exc}")

    # (b) Heuristic pKa — always available
    try:
        props["pka_groups"] = estimate_pka_groups(mol)
    except Exception as exc:
        props["warnings"].append(f"pKa estimation failed: {exc}")

    # (c) xTB electronic + Fukui
    xtb_bin_env = os.environ.get("EASYSPECTRA_XTB", "xtb")
    xtb_bin = None if xtb_bin_env == "__none__" else shutil.which(xtb_bin_env)
    timeout_s = int(os.environ.get("EASYSPECTRA_XTB_TIMEOUT", "25"))

    if not xtb_bin:
        props["warnings"].append(
            "xTB not available — skipping electronic properties and Fukui indices."
        )
        warnings.append("xTB not available; electronic properties and Fukui indices skipped.")
        return props

    target_conf = selected_confs[0]
    conf_id = target_conf.conf_id
    formal_charge = Chem.GetFormalCharge(mol)

    sp_dir = output_dir / "xtb_properties"
    sp_dir.mkdir(parents=True, exist_ok=True)

    if progress_path:
        write_progress(progress_path, "properties_sp", "Running xTB SP (N electrons)", 0.87)

    result_n = run_xtb_sp(xtb_bin, mol, conf_id, solvent, sp_dir, formal_charge, timeout_s, label="spN")
    if result_n is None:
        props["warnings"].append("xTB SP (N electrons) failed — skipping electronic properties and Fukui.")
        warnings.append("xTB SP run failed for properties; electronic properties skipped.")
        return props

    props["used_xtb"] = True
    props["electronic"] = {
        "homo_ev": round(result_n["homo_ev"], 4),
        "lumo_ev": round(result_n["lumo_ev"], 4),
        "gap_ev": round(result_n["gap_ev"], 4),
        "dipole_debye": round(result_n["dipole_debye"], 3),
    }

    # Redox estimates from HOMO/LUMO
    if result_n["homo_ev"] != 0.0 or result_n["lumo_ev"] != 0.0:
        try:
            props["redox"] = estimate_redox_potentials(result_n["homo_ev"], result_n["lumo_ev"])
        except Exception as exc:
            props["warnings"].append(f"Redox estimation failed: {exc}")

    # Fukui: N+1 and N-1 SP runs
    if progress_path:
        write_progress(progress_path, "properties_sp", "Running xTB SP (N+1 electrons)", 0.90)
    result_np1 = run_xtb_sp(xtb_bin, mol, conf_id, solvent, sp_dir, formal_charge - 1, timeout_s, label="spNp1")

    if progress_path:
        write_progress(progress_path, "properties_sp", "Running xTB SP (N-1 electrons)", 0.93)
    result_nm1 = run_xtb_sp(xtb_bin, mol, conf_id, solvent, sp_dir, formal_charge + 1, timeout_s, label="spNm1")

    if result_np1 is not None and result_nm1 is not None:
        try:
            fukui = compute_fukui_indices(
                result_n["charges"],
                result_np1["charges"],
                result_nm1["charges"],
            )
            mol_no_h = Chem.RemoveHs(mol)
            # Build element list from heavy-atom mol (indices may differ from full mol)
            # Use full mol atom ordering — xTB works on full mol with H
            for i, f in enumerate(fukui):
                atom = mol.GetAtomWithIdx(i)
                props["fukui"].append({
                    "atom_idx": i + 1,  # 1-based
                    "element": atom.GetSymbol(),
                    "f_plus": f["f_plus"],
                    "f_minus": f["f_minus"],
                    "f_zero": f["f_zero"],
                })
        except Exception as exc:
            props["warnings"].append(f"Fukui index computation failed: {exc}")
    else:
        props["warnings"].append(
            "xTB SP (N±1 electrons) failed — Fukui indices not available."
        )

    return props


def write_properties_json(path: Path, props: dict) -> None:
    """Write the properties dict to a structured JSON file."""
    path.write_text(json.dumps(props, indent=2, ensure_ascii=False), encoding="utf-8")


def write_spectrum_csv(path: Path, data: Sequence[Tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ppm", "intensity"])
        for ppm, intensity in data:
            writer.writerow([f"{ppm:.6f}", f"{intensity:.8f}"])


def write_cd_spectrum_csv(path: Path, data: Sequence[Tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "intensity"])
        for wavelength_nm, intensity in data:
            writer.writerow([f"{wavelength_nm:.6f}", f"{intensity:.8f}"])


def write_peaks_csv(path: Path, groups: Sequence[GroupPrediction], nucleus_symbol: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_id", "shift_ppm", "shift_hz", "multiplicity", "j_hz", "relative_integral", "assignment", "atom_indices"])
        for group in groups:
            atom_text = ",".join(str(idx) for idx in group.atom_indices)
            assignment = ",".join(f"{nucleus_symbol}{idx}" for idx in group.atom_indices)
            writer.writerow(
                [
                    group.group_id,
                    f"{group.shift_ppm:.4f}",
                    f"{group.shift_hz:.2f}",
                    group.multiplicity,
                    f"{group.j_hz:.2f}",
                    f"{group.integral:.2f}",
                    assignment,
                    atom_text,
                ]
            )


def write_cd_peaks_csv(path: Path, data: Sequence[Tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band_id", "wavelength_nm", "relative_intensity"])
        for idx, (wavelength_nm, intensity) in enumerate(_find_cd_band_extrema(data), start=1):
            writer.writerow([idx, f"{wavelength_nm:.4f}", f"{intensity:.5f}"])


def write_assignments(path_json: Path, path_csv: Path, groups: Sequence[GroupPrediction], nucleus_symbol: str) -> None:
    payload = {
        "mapping_mode": "grouped",
        "nucleus": nucleus_symbol,
        "groups": [
            {
                "group_id": group.group_id,
                "atom_indices": group.atom_indices,
                "center_ppm": round(group.shift_ppm, 4),
                "multiplicity": group.multiplicity,
            }
            for group in groups
        ],
    }
    path_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with path_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_id", "center_ppm", "multiplicity", "atom_indices"])
        for group in groups:
            writer.writerow(
                [
                    group.group_id,
                    f"{group.shift_ppm:.4f}",
                    group.multiplicity,
                    ",".join(str(idx) for idx in group.atom_indices),
                ]
            )


def write_empty_assignments(path_json: Path, path_csv: Path, label: str) -> None:
    payload = {
        "mapping_mode": "none",
        "label": label,
        "groups": [],
    }
    path_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    with path_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_id", "center", "label", "atom_indices"])


def write_spectra_manifest(
    path_csv: Path,
    rows: Sequence[Tuple[str, Path, Path, Path]],
) -> None:
    with path_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["spectrum_label", "spectrum_csv", "peaks_csv", "assignments_csv"])
        for nucleus_label, spectrum_csv, peaks_csv, assignments_csv in rows:
            writer.writerow([nucleus_label, str(spectrum_csv), str(peaks_csv), str(assignments_csv)])


# ---------------------------------------------------------------------------
# Reaction comparison helpers
# ---------------------------------------------------------------------------

def _resample_onto_axis(
    data: Sequence[Tuple[float, float]],
    ppm_axis: np.ndarray,
) -> np.ndarray:
    """Linearly interpolate a spectrum onto a uniform ppm_axis."""
    if not data:
        return np.zeros_like(ppm_axis)
    x = np.array([p[0] for p in data], dtype=float)
    y = np.array([p[1] for p in data], dtype=float)
    # np.interp requires ascending x; spectra are stored high→low, so flip.
    if len(x) > 1 and x[0] > x[-1]:
        x = x[::-1]
        y = y[::-1]
    return np.interp(ppm_axis, x, y, left=0.0, right=0.0)


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity in [0, 1]; 1 = identical spectra."""
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a <= 0.0 or norm_b <= 0.0:
        return 0.0
    return max(0.0, min(1.0, float(np.dot(a, b) / (norm_a * norm_b))))


def _match_peaks_greedy(
    reactant_groups: List[GroupPrediction],
    product_groups: List[GroupPrediction],
    tolerance_ppm: float,
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Greedy nearest-peak matching within tolerance.

    Returns (matched, lost, new_peaks):
      matched   — peaks present in both; status = "retained" or "shifted"
      lost      — reactant peaks with no product match within tolerance
      new_peaks — product peaks with no reactant match within tolerance
    """
    r_peaks = sorted(reactant_groups, key=lambda g: g.shift_ppm)
    p_peaks = sorted(product_groups, key=lambda g: g.shift_ppm)

    r_used = [False] * len(r_peaks)
    p_used = [False] * len(p_peaks)
    matched: List[dict] = []

    for r_idx, r_group in enumerate(r_peaks):
        best_p_idx: Optional[int] = None
        best_delta = float("inf")
        for p_idx, p_group in enumerate(p_peaks):
            if p_used[p_idx]:
                continue
            delta = abs(r_group.shift_ppm - p_group.shift_ppm)
            if delta <= tolerance_ppm and delta < best_delta:
                best_delta = delta
                best_p_idx = p_idx
        if best_p_idx is not None:
            r_used[r_idx] = True
            p_used[best_p_idx] = True
            signed_delta = r_peaks[r_idx].shift_ppm - p_peaks[best_p_idx].shift_ppm
            status = "retained" if abs(signed_delta) <= tolerance_ppm * 0.3 else "shifted"
            matched.append({
                "reactant_ppm": round(r_peaks[r_idx].shift_ppm, 4),
                "product_ppm": round(p_peaks[best_p_idx].shift_ppm, 4),
                "delta_ppm": round(signed_delta, 4),
                "status": status,
                "reactant_multiplicity": r_peaks[r_idx].multiplicity,
                "product_multiplicity": p_peaks[best_p_idx].multiplicity,
                "reactant_integral": round(r_peaks[r_idx].integral, 2),
                "product_integral": round(p_peaks[best_p_idx].integral, 2),
            })

    lost = [
        {
            "ppm": round(g.shift_ppm, 4),
            "multiplicity": g.multiplicity,
            "integral": round(g.integral, 2),
            "status": "lost",
        }
        for i, g in enumerate(r_peaks) if not r_used[i]
    ]
    new_peaks = [
        {
            "ppm": round(g.shift_ppm, 4),
            "multiplicity": g.multiplicity,
            "integral": round(g.integral, 2),
            "status": "new",
        }
        for i, g in enumerate(p_peaks) if not p_used[i]
    ]
    return matched, lost, new_peaks


def compute_spectral_comparison(
    reactant_groups: List[GroupPrediction],
    product_groups: List[GroupPrediction],
    reactant_spectrum: List[Tuple[float, float]],
    product_spectrum: List[Tuple[float, float]],
    nucleus: str,
) -> dict:
    """Compare two NMR predictions; return similarity score and peak diff."""
    tolerance = COMPARE_TOLERANCE_PPM.get(nucleus, 0.5)
    matched, lost, new_peaks = _match_peaks_greedy(reactant_groups, product_groups, tolerance)

    # Build a common ppm axis for spectral similarity.
    all_ppm = [p[0] for p in reactant_spectrum] + [p[0] for p in product_spectrum]
    if all_ppm:
        ppm_axis = np.linspace(max(all_ppm), min(all_ppm), 3000)
        r_vec = _resample_onto_axis(reactant_spectrum, ppm_axis)
        p_vec = _resample_onto_axis(product_spectrum, ppm_axis)
        similarity = _cosine_similarity(r_vec, p_vec)
    else:
        similarity = 0.0

    n_retained = sum(1 for m in matched if m["status"] == "retained")
    n_shifted = sum(1 for m in matched if m["status"] == "shifted")
    indicator = _reaction_indicator(similarity, len(lost), len(new_peaks), n_shifted)

    return {
        "spectral_similarity_pct": round(similarity * 100.0, 1),
        "retained_peaks": n_retained,
        "shifted_peaks": n_shifted,
        "lost_peaks": len(lost),
        "new_peaks": len(new_peaks),
        "reaction_indicator": indicator,
        "matched_peaks": matched,
        "lost_peak_details": lost,
        "new_peak_details": new_peaks,
    }


def _reaction_indicator(similarity: float, n_lost: int, n_new: int, n_shifted: int) -> str:
    if similarity >= 0.92 and n_lost == 0 and n_new == 0 and n_shifted == 0:
        return "no_reaction"
    if similarity < 0.60 or (n_lost + n_new) >= 2:
        return "reaction_likely"
    return "inconclusive"


def _overall_reaction_indicator(results_by_nucleus: Dict[str, dict]) -> str:
    if not results_by_nucleus:
        return "inconclusive"
    indicators = [r.get("reaction_indicator", "inconclusive") for r in results_by_nucleus.values()]
    if any(i == "reaction_likely" for i in indicators):
        return "reaction_likely"
    if all(i == "no_reaction" for i in indicators):
        return "no_reaction"
    return "inconclusive"


def _reaction_description(indicator: str) -> str:
    return {
        "reaction_likely": (
            "Significant spectral differences detected — reaction likely occurred."
        ),
        "no_reaction": (
            "Spectra are highly similar — reaction likely did not occur."
        ),
        "inconclusive": (
            "Moderate spectral differences — result inconclusive; inspect overlays manually."
        ),
    }.get(indicator, "Unknown indicator.")


def write_compare_spectrum_csv(
    path: Path,
    reactant_data: Sequence[Tuple[float, float]],
    product_data: Sequence[Tuple[float, float]],
) -> None:
    """Write a 3-column overlay CSV: ppm, reactant_intensity, product_intensity."""
    if not reactant_data and not product_data:
        return
    all_ppm = [p[0] for p in reactant_data] + [p[0] for p in product_data]
    ppm_axis = np.linspace(max(all_ppm), min(all_ppm), 6000)
    r_vals = _resample_onto_axis(reactant_data, ppm_axis)
    p_vals = _resample_onto_axis(product_data, ppm_axis)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ppm", "reactant_intensity", "product_intensity"])
        for ppm, r, p in zip(ppm_axis.tolist(), r_vals.tolist(), p_vals.tolist()):
            writer.writerow([f"{ppm:.6f}", f"{r:.8f}", f"{p:.8f}"])


def write_peaks_diff_csv(
    path: Path,
    matched: List[dict],
    lost: List[dict],
    new_peaks: List[dict],
) -> None:
    """Write peak-diff table with status annotations."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "status", "reactant_ppm", "product_ppm", "delta_ppm",
            "reactant_multiplicity", "product_multiplicity",
            "reactant_integral", "product_integral",
        ])
        for m in matched:
            writer.writerow([
                m["status"], m["reactant_ppm"], m["product_ppm"], m["delta_ppm"],
                m["reactant_multiplicity"], m["product_multiplicity"],
                m["reactant_integral"], m["product_integral"],
            ])
        for row in lost:
            writer.writerow([
                "lost", row["ppm"], "", "",
                row.get("multiplicity", ""), "",
                row.get("integral", ""), "",
            ])
        for row in new_peaks:
            writer.writerow([
                "new", "", row["ppm"], "",
                "", row.get("multiplicity", ""),
                "", row.get("integral", ""),
            ])


def main() -> int:
    args = parse_args()

    request_path = Path(args.request)
    response_path = Path(args.response)
    request = json.loads(request_path.read_text(encoding="utf-8"))

    job_id = str(request.get("job_id", "job-unknown"))
    mode = str(request.get("mode", "predict")).lower()
    workflow_obj = request.get("workflow", {})
    if isinstance(workflow_obj, dict):
        workflow_kind = str(workflow_obj.get("kind", "all")).strip().lower()
    else:
        workflow_kind = str(workflow_obj).strip().lower()
    if not workflow_kind:
        workflow_kind = "all"
    input_obj = request.get("input", {})
    settings = request.get("settings", {})

    input_format = str(input_obj.get("format", "unknown")).lower()
    input_value = str(input_obj.get("value", ""))

    frequency_mhz = float(settings.get("frequency_mhz", 400.0))
    line_shape = str(settings.get("line_shape", "lorentzian")).lower()
    fwhm_hz = float(settings.get("fwhm_hz", 1.0))
    max_conformers = int(settings.get("max_conformers", 20))
    boltzmann_cutoff = float(settings.get("boltzmann_cutoff", 0.99))
    energy_window_kcal = float(settings.get("energy_window_kcal", 6.0))
    nucleus_request = normalize_nucleus(str(settings.get("nucleus", "auto")))
    need_editable_xyz = parse_bool(settings.get("need_editable_xyz", mode != "preview"), mode != "preview")

    solvent_input = str(settings.get("solvent", "cdcl3")).lower()
    solvent = SOLVENT_ALIASES.get(solvent_input, solvent_input)

    warnings: List[str] = []
    output_dir = Path(str(request.get("output_dir", request_path.parent)))
    output_dir.mkdir(parents=True, exist_ok=True)
    progress_path = Path(str(request.get("progress_json", output_dir / "progress.json")))
    write_progress(progress_path, "initializing", "Starting easySpectra workflow", 0.01)

    if workflow_kind not in {"all", "nmr", "cd", "ir", "uvvis", "compare", "properties"}:
        message = (
            f"Unsupported workflow kind '{workflow_kind}'. "
            "Currently available workflows: all, nmr, cd, ir, uvvis, compare, properties."
        )
        write_progress(progress_path, "failed", message, 1.0)
        response = {
            "status": "failed",
            "message": message,
            "workflow": {"kind": workflow_kind},
            "output_dir": str(output_dir),
            "progress_json": str(progress_path),
            "warnings": [],
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    start_time = time.time()

    write_progress(progress_path, "parse_input", "Parsing input and validating structure", 0.06)
    try:
        mol = load_molecule(input_format, input_value, warnings)
    except Exception as exc:
        write_progress(progress_path, "failed", f"Input parsing failed: {exc}", 1.0)
        response = {
            "status": "failed",
            "message": f"Input parsing failed: {exc}",
            "output_dir": str(output_dir),
            "progress_json": str(progress_path),
            "warnings": warnings,
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    structure_svg = output_dir / "structure.svg"
    structure_atoms_csv = output_dir / "structure_atoms.csv"
    structure_bonds_csv = output_dir / "structure_bonds.csv"
    structure_xyz: Optional[Path] = None
    write_progress(progress_path, "structure_2d", "Generating 2D depiction and atom mapping", 0.14)
    try:
        write_structure_svg(mol, structure_svg)
        write_structure_geometry_csv(mol, structure_atoms_csv, structure_bonds_csv)
    except Exception:
        warnings.append("Could not generate 2D structure visualization files.")
    if need_editable_xyz:
        structure_xyz = output_dir / "structure_edit.xyz"
        try:
            write_editable_xyz(mol, structure_xyz, warnings)
        except Exception:
            warnings.append("Could not generate editable XYZ export for structure editor.")

    if mode == "preview":
        write_progress(progress_path, "done", "Preview complete", 1.0)
        response = {
            "status": "ok",
            "message": "Preview ready",
            "workflow": {"kind": workflow_kind},
            "output_dir": str(output_dir),
            "structure_svg": str(structure_svg),
            "structure_atoms_csv": str(structure_atoms_csv),
            "structure_bonds_csv": str(structure_bonds_csv),
            "structure_xyz": str(structure_xyz) if structure_xyz is not None else "",
            "progress_json": str(progress_path),
            "warnings": warnings,
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    # ── Compare workflow ─────────────────────────────────────────────────────
    if workflow_kind == "compare":
        compare_input_obj = request.get("compare_input", {})
        compare_format = str(compare_input_obj.get("format", "smiles")).lower()
        compare_value = str(compare_input_obj.get("value", ""))

        if not compare_value:
            message = "compare workflow requires 'compare_input.value' in the request."
            write_progress(progress_path, "failed", message, 1.0)
            response_path.write_text(
                json.dumps({"status": "failed", "message": message,
                            "output_dir": str(output_dir), "warnings": warnings}, indent=2),
                encoding="utf-8",
            )
            return 0

        write_progress(progress_path, "parse_input", "Parsing product structure", 0.08)
        warnings_product: List[str] = []
        try:
            mol_product = load_molecule(compare_format, compare_value, warnings_product)
        except Exception as exc:
            message = f"Product input parsing failed: {exc}"
            write_progress(progress_path, "failed", message, 1.0)
            response_path.write_text(
                json.dumps({"status": "failed", "message": message,
                            "output_dir": str(output_dir), "warnings": warnings}, indent=2),
                encoding="utf-8",
            )
            return 0

        # Write product 2D structure.
        struct_svg_product = output_dir / "structure_product.svg"
        struct_atoms_product = output_dir / "structure_atoms_product.csv"
        struct_bonds_product = output_dir / "structure_bonds_product.csv"
        write_progress(progress_path, "structure_2d", "Generating product 2D depiction", 0.12)
        try:
            write_structure_svg(mol_product, struct_svg_product)
            write_structure_geometry_csv(mol_product, struct_atoms_product, struct_bonds_product)
        except Exception:
            warnings_product.append("Could not generate 2D structure visualization for product.")

        # Determine which nuclei to compare.
        nuclei_order = ["1h", "13c", "19f", "31p"]
        available_reactant = set(nuclei_present(mol))
        available_product = set(nuclei_present(mol_product))
        if nucleus_request == "auto":
            all_nuclei = sorted(
                available_reactant | available_product,
                key=lambda n: nuclei_order.index(n) if n in nuclei_order else 99,
            )
        else:
            all_nuclei = [nucleus_request]

        # Conformer generation and optimisation — reactant.
        seed_r = hash_seed(job_id, input_value)
        seed_p = hash_seed(job_id, compare_value)
        write_progress(progress_path, "conformer_generation", "Embedding reactant conformers", 0.18)
        conf_ids_r = embed_conformers(mol, max_conformers=max_conformers, seed=seed_r, warnings=warnings)
        results_r = optimize_conformers(mol, conf_ids_r, solvent, output_dir / "xtb_reactant",
                                        warnings, progress_path=None)
        sel_r, w_r, rel_r = select_conformers(results_r, energy_window_kcal, boltzmann_cutoff)

        # Conformer generation and optimisation — product.
        write_progress(progress_path, "conformer_generation", "Embedding product conformers", 0.45)
        conf_ids_p = embed_conformers(mol_product, max_conformers=max_conformers, seed=seed_p,
                                      warnings=warnings_product)
        results_p = optimize_conformers(mol_product, conf_ids_p, solvent, output_dir / "xtb_product",
                                        warnings_product, progress_path=None)
        sel_p, w_p, rel_p = select_conformers(results_p, energy_window_kcal, boltzmann_cutoff)

        spectra_rows: List[Tuple[str, Path, Path, Path]] = []
        comparison_results: Dict[str, dict] = {}

        write_progress(progress_path, "nmr_parameter_estimation",
                       "Predicting and comparing NMR spectra", 0.70)

        for nucleus in all_nuclei:
            has_r = nucleus in available_reactant
            has_p = nucleus in available_product
            nucleus_label = NUCLEUS_LABEL.get(nucleus, nucleus.upper())
            nucleus_symbol = NUCLEUS_SYMBOL.get(nucleus, "H")

            groups_r: List[GroupPrediction] = []
            groups_p: List[GroupPrediction] = []
            spectrum_r: List[Tuple[float, float]] = []
            spectrum_p: List[Tuple[float, float]] = []

            if has_r:
                groups_r = build_group_predictions(mol, sel_r, w_r, frequency_mhz, nucleus)
                spectrum_r = simulate_spectrum(groups_r, frequency_mhz, line_shape, fwhm_hz, nucleus)
                spec_csv_r = output_dir / f"spectrum_{nucleus}_reactant.csv"
                peaks_csv_r = output_dir / f"peaks_{nucleus}_reactant.csv"
                asgn_json_r = output_dir / f"assignments_{nucleus}_reactant.json"
                asgn_csv_r = output_dir / f"assignments_{nucleus}_reactant.csv"
                write_spectrum_csv(spec_csv_r, spectrum_r)
                write_peaks_csv(peaks_csv_r, groups_r, nucleus_symbol)
                write_assignments(asgn_json_r, asgn_csv_r, groups_r, nucleus_symbol)
                spectra_rows.append((f"{nucleus_label} (Reactant)", spec_csv_r,
                                     peaks_csv_r, asgn_csv_r))

            if has_p:
                groups_p = build_group_predictions(mol_product, sel_p, w_p, frequency_mhz, nucleus)
                spectrum_p = simulate_spectrum(groups_p, frequency_mhz, line_shape, fwhm_hz, nucleus)
                spec_csv_p = output_dir / f"spectrum_{nucleus}_product.csv"
                peaks_csv_p = output_dir / f"peaks_{nucleus}_product.csv"
                asgn_json_p = output_dir / f"assignments_{nucleus}_product.json"
                asgn_csv_p = output_dir / f"assignments_{nucleus}_product.csv"
                write_spectrum_csv(spec_csv_p, spectrum_p)
                write_peaks_csv(peaks_csv_p, groups_p, nucleus_symbol)
                write_assignments(asgn_json_p, asgn_csv_p, groups_p, nucleus_symbol)
                spectra_rows.append((f"{nucleus_label} (Product)", spec_csv_p,
                                     peaks_csv_p, asgn_csv_p))

            if has_r and has_p and spectrum_r and spectrum_p:
                cmp = compute_spectral_comparison(groups_r, groups_p, spectrum_r, spectrum_p, nucleus)
                compare_spec_csv = output_dir / f"spectrum_compare_{nucleus}.csv"
                peaks_diff_csv = output_dir / f"peaks_diff_{nucleus}.csv"
                write_compare_spectrum_csv(compare_spec_csv, spectrum_r, spectrum_p)
                write_peaks_diff_csv(peaks_diff_csv, cmp["matched_peaks"],
                                     cmp["lost_peak_details"], cmp["new_peak_details"])
                comparison_results[nucleus_label] = cmp
            elif has_r and not has_p:
                comparison_results[nucleus_label] = {
                    "spectral_similarity_pct": 0.0,
                    "retained_peaks": 0, "shifted_peaks": 0,
                    "lost_peaks": len(groups_r), "new_peaks": 0,
                    "reaction_indicator": "reaction_likely",
                    "note": f"{nucleus_label} present in reactant but absent in product.",
                }
            elif has_p and not has_r:
                comparison_results[nucleus_label] = {
                    "spectral_similarity_pct": 0.0,
                    "retained_peaks": 0, "shifted_peaks": 0,
                    "lost_peaks": 0, "new_peaks": len(groups_p),
                    "reaction_indicator": "reaction_likely",
                    "note": f"{nucleus_label} absent in reactant but present in product.",
                }

        spectra_manifest_csv = output_dir / "spectra_manifest.csv"
        write_spectra_manifest(spectra_manifest_csv, spectra_rows)

        overall = _overall_reaction_indicator(comparison_results)
        reaction_summary = {
            "workflow": "compare",
            "reactant_formula": Chem.rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(mol)),
            "product_formula": Chem.rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(mol_product)),
            "nuclei": list(comparison_results.keys()),
            "results": comparison_results,
            "overall_indicator": overall,
            "overall_description": _reaction_description(overall),
        }
        reaction_summary_json = output_dir / "reaction_summary.json"
        reaction_summary_json.write_text(json.dumps(reaction_summary, indent=2), encoding="utf-8")

        write_progress(progress_path, "done", "Comparison complete", 1.0)
        response = {
            "status": "ok",
            "message": "Comparison complete",
            "workflow": {"kind": "compare"},
            "output_dir": str(output_dir),
            "spectra_manifest_csv": str(spectra_manifest_csv),
            "reaction_summary_json": str(reaction_summary_json),
            "structure_svg": str(structure_svg),
            "structure_atoms_csv": str(structure_atoms_csv),
            "structure_bonds_csv": str(structure_bonds_csv),
            "structure_product_svg": str(struct_svg_product),
            "structure_atoms_product_csv": str(struct_atoms_product),
            "structure_bonds_product_csv": str(struct_bonds_product),
            "progress_json": str(progress_path),
            "warnings": warnings + warnings_product,
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    run_nmr = workflow_kind in {"all", "nmr"}
    run_cd = workflow_kind in {"all", "cd"}
    run_ir = workflow_kind in {"all", "ir"}
    run_uvvis = workflow_kind in {"all", "uvvis"}
    run_properties = workflow_kind in {"all", "properties"}

    target_nuclei: List[str] = []
    if run_nmr:
        available_nuclei = nuclei_present(mol)
        if nucleus_request == "auto":
            target_nuclei = available_nuclei
            if not target_nuclei:
                if workflow_kind == "nmr":
                    write_progress(progress_path, "failed", "No supported nuclei (1H/13C/19F/31P) found", 1.0)
                    response = {
                        "status": "failed",
                        "message": "No supported nuclei (1H/13C/19F/31P) present in the molecule.",
                        "output_dir": str(output_dir),
                        "structure_svg": str(structure_svg),
                        "structure_atoms_csv": str(structure_atoms_csv),
                        "structure_bonds_csv": str(structure_bonds_csv),
                        "structure_xyz": str(structure_xyz) if structure_xyz is not None else "",
                        "progress_json": str(progress_path),
                        "warnings": warnings,
                    }
                    response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
                    return 0
                warnings.append(
                    "No supported nuclei (1H/13C/19F/31P) present; skipped NMR products for workflow='all'."
                )
                run_nmr = False
        else:
            target_nuclei = [nucleus_request]
            target_atomic_number = {"1h": 1, "13c": 6, "19f": 9, "31p": 15}.get(nucleus_request, 1)
            target_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == target_atomic_number)
            if target_count == 0:
                label = NUCLEUS_LABEL.get(nucleus_request, nucleus_request.upper())
                if workflow_kind == "nmr":
                    write_progress(progress_path, "failed", f"No {label} nuclei found in this molecule", 1.0)
                    response = {
                        "status": "failed",
                        "message": f"No {label} nuclei present in the molecule.",
                        "output_dir": str(output_dir),
                        "structure_svg": str(structure_svg),
                        "structure_atoms_csv": str(structure_atoms_csv),
                        "structure_bonds_csv": str(structure_bonds_csv),
                        "structure_xyz": str(structure_xyz) if structure_xyz is not None else "",
                        "progress_json": str(progress_path),
                        "warnings": warnings,
                    }
                    response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
                    return 0
                warnings.append(
                    f"No {label} nuclei present; skipped NMR products for workflow='all'."
                )
                run_nmr = False
                target_nuclei = []

    seed = hash_seed(job_id, input_value)
    write_progress(
        progress_path,
        "conformer_generation",
        f"Embedding up to {max_conformers} conformers using ETKDG + MMFF pre-optimization",
        0.24,
    )
    conformer_ids = embed_conformers(mol, max_conformers=max_conformers, seed=seed, warnings=warnings)
    write_progress(progress_path, "conformer_generation", f"Generated {len(conformer_ids)} conformers", 0.32)

    write_progress(
        progress_path,
        "conformer_optimization",
        f"Starting conformer optimization with xTB GFN2-xTB + ALPB({solvent}), MMFF fallback",
        0.36,
    )
    conformer_results = optimize_conformers(
        mol,
        conformer_ids,
        solvent=solvent,
        output_dir=output_dir,
        warnings=warnings,
        progress_path=progress_path,
    )
    write_progress(progress_path, "conformer_optimization", "Conformer optimization complete", 0.70)

    write_progress(
        progress_path,
        "boltzmann_selection",
        (
            f"Applying Boltzmann weighting (cutoff={boltzmann_cutoff:.2f}, "
            f"energy window={energy_window_kcal:.1f} kcal/mol)"
        ),
        0.76,
    )
    selected_confs, conf_weights, rel_energies = select_conformers(
        conformer_results,
        energy_window_kcal=energy_window_kcal,
        boltzmann_cutoff=boltzmann_cutoff,
    )

    spectra_rows: List[Tuple[str, Path, Path, Path]] = []
    output_by_label: Dict[str, Dict[str, Path]] = {}
    produced_products: List[str] = []
    default_label = (
        "CD" if run_cd and not run_nmr and not run_ir and not run_uvvis
        else "UV-Vis" if run_uvvis and not run_nmr and not run_cd and not run_ir
        else "1H"
    )

    if run_nmr:
        total_targets = max(1, len(target_nuclei))
        nmr_progress_start = 0.80 if run_cd else 0.82
        nmr_progress_span = 0.08 if run_cd else 0.06
        for idx, target_nucleus in enumerate(target_nuclei, start=1):
            nucleus_label = NUCLEUS_LABEL.get(target_nucleus, target_nucleus.upper())
            nucleus_symbol = NUCLEUS_SYMBOL.get(target_nucleus, "H")
            frac_base = float(idx - 1) / float(total_targets)
            frac_next = float(idx) / float(total_targets)

            write_progress(
                progress_path,
                "nmr_parameter_estimation",
                (
                    f"[{idx}/{total_targets}] Estimating {nucleus_label} shifts"
                    " and couplings in first-order approximation"
                ),
                nmr_progress_start + nmr_progress_span * frac_base,
            )
            group_predictions = build_group_predictions(
                mol,
                selected_confs,
                conf_weights,
                frequency_mhz=frequency_mhz,
                nucleus=target_nucleus,
            )

            write_progress(
                progress_path,
                "spectrum_simulation",
                f"[{idx}/{total_targets}] Simulating {nucleus_label} spectrum",
                nmr_progress_start + nmr_progress_span * (frac_base + frac_next) / 2.0,
            )
            spectrum = simulate_spectrum(group_predictions, frequency_mhz, line_shape, fwhm_hz, target_nucleus)

            suffix = target_nucleus
            spectrum_csv = output_dir / f"spectrum_{suffix}.csv"
            peaks_csv = output_dir / f"peaks_{suffix}.csv"
            assignments_json = output_dir / f"assignments_{suffix}.json"
            assignments_csv = output_dir / f"assignments_{suffix}.csv"

            write_spectrum_csv(spectrum_csv, spectrum)
            write_peaks_csv(peaks_csv, group_predictions, nucleus_symbol)
            write_assignments(assignments_json, assignments_csv, group_predictions, nucleus_symbol)

            spectra_rows.append((nucleus_label, spectrum_csv, peaks_csv, assignments_csv))
            output_by_label[nucleus_label] = {
                "spectrum_csv": spectrum_csv,
                "peaks_csv": peaks_csv,
                "assignments_json": assignments_json,
                "assignments_csv": assignments_csv,
            }
            produced_products.append(nucleus_label)
        if "1H" in output_by_label:
            default_label = "1H"
        elif output_by_label:
            default_label = spectra_rows[0][0]

    if run_cd:
        write_progress(
            progress_path,
            "cd_parameter_estimation",
            "Estimating CD band strengths and signs from conformer ensemble",
            0.90 if run_nmr else 0.84,
        )
        spectrum = simulate_cd_spectrum(mol, selected_confs, conf_weights)
        write_progress(
            progress_path,
            "spectrum_simulation",
            "Simulating CD spectrum",
            0.93 if run_nmr else 0.89,
        )

        spectrum_csv = output_dir / "spectrum_cd.csv"
        peaks_csv = output_dir / "peaks_cd.csv"
        assignments_json = output_dir / "assignments_cd.json"
        assignments_csv = output_dir / "assignments_cd.csv"

        write_cd_spectrum_csv(spectrum_csv, spectrum)
        write_cd_peaks_csv(peaks_csv, spectrum)
        write_empty_assignments(assignments_json, assignments_csv, "CD")

        spectra_rows.append(("CD", spectrum_csv, peaks_csv, assignments_csv))
        output_by_label["CD"] = {
            "spectrum_csv": spectrum_csv,
            "peaks_csv": peaks_csv,
            "assignments_json": assignments_json,
            "assignments_csv": assignments_csv,
        }
        produced_products.append("CD")
        if not run_nmr or default_label not in output_by_label:
            default_label = "CD"

    if run_ir:
        write_progress(
            progress_path,
            "ir_hessian",
            "Running xTB Hessian for IR prediction (Boltzmann-weighted conformers)",
            0.84 if not run_nmr and not run_cd else 0.91,
        )
        ir_spectrum, ir_labeled_peaks, ir_used_heuristic = predict_ir_spectrum(
            mol,
            selected_confs,
            conf_weights,
            solvent=solvent,
            output_dir=output_dir,
            warnings=warnings,
            progress_path=progress_path,
        )
        write_progress(progress_path, "ir_spectrum_simulation", "Simulating IR spectrum", 0.92 if not run_nmr and not run_cd else 0.95)

        spectrum_csv_ir = output_dir / "spectrum_ir.csv"
        peaks_csv_ir = output_dir / "peaks_ir.csv"
        assignments_json_ir = output_dir / "assignments_ir.json"
        assignments_csv_ir = output_dir / "assignments_ir.csv"

        write_ir_spectrum_csv(spectrum_csv_ir, ir_spectrum)
        write_ir_peaks_csv(peaks_csv_ir, ir_labeled_peaks)
        write_empty_assignments(assignments_json_ir, assignments_csv_ir, "IR")

        spectra_rows.append(("IR", spectrum_csv_ir, peaks_csv_ir, assignments_csv_ir))
        output_by_label["IR"] = {
            "spectrum_csv": spectrum_csv_ir,
            "peaks_csv": peaks_csv_ir,
            "assignments_json": assignments_json_ir,
            "assignments_csv": assignments_csv_ir,
        }
        produced_products.append("IR")
        if not run_nmr and not run_cd or default_label not in output_by_label:
            default_label = "IR"

    # ── UV-Vis workflow ───────────────────────────────────────────────────────
    if run_uvvis:
        write_progress(
            progress_path,
            "uvvis_stda",
            "Running xTB sTDA for UV-Vis prediction (Boltzmann-weighted conformers)",
            0.84 if not run_nmr and not run_cd and not run_ir else 0.93,
        )
        uvvis_spectrum, uvvis_labeled_peaks, uvvis_used_heuristic = predict_uvvis_spectrum(
            mol,
            selected_confs,
            conf_weights,
            solvent=solvent,
            output_dir=output_dir,
            warnings=warnings,
            progress_path=progress_path,
        )
        write_progress(progress_path, "uvvis_spectrum_simulation", "Simulating UV-Vis spectrum", 0.94 if not run_nmr and not run_cd and not run_ir else 0.96)

        spectrum_csv_uvvis = output_dir / "spectrum_uvvis.csv"
        peaks_csv_uvvis = output_dir / "peaks_uvvis.csv"
        assignments_json_uvvis = output_dir / "assignments_uvvis.json"
        assignments_csv_uvvis = output_dir / "assignments_uvvis.csv"

        write_uvvis_spectrum_csv(spectrum_csv_uvvis, uvvis_spectrum)
        write_uvvis_peaks_csv(peaks_csv_uvvis, uvvis_labeled_peaks)
        write_empty_assignments(assignments_json_uvvis, assignments_csv_uvvis, "UV-Vis")

        spectra_rows.append(("UV-Vis", spectrum_csv_uvvis, peaks_csv_uvvis, assignments_csv_uvvis))
        output_by_label["UV-Vis"] = {
            "spectrum_csv": spectrum_csv_uvvis,
            "peaks_csv": peaks_csv_uvvis,
            "assignments_json": assignments_json_uvvis,
            "assignments_csv": assignments_csv_uvvis,
        }
        produced_products.append("UV-Vis")
        if not run_nmr and not run_cd and not run_ir or default_label not in output_by_label:
            default_label = "UV-Vis"

    # ── Properties workflow ───────────────────────────────────────────────────
    properties_json_path: Optional[Path] = None
    if run_properties:
        write_progress(
            progress_path,
            "properties_sp",
            "Computing molecular properties (RDKit + xTB SP + Fukui + pKa)",
            0.84 if not run_nmr and not run_cd and not run_ir else 0.95,
        )
        try:
            props = predict_molecular_properties(
                mol,
                selected_confs,
                conf_weights,
                solvent=solvent,
                output_dir=output_dir,
                warnings=warnings,
                progress_path=progress_path,
            )
            properties_json_path = output_dir / "properties.json"
            write_properties_json(properties_json_path, props)
        except Exception as exc:
            warnings.append(f"Properties calculation failed: {exc}")
            properties_json_path = None

    # For a properties-only workflow, skip the spectra output check
    if workflow_kind == "properties":
        write_progress(progress_path, "done", "Properties prediction complete", 1.0)
        response = {
            "status": "ok",
            "message": "Properties prediction complete",
            "workflow": {"kind": workflow_kind},
            "output_dir": str(output_dir),
            "properties_json": str(properties_json_path) if properties_json_path else "",
            "structure_svg": str(structure_svg),
            "structure_atoms_csv": str(structure_atoms_csv),
            "structure_bonds_csv": str(structure_bonds_csv),
            "structure_xyz": str(structure_xyz) if structure_xyz is not None else "",
            "progress_json": str(progress_path),
            "warnings": warnings,
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    if not output_by_label:
        message = "No spectra products were generated for this workflow configuration."
        write_progress(progress_path, "failed", message, 1.0)
        response = {
            "status": "failed",
            "message": message,
            "workflow": {"kind": workflow_kind},
            "output_dir": str(output_dir),
            "progress_json": str(progress_path),
            "warnings": warnings,
        }
        response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
        return 0

    spectrum_csv = output_by_label[default_label]["spectrum_csv"]
    peaks_csv = output_by_label[default_label]["peaks_csv"]
    assignments_json = output_by_label[default_label]["assignments_json"]
    assignments_csv = output_by_label[default_label]["assignments_csv"]
    spectra_manifest_csv = output_dir / "spectra_manifest.csv"
    audit_json = output_dir / "audit.json"

    write_progress(progress_path, "write_outputs", "Writing spectrum, peak groups, assignments, and audit", 0.96)
    write_spectra_manifest(spectra_manifest_csv, spectra_rows)

    runtime_s = time.time() - start_time
    if workflow_kind == "cd":
        pipeline_summary = {
            "workflow_kind": workflow_kind,
            "geometry": "ETKDG + xTB(opt) with MMFF fallback",
            "simulation_mode": "ensemble CD scaffold",
            "products": ["CD"],
        }
    elif workflow_kind == "nmr":
        pipeline_summary = {
            "workflow_kind": workflow_kind,
            "geometry": "ETKDG + xTB(opt) with MMFF fallback",
            "simulation_mode": "first-order",
            "nucleus": default_label,
            "nuclei": [NUCLEUS_LABEL.get(n, n.upper()) for n in target_nuclei],
            "products": [NUCLEUS_LABEL.get(n, n.upper()) for n in target_nuclei],
        }
    elif workflow_kind == "ir":
        pipeline_summary = {
            "workflow_kind": workflow_kind,
            "geometry": "ETKDG + xTB(opt) with MMFF fallback",
            "simulation_mode": "xTB Hessian / force-field heuristic",
            "products": ["IR"],
        }
    elif workflow_kind == "uvvis":
        pipeline_summary = {
            "workflow_kind": workflow_kind,
            "geometry": "ETKDG + xTB(opt) with MMFF fallback",
            "simulation_mode": "xTB sTDA / chromophore heuristic",
            "products": ["UV-Vis"],
        }
    else:
        nmr_products = [label for label in produced_products if label in {"1H", "13C", "19F", "31P"}]
        pipeline_summary = {
            "workflow_kind": workflow_kind,
            "geometry": "ETKDG + xTB(opt) with MMFF fallback",
            "simulation_mode": "first-order (NMR) + ensemble CD scaffold",
            "nuclei": nmr_products,
            "products": produced_products,
        }
    audit = {
        "tool": "easySpectra",
        "backend_version": "0.0.1",
        "job_id": job_id,
        "timestamp_unix": int(time.time()),
        "input": {
            "format": input_format,
            "preview": input_value[:120],
            "formula": Chem.rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(mol)),
            "formal_charge": Chem.GetFormalCharge(mol),
        },
        "settings": settings,
        "workflow": {"kind": workflow_kind},
        "solvent_xtb": solvent,
        "runtime_seconds": round(runtime_s, 4),
        "conformers": {
            "generated": len(conformer_ids),
            "optimized": len(conformer_results),
            "selected": len(selected_confs),
            "selected_relative_kcal": [round(e, 4) for e in rel_energies],
            "selected_weights": [round(float(w), 6) for w in conf_weights],
            "methods": [c.method for c in selected_confs],
        },
        "warnings": warnings,
        "pipeline": pipeline_summary,
    }
    audit_json.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    write_progress(progress_path, "done", "Prediction complete", 1.0)

    response = {
        "status": "ok",
        "message": "Prediction complete",
        "workflow": {"kind": workflow_kind},
        "output_dir": str(output_dir),
        "spectrum_csv": str(spectrum_csv),
        "peaks_csv": str(peaks_csv),
        "assignments_json": str(assignments_json),
        "assignments_csv": str(assignments_csv),
        "spectra_manifest_csv": str(spectra_manifest_csv),
        "audit_json": str(audit_json),
        "structure_svg": str(structure_svg),
        "structure_atoms_csv": str(structure_atoms_csv),
        "structure_bonds_csv": str(structure_bonds_csv),
        "structure_xyz": str(structure_xyz) if structure_xyz is not None else "",
        "properties_json": str(properties_json_path) if properties_json_path else "",
        "progress_json": str(progress_path),
        "warnings": warnings,
    }
    response_path.write_text(json.dumps(response, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
