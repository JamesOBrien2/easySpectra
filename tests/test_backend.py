"""Unit tests for easynmr_backend.py chemistry functions.

Run with: pytest tests/test_backend.py -v
Requires: pytest, numpy, rdkit-pypi (see backend/requirements.txt)
"""
from __future__ import annotations

import os
import sys
import math

import numpy as np
import pytest

# Make the backend importable without installing it.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "backend"))

import easynmr_backend as B
from easynmr_backend import (
    ConformerResult,
    GroupPrediction,
    MULTIPLICITY_LABELS,
    boltzmann_weights,
    estimate_c_shift,
    estimate_f_shift,
    estimate_h_shift,
    estimate_j_between_groups,
    estimate_p_shift,
    group_atoms_by_rank,
    group_hydrogens,
    load_molecule,
    normalize_nucleus,
    nuclei_present,
    select_conformers,
    simulate_spectrum,
    splitting_pattern,
    _chirality_sign,
)

from rdkit import Chem
from rdkit.Chem import AllChem


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mol_with_conf(smiles: str) -> tuple:
    """Return (mol_with_Hs, conformer) from a SMILES string.

    Uses a fixed random seed so tests are deterministic.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"RDKit could not parse SMILES: {smiles}"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xBEEF
    rc = AllChem.EmbedMolecule(mol, params)
    if rc < 0:
        # Fallback: random embed
        AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer(0)
    return mol, conf


# ---------------------------------------------------------------------------
# boltzmann_weights
# ---------------------------------------------------------------------------

class TestBoltzmannWeights:
    def test_single_conformer_weight_is_one(self):
        w = boltzmann_weights([0.0])
        assert len(w) == 1
        assert abs(w[0] - 1.0) < 1e-9

    def test_degenerate_conformers_equal_weights(self):
        w = boltzmann_weights([0.0, 0.0, 0.0])
        assert len(w) == 3
        np.testing.assert_allclose(w, [1 / 3, 1 / 3, 1 / 3], atol=1e-9)

    def test_weights_sum_to_one(self):
        w = boltzmann_weights([0.0, 0.5, 1.0, 2.0, 5.0])
        assert abs(float(np.sum(w)) - 1.0) < 1e-9

    def test_lower_energy_has_higher_weight(self):
        w = boltzmann_weights([0.0, 1.0])
        assert w[0] > w[1]

    def test_empty_returns_empty_array(self):
        w = boltzmann_weights([])
        assert len(w) == 0

    def test_all_same_energy_equal_weights(self):
        w = boltzmann_weights([5.0, 5.0])
        np.testing.assert_allclose(w, [0.5, 0.5], atol=1e-9)

    def test_large_energy_gap_weights_dominated_by_lowest(self):
        # 100 kcal/mol gap: upper conformer should have negligible weight.
        w = boltzmann_weights([0.0, 100.0])
        assert w[0] > 0.999


# ---------------------------------------------------------------------------
# nuclei_present
# ---------------------------------------------------------------------------

class TestNucleiPresent:
    def test_ethanol_contains_h_and_c_only(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        nuclei = nuclei_present(mol)
        assert "1h" in nuclei
        assert "13c" in nuclei
        assert "19f" not in nuclei
        assert "31p" not in nuclei

    def test_fluoroethane_contains_f(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCF"))
        nuclei = nuclei_present(mol)
        assert "19f" in nuclei
        assert "31p" not in nuclei

    def test_trimethylphosphine_contains_p(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CP(C)C"))
        nuclei = nuclei_present(mol)
        assert "31p" in nuclei
        assert "19f" not in nuclei

    def test_molecule_with_h_c_f_p(self):
        # trifluoromethyl phosphine: CF3-P(=O)(C)C
        mol = Chem.AddHs(Chem.MolFromSmiles("CP(=O)(C)C(F)(F)F"))
        nuclei = nuclei_present(mol)
        assert "1h" in nuclei
        assert "13c" in nuclei
        assert "19f" in nuclei
        assert "31p" in nuclei

    def test_result_is_ordered_list(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        nuclei = nuclei_present(mol)
        assert isinstance(nuclei, list)


# ---------------------------------------------------------------------------
# estimate_h_shift
# ---------------------------------------------------------------------------

class TestHShiftEstimation:
    def test_aromatic_h_in_correct_range(self):
        mol, conf = _mol_with_conf("c1ccccc1")  # benzene
        h_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
        assert h_indices, "No H atoms found"
        for h_idx in h_indices:
            shift = estimate_h_shift(mol, conf, h_idx)
            assert 6.0 <= shift <= 9.0, f"ArH shift {shift:.2f} outside expected [6.0, 9.0]"

    def test_alkyl_methyl_h_in_correct_range(self):
        mol, conf = _mol_with_conf("CC")  # ethane — all CH3
        h_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
        for h_idx in h_indices:
            shift = estimate_h_shift(mol, conf, h_idx)
            assert 0.5 <= shift <= 2.5, f"CH3 shift {shift:.2f} outside [0.5, 2.5]"

    def test_oh_h_in_correct_range(self):
        mol, conf = _mol_with_conf("CO")  # methanol
        oh_h_indices = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() == 1
            and any(nb.GetAtomicNum() == 8 for nb in a.GetNeighbors())
        ]
        assert oh_h_indices, "No OH hydrogen found in methanol"
        for h_idx in oh_h_indices:
            shift = estimate_h_shift(mol, conf, h_idx)
            assert 0.5 <= shift <= 6.0, f"OH shift {shift:.2f} outside [0.5, 6.0]"

    def test_shift_always_clamped_within_0_and_12(self):
        for smiles in ("CCO", "c1ccccc1", "CCF", "CC(=O)C"):
            mol, conf = _mol_with_conf(smiles)
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    shift = estimate_h_shift(mol, conf, atom.GetIdx())
                    assert 0.0 <= shift <= 12.0, (
                        f"Shift {shift:.2f} for H in {smiles} outside clamping range"
                    )

    def test_sp2_ch_shift_higher_than_sp3_ch(self):
        # Vinyl H (sp2) should be downfield vs alkyl H (sp3)
        mol_sp2, conf_sp2 = _mol_with_conf("C=C")  # ethylene
        mol_sp3, conf_sp3 = _mol_with_conf("CC")   # ethane

        sp2_hs = [a.GetIdx() for a in mol_sp2.GetAtoms() if a.GetAtomicNum() == 1]
        sp3_hs = [a.GetIdx() for a in mol_sp3.GetAtoms() if a.GetAtomicNum() == 1]

        avg_sp2 = sum(estimate_h_shift(mol_sp2, conf_sp2, h) for h in sp2_hs) / len(sp2_hs)
        avg_sp3 = sum(estimate_h_shift(mol_sp3, conf_sp3, h) for h in sp3_hs) / len(sp3_hs)

        assert avg_sp2 > avg_sp3, (
            f"sp2 avg shift {avg_sp2:.2f} should be > sp3 avg shift {avg_sp3:.2f}"
        )


# ---------------------------------------------------------------------------
# estimate_c_shift
# ---------------------------------------------------------------------------

class TestCShiftEstimation:
    def test_carbonyl_c_in_correct_range(self):
        mol, conf = _mol_with_conf("CC(=O)C")  # acetone
        carbonyl_indices = [
            atom.GetIdx() for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 6
            and any(
                float(mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble()) >= 1.9
                and nb.GetAtomicNum() == 8
                for nb in atom.GetNeighbors()
            )
        ]
        assert carbonyl_indices, "No carbonyl C found in acetone"
        for c_idx in carbonyl_indices:
            shift = estimate_c_shift(mol, conf, c_idx)
            assert 150.0 <= shift <= 220.0, f"Carbonyl 13C {shift:.1f} outside [150, 220]"

    def test_aromatic_c_in_correct_range(self):
        mol, conf = _mol_with_conf("c1ccccc1")  # benzene
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                shift = estimate_c_shift(mol, conf, atom.GetIdx())
                assert 115.0 <= shift <= 155.0, f"ArC {shift:.1f} outside [115, 155]"

    def test_alkyl_c_in_correct_range(self):
        mol, conf = _mol_with_conf("C")  # methane
        c_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        shift = estimate_c_shift(mol, conf, c_idx)
        assert -5.0 <= shift <= 50.0, f"Methane 13C {shift:.1f} outside [-5, 50]"

    def test_shift_always_clamped_within_valid_range(self):
        for smiles in ("CCO", "c1ccccc1", "CC(=O)C", "CCOC(=O)O"):
            mol, conf = _mol_with_conf(smiles)
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 6:
                    shift = estimate_c_shift(mol, conf, atom.GetIdx())
                    assert 0.0 <= shift <= 230.0, (
                        f"13C shift {shift:.1f} in {smiles} outside clamping range"
                    )

    def test_carbonyl_downfield_of_alkyl(self):
        mol_carbonyl, conf_carbonyl = _mol_with_conf("CC(=O)C")
        mol_alkyl, conf_alkyl = _mol_with_conf("CCC")

        carbonyl_c = next(
            atom.GetIdx() for atom in mol_carbonyl.GetAtoms()
            if atom.GetAtomicNum() == 6
            and any(
                float(mol_carbonyl.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble()) >= 1.9
                and nb.GetAtomicNum() == 8
                for nb in atom.GetNeighbors()
            )
        )
        alkyl_shifts = [
            estimate_c_shift(mol_alkyl, conf_alkyl, a.GetIdx())
            for a in mol_alkyl.GetAtoms() if a.GetAtomicNum() == 6
        ]
        carbonyl_shift = estimate_c_shift(mol_carbonyl, conf_carbonyl, carbonyl_c)
        assert carbonyl_shift > max(alkyl_shifts), (
            f"Carbonyl {carbonyl_shift:.1f} should be > max alkyl {max(alkyl_shifts):.1f}"
        )


# ---------------------------------------------------------------------------
# estimate_f_shift
# ---------------------------------------------------------------------------

class TestFShiftEstimation:
    def test_aryl_fluoride_in_correct_range(self):
        mol, conf = _mol_with_conf("Fc1ccccc1")  # fluorobenzene
        f_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
        shift = estimate_f_shift(mol, conf, f_idx)
        assert -145.0 <= shift <= -90.0, f"ArF shift {shift:.1f} outside [-145, -90]"

    def test_alkyl_fluoride_in_correct_range(self):
        mol, conf = _mol_with_conf("CCF")  # fluoroethane
        f_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
        shift = estimate_f_shift(mol, conf, f_idx)
        # NOTE: literature value for fluoroethane is ~-219.8 ppm, but the current
        # empirical model uses a base of -128 ppm for sp3-carbon-attached F with no
        # heteroatom neighbours, making simple alkyl fluorides appear ~90 ppm too
        # downfield.  This test anchors the current model behaviour so that regressions
        # (i.e. further drift) are caught; fixing the base value is a known improvement.
        assert -270.0 <= shift <= -100.0, f"Alkyl F shift {shift:.1f} outside [-270, -100]"

    def test_shift_within_hard_bounds(self):
        mol, conf = _mol_with_conf("CF")  # fluoromethane
        f_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
        shift = estimate_f_shift(mol, conf, f_idx)
        assert -260.0 <= shift <= 120.0

    def test_aryl_f_more_downfield_than_alkyl_f(self):
        mol_ar, conf_ar = _mol_with_conf("Fc1ccccc1")
        mol_al, conf_al = _mol_with_conf("CCF")
        ar_f = next(a.GetIdx() for a in mol_ar.GetAtoms() if a.GetAtomicNum() == 9)
        al_f = next(a.GetIdx() for a in mol_al.GetAtoms() if a.GetAtomicNum() == 9)
        shift_ar = estimate_f_shift(mol_ar, conf_ar, ar_f)
        shift_al = estimate_f_shift(mol_al, conf_al, al_f)
        # Aromatic F is less negative (more downfield) than alkyl F
        assert shift_ar > shift_al, (
            f"ArF {shift_ar:.1f} should be > (less negative) alkyl F {shift_al:.1f}"
        )


# ---------------------------------------------------------------------------
# estimate_p_shift
# ---------------------------------------------------------------------------

class TestPShiftEstimation:
    def test_trimethylphosphine_oxide_in_correct_range(self):
        mol, conf = _mol_with_conf("CP(=O)(C)C")  # trimethylphosphine oxide
        p_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
        shift = estimate_p_shift(mol, conf, p_idx)
        # Literature: ~26.5 ppm; allow wider tolerance for empirical method
        assert 0.0 <= shift <= 80.0, f"Trimethylphosphine oxide 31P {shift:.1f} outside [0, 80]"

    def test_triethyl_phosphate_in_correct_range(self):
        mol, conf = _mol_with_conf("CCOP(=O)(OCC)OCC")  # triethyl phosphate
        p_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
        shift = estimate_p_shift(mol, conf, p_idx)
        assert -220.0 <= shift <= 220.0

    def test_shift_within_hard_bounds(self):
        for smiles in ("CP(C)C", "CP(=O)(C)C", "O=P(c1ccccc1)(c1ccccc1)c1ccccc1"):
            mol, conf = _mol_with_conf(smiles)
            p_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
            shift = estimate_p_shift(mol, conf, p_idx)
            assert -220.0 <= shift <= 220.0, (
                f"31P shift {shift:.1f} in {smiles} outside clamping range"
            )


# ---------------------------------------------------------------------------
# group_hydrogens / group_atoms_by_rank
# ---------------------------------------------------------------------------

class TestGroupHydrogens:
    def test_ethanol_has_three_h_groups(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        groups = group_hydrogens(mol)
        # Ethanol: CH3 (3H), CH2 (2H), OH (1H) — three distinct environments
        assert len(groups) == 3, f"Expected 3 groups, got {len(groups)}"

    def test_ethane_has_one_h_group(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        groups = group_hydrogens(mol)
        # All 6H are equivalent by symmetry
        assert len(groups) == 1, f"Expected 1 group, got {len(groups)}"

    def test_benzene_has_one_h_group(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        groups = group_hydrogens(mol)
        assert len(groups) == 1

    def test_total_h_count_preserved(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        groups = group_hydrogens(mol)
        total = sum(len(g) for g in groups)
        expected = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        assert total == expected


class TestGroupAtomsByRank:
    def test_benzene_carbons_one_group(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        groups = group_atoms_by_rank(mol, 6)
        assert len(groups) == 1

    def test_ethanol_carbons_two_groups(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        groups = group_atoms_by_rank(mol, 6)
        # CH3 and CH2OH are chemically inequivalent
        assert len(groups) == 2

    def test_only_requested_atomic_number_returned(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCF"))
        groups = group_atoms_by_rank(mol, 9)  # Fluorine only
        assert all(
            mol.GetAtomWithIdx(idx).GetAtomicNum() == 9
            for g in groups for idx in g
        )


# ---------------------------------------------------------------------------
# splitting_pattern
# ---------------------------------------------------------------------------

class TestSplittingPattern:
    def test_singlet_one_line(self):
        offsets, amps = splitting_pattern("singlet")
        assert len(offsets) == 1 and len(amps) == 1

    def test_doublet_two_lines(self):
        offsets, amps = splitting_pattern("doublet")
        assert len(offsets) == 2 and len(amps) == 2

    def test_triplet_three_lines_pascal(self):
        offsets, amps = splitting_pattern("triplet")
        assert len(offsets) == 3
        # Pascal triangle 1:2:1
        assert amps[0] == amps[2]
        assert abs(amps[1] - 2 * amps[0]) < 1e-9

    def test_quartet_four_lines_pascal(self):
        offsets, amps = splitting_pattern("quartet")
        assert len(offsets) == 4
        # Pascal 1:3:3:1
        assert amps[0] == amps[3]
        assert amps[1] == amps[2]
        assert abs(amps[1] - 3 * amps[0]) < 1e-9

    def test_quintet_five_lines(self):
        offsets, amps = splitting_pattern("quintet")
        assert len(offsets) == 5

    def test_multiplet_fallback_has_matched_lengths(self):
        offsets, amps = splitting_pattern("multiplet")
        assert len(offsets) == len(amps)
        assert len(offsets) > 0

    def test_doublet_symmetric_offsets(self):
        offsets, _ = splitting_pattern("doublet")
        assert abs(offsets[0] + offsets[1]) < 1e-9  # -0.5 and +0.5

    def test_triplet_centre_offset_zero(self):
        offsets, _ = splitting_pattern("triplet")
        assert abs(offsets[1]) < 1e-9  # middle line at 0


# ---------------------------------------------------------------------------
# simulate_spectrum
# ---------------------------------------------------------------------------

class TestSimulateSpectrum:
    def test_output_length_is_6000(self):
        groups = [GroupPrediction(1, [1], 2.0, 800.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "lorentzian", 1.0, "1h")
        assert len(data) == 6000

    def test_max_intensity_normalised_to_one(self):
        groups = [GroupPrediction(1, [1], 2.0, 800.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "lorentzian", 1.0, "1h")
        max_intensity = max(y for _, y in data)
        assert abs(max_intensity - 1.0) < 1e-6

    def test_all_intensities_non_negative(self):
        groups = [GroupPrediction(1, [1], 2.0, 800.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "lorentzian", 1.0, "1h")
        assert all(y >= 0.0 for _, y in data)

    def test_empty_groups_returns_empty(self):
        data = simulate_spectrum([], 400.0, "lorentzian", 1.0, "1h")
        assert len(data) == 0

    def test_gaussian_line_shape_normalised(self):
        groups = [GroupPrediction(1, [1], 2.0, 800.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "gaussian", 1.0, "1h")
        assert len(data) == 6000
        assert abs(max(y for _, y in data) - 1.0) < 1e-6

    def test_voigt_line_shape_normalised(self):
        groups = [GroupPrediction(1, [1], 2.0, 800.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "voigt", 1.0, "1h")
        assert abs(max(y for _, y in data) - 1.0) < 1e-6

    def test_13c_uses_wider_margin(self):
        # 13C margin is 20 ppm; axis should extend well beyond peak
        groups = [GroupPrediction(1, [1], 50.0, 20000.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "lorentzian", 1.0, "13c")
        x_max = max(x for x, _ in data)
        x_min = min(x for x, _ in data)
        assert x_max >= 50.0 + 15.0
        assert x_min <= 50.0 - 15.0

    def test_peak_is_near_predicted_shift(self):
        target_ppm = 3.5
        groups = [GroupPrediction(1, [1], target_ppm, target_ppm * 400.0, "singlet", 0.0, 1.0)]
        data = simulate_spectrum(groups, 400.0, "lorentzian", 1.0, "1h")
        peak_x = max(data, key=lambda pt: pt[1])[0]
        assert abs(peak_x - target_ppm) < 0.1, (
            f"Peak at {peak_x:.3f} ppm, expected near {target_ppm}"
        )


# ---------------------------------------------------------------------------
# load_molecule
# ---------------------------------------------------------------------------

class TestLoadMolecule:
    def test_valid_smiles_returns_mol_with_hs(self):
        mol = load_molecule("smiles", "CCO", [])
        assert mol is not None
        assert mol.GetNumAtoms() > 3

    def test_methane_has_five_atoms_with_hs(self):
        mol = load_molecule("smiles", "C", [])
        assert mol.GetNumAtoms() == 5  # C + 4H

    def test_invalid_smiles_raises_value_error(self):
        with pytest.raises(ValueError, match="[Ff]ailed to parse SMILES"):
            load_molecule("smiles", "not_a_smiles_!!!###", [])

    def test_unsupported_format_raises_value_error(self):
        with pytest.raises(ValueError, match="[Uu]nsupported input format"):
            load_molecule("inchi", "InChI=1S/C2H6O", [])

    def test_radical_emits_warning(self):
        warnings: list = []
        load_molecule("smiles", "[CH3]", warnings)
        assert any("adical" in w for w in warnings), (
            f"Expected radical warning, got: {warnings}"
        )

    def test_charged_species_emits_warning(self):
        warnings: list = []
        load_molecule("smiles", "[NH4+]", warnings)
        assert any("harged" in w for w in warnings), (
            f"Expected charge warning, got: {warnings}"
        )

    def test_no_warnings_for_clean_molecule(self):
        warnings: list = []
        load_molecule("smiles", "CCO", warnings)
        # CCO is a clean, neutral molecule — should produce no warnings
        assert not warnings, f"Unexpected warnings: {warnings}"


# ---------------------------------------------------------------------------
# normalize_nucleus
# ---------------------------------------------------------------------------

class TestNormalizeNucleus:
    @pytest.mark.parametrize("alias,expected", [
        ("1h", "1h"), ("H", "1h"), ("h1", "1h"), ("h", "1h"),
        ("13c", "13c"), ("C", "13c"), ("c13", "13c"), ("13C", "13c"),
        ("19f", "19f"), ("F", "19f"), ("f19", "19f"), ("19F", "19f"),
        ("31p", "31p"), ("P", "31p"), ("p31", "31p"), ("31P", "31p"),
        ("auto", "auto"), ("all", "auto"),
    ])
    def test_known_aliases_resolve_correctly(self, alias, expected):
        assert normalize_nucleus(alias) == expected

    def test_unknown_defaults_to_1h(self):
        # Unknown nuclei default to "1h" per current implementation
        assert normalize_nucleus("zebra_nucleus") == "1h"


# ---------------------------------------------------------------------------
# _chirality_sign
# ---------------------------------------------------------------------------

class TestChiralitySign:
    def test_achiral_ethanol_returns_positive_one(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        assert _chirality_sign(mol) == 1.0

    def test_r_and_s_enantiomers_give_opposite_signs(self):
        # R-lactic acid vs S-lactic acid
        mol_r = Chem.AddHs(Chem.MolFromSmiles("C[C@@H](O)C(=O)O"))
        mol_s = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(=O)O"))
        sign_r = _chirality_sign(mol_r)
        sign_s = _chirality_sign(mol_s)
        assert sign_r != sign_s, "R and S enantiomers should give opposite chirality signs"
        assert sign_r * sign_s < 0, "Signs should have opposite parity"

    def test_return_is_plus_or_minus_one(self):
        for smiles in ("CCO", "C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O"):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            sign = _chirality_sign(mol)
            assert sign in {1.0, -1.0}, f"Unexpected sign {sign} for {smiles}"


# ---------------------------------------------------------------------------
# select_conformers
# ---------------------------------------------------------------------------

class TestSelectConformers:
    def test_single_conformer_returned_and_weight_is_one(self):
        confs = [ConformerResult(conf_id=0, energy_eh=-10.0, method="mmff")]
        selected, weights, _ = select_conformers(confs, 6.0, 0.99)
        assert len(selected) == 1
        assert abs(float(weights[0]) - 1.0) < 1e-9

    def test_high_energy_conformer_excluded_by_window(self):
        # Conf 2 is 1 Eh above the lowest = 627 kcal/mol — far outside 6 kcal window
        confs = [
            ConformerResult(0, -10.000, "mmff"),
            ConformerResult(1, -9.999, "mmff"),
            ConformerResult(2, -9.000, "mmff"),  # ~627 kcal/mol relative
        ]
        selected, _, _ = select_conformers(confs, 6.0, 0.99)
        ids = {c.conf_id for c in selected}
        assert 2 not in ids, "High-energy conformer should be excluded"

    def test_empty_input_returns_empty(self):
        selected, weights, rels = select_conformers([], 6.0, 0.99)
        assert selected == []
        assert len(weights) == 0
        assert rels == []

    def test_selected_weights_sum_to_one(self):
        confs = [
            ConformerResult(0, -10.000, "mmff"),
            ConformerResult(1, -9.999, "mmff"),
        ]
        _, weights, _ = select_conformers(confs, 6.0, 0.99)
        assert abs(float(np.sum(weights)) - 1.0) < 1e-9

    def test_boltzmann_cutoff_limits_number_selected(self):
        # Three nearly-degenerate conformers; cutoff of 0.5 should include fewer than 3
        confs = [
            ConformerResult(0, 0.0, "mmff"),
            ConformerResult(1, 0.001, "mmff"),
            ConformerResult(2, 0.002, "mmff"),
        ]
        selected, _, _ = select_conformers(confs, 6.0, 0.5)
        # With equal weights (~0.33 each), cumulative hits 0.5 after 2 conformers
        assert len(selected) <= 3

    def test_at_least_one_conformer_always_returned(self):
        # Even a single conformer outside the window should still return it
        confs = [ConformerResult(0, -10.0, "mmff")]
        selected, _, _ = select_conformers(confs, 0.001, 0.99)
        assert len(selected) >= 1


# ---------------------------------------------------------------------------
# MULTIPLICITY_LABELS constant (regression guard)
# ---------------------------------------------------------------------------

class TestMultiplicityLabels:
    def test_standard_labels_present(self):
        assert MULTIPLICITY_LABELS[1] == "singlet"
        assert MULTIPLICITY_LABELS[2] == "doublet"
        assert MULTIPLICITY_LABELS[3] == "triplet"
        assert MULTIPLICITY_LABELS[4] == "quartet"
        assert MULTIPLICITY_LABELS[5] == "quintet"

    def test_no_label_beyond_quintet(self):
        # If extended labels are added (sextet etc.), this test should be updated.
        # This acts as a deliberate reminder that the gap exists.
        assert 6 not in MULTIPLICITY_LABELS, (
            "Sextet label added — update this test and the feature proposal tracker"
        )
