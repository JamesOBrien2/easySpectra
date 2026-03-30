#pragma once

#include <map>
#include <string>
#include <vector>

namespace easynmr {

struct FukuiAtom {
    int atom_index = 0;   // 1-based
    std::string element;
    double f_plus  = 0.0;
    double f_minus = 0.0;
    double f_zero  = 0.0;
};

struct PkaGroup {
    int atom_index  = 0;  // 1-based
    double pka_est  = 0.0;
    double pka_low  = 0.0;
    double pka_high = 0.0;
    std::string group_name;
    std::string site_type;  // "acidic" or "basic"
};

struct MsAdduct {
    std::string name;
    double mz   = 0.0;
    int charge  = 1;
    std::string mode;  // "+" or "-"
};

struct IsotopePeak {
    int    mass_offset        = 0;    // 0 = M, 1 = M+1, etc.
    double relative_abundance = 0.0;  // 0–100 (100 = base peak)
};

struct MolecularProperties {
    // RDKit descriptors
    std::string formula;
    double mw            = 0.0;
    double exact_mw      = 0.0;
    double logp          = 0.0;
    double tpsa          = 0.0;
    int    hbd           = 0;
    int    hba           = 0;
    int    rotbonds      = 0;
    int    ar_rings      = 0;
    int    heavy_atoms   = 0;
    int    formal_charge = 0;

    // xTB electronic
    double homo_ev       = 0.0;
    double lumo_ev       = 0.0;
    double gap_ev        = 0.0;
    double dipole_debye  = 0.0;
    bool   has_electronic = false;

    // Fukui indices (per atom, including H)
    std::vector<FukuiAtom> fukui;
    bool has_fukui = false;

    // pKa estimates
    std::vector<PkaGroup> pka_groups;

    // Redox estimates (V)
    double e_ox_nhe  = 0.0;
    double e_red_nhe = 0.0;
    double e_ox_fc   = 0.0;
    double e_red_fc  = 0.0;
    bool   has_redox = false;

    // MS adducts + isotope pattern
    double ms_monoisotopic = 0.0;
    std::vector<MsAdduct> ms_adducts;
    std::vector<IsotopePeak> isotope_peaks;
    bool has_ms = false;

    std::vector<std::string> warnings;
    bool valid = false;
};

MolecularProperties load_properties_json(const std::string &path);

} // namespace easynmr
