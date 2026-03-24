#pragma once

#include <string>
#include <vector>

namespace easynmr {

enum class InputFormat {
    Smiles,
    Mol,
    Sdf,
    Xyz,
    Unknown
};

inline std::string to_string(InputFormat format) {
    switch (format) {
    case InputFormat::Smiles:
        return "smiles";
    case InputFormat::Mol:
        return "mol";
    case InputFormat::Sdf:
        return "sdf";
    case InputFormat::Xyz:
        return "xyz";
    default:
        return "unknown";
    }
}

inline InputFormat input_format_from_string(const std::string &value) {
    if (value == "smiles") {
        return InputFormat::Smiles;
    }
    if (value == "mol") {
        return InputFormat::Mol;
    }
    if (value == "sdf") {
        return InputFormat::Sdf;
    }
    if (value == "xyz") {
        return InputFormat::Xyz;
    }
    return InputFormat::Unknown;
}

struct JobConfig {
    std::string job_name = "untitled";
    std::string input_value;
    InputFormat input_format = InputFormat::Smiles;
    std::string nucleus = "auto";
    std::string solvent = "cdcl3";
    double frequency_mhz = 400.0;
    bool auto_protomer = true;
    bool manual_lock = false;
    double ph = 7.0;
    std::size_t max_conformers = 50;
    double boltzmann_cutoff = 0.99;
    double energy_window_kcal = 6.0;
    std::string line_shape = "lorentzian";
    double fwhm_hz = 1.0;
    bool need_editable_xyz = false;
    std::string output_root = "output";
    std::string python_executable = EASYNMR_DEFAULT_PYTHON;
    std::string backend_script = EASYNMR_BACKEND_SCRIPT;
};

struct JobOutputs {
    std::string job_id;
    std::string status;
    std::string message;
    std::string output_dir;
    std::string spectrum_csv;
    std::string peaks_csv;
    std::string assignments_json;
    std::string assignments_csv;
    std::string structure_svg;
    std::string structure_atoms_csv;
    std::string structure_bonds_csv;
    std::string structure_xyz;
    std::string spectra_manifest_csv;
    std::string audit_json;
    std::vector<std::string> warnings;
};

} // namespace easynmr
