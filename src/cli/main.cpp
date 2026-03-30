#include "core/batch_parser.h"
#include "core/job.h"
#include "core/pipeline.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void print_help() {
    std::cout << "easySpectra CLI\n"
              << "Usage:\n"
              << "  easynmr --input <value> [options]\n"
              << "  easynmr --input <reactant> --compare <product> [options]\n"
              << "  easynmr --batch-smi <file> [options]\n\n"
              << "Options:\n"
              << "  --input <value>                Input content (e.g. SMILES)\n"
              << "  --input-format <smiles|mol|sdf|xyz>\n"
              << "  --compare <value>              Second structure for reaction comparison\n"
              << "                                 (implies --workflow compare)\n"
              << "  --compare-format <smiles|mol|sdf|xyz>\n"
              << "  --name <job_name>\n"
              << "  --workflow <all|nmr|cd|ir|uvvis|compare|properties>\n"
              << "  --level <low|medium|high>      Level of theory (default: medium)\n"
              << "  --output-dir <path>\n"
              << "  --solvent <cdcl3|dmso|h2o>\n"
              << "  --nucleus <auto|1H|13C|19F|31P>\n"
              << "  --frequency-mhz <number>\n"
              << "  --line-shape <lorentzian|gaussian|voigt>\n"
              << "  --fwhm-hz <number>\n"
              << "  --max-conformers <n>\n"
              << "  --energy-window-kcal <number>\n"
              << "  --boltzmann-cutoff <0-1>\n"
              << "  --ph <number>\n"
              << "  --manual-lock\n"
              << "  --python <path_to_python>\n"
              << "  --batch-sdf <file>             Queue all molecules from an SDF file\n"
              << "  --batch-csv <file>             Queue molecules from a CSV file (smiles + name columns)\n"
              << "  --batch-smi <file>             Queue molecules from a SMILES file (one per line)\n"
              << "  --batch-xyz <file>             Queue molecules from a multi-XYZ file\n"
              << "  --batch-cdxml <file>           Queue molecules from a ChemDraw XML file (requires RDKit)\n"
              << "  -h, --help\n";
}

} // namespace

int main(int argc, char **argv) {
    easynmr::JobConfig config;
    std::vector<easynmr::BatchEntry> batch_entries;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        auto take_value = [&](const std::string &flag) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << flag << "\n";
                std::exit(1);
            }
            ++i;
            return argv[i];
        };

        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
        if (arg == "--input") {
            config.input_value = take_value(arg);
            continue;
        }
        if (arg == "--input-format") {
            config.input_format = easynmr::input_format_from_string(take_value(arg));
            continue;
        }
        if (arg == "--name") {
            config.job_name = take_value(arg);
            continue;
        }
        if (arg == "--workflow") {
            config.workflow_kind = easynmr::workflow_kind_from_string(take_value(arg));
            if (config.workflow_kind == easynmr::WorkflowKind::Unknown) {
                std::cerr << "Unsupported --workflow value (supported: all, nmr, cd, compare)\n";
                return 1;
            }
            continue;
        }
        if (arg == "--compare") {
            config.compare_input_value = take_value(arg);
            // Implicitly enable compare workflow unless user specified another kind.
            if (config.workflow_kind == easynmr::WorkflowKind::All) {
                config.workflow_kind = easynmr::WorkflowKind::Compare;
            }
            continue;
        }
        if (arg == "--compare-format") {
            config.compare_input_format = easynmr::input_format_from_string(take_value(arg));
            continue;
        }
        if (arg == "--output-dir") {
            config.output_root = take_value(arg);
            continue;
        }
        if (arg == "--solvent") {
            config.solvent = take_value(arg);
            continue;
        }
        if (arg == "--nucleus") {
            config.nucleus = take_value(arg);
            continue;
        }
        if (arg == "--frequency-mhz") {
            config.frequency_mhz = std::stod(take_value(arg));
            continue;
        }
        if (arg == "--line-shape") {
            config.line_shape = take_value(arg);
            continue;
        }
        if (arg == "--fwhm-hz") {
            config.fwhm_hz = std::stod(take_value(arg));
            continue;
        }
        if (arg == "--max-conformers") {
            config.max_conformers = static_cast<std::size_t>(std::stoul(take_value(arg)));
            continue;
        }
        if (arg == "--energy-window-kcal") {
            config.energy_window_kcal = std::stod(take_value(arg));
            continue;
        }
        if (arg == "--boltzmann-cutoff") {
            config.boltzmann_cutoff = std::stod(take_value(arg));
            continue;
        }
        if (arg == "--ph") {
            config.ph = std::stod(take_value(arg));
            continue;
        }
        if (arg == "--manual-lock") {
            config.manual_lock = true;
            config.auto_protomer = false;
            continue;
        }
        if (arg == "--python") {
            config.python_executable = take_value(arg);
            continue;
        }
        if (arg == "--level") {
            config.level_of_theory = take_value(arg);
            continue;
        }
        if (arg == "--batch-sdf" || arg == "--batch-csv" || arg == "--batch-smi" ||
            arg == "--batch-xyz" || arg == "--batch-cdxml") {
            auto new_entries = easynmr::parse_batch_file(take_value(arg), config.python_executable);
            for (auto &e : new_entries) batch_entries.push_back(std::move(e));
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        print_help();
        return 1;
    }

    // Batch mode: if any --batch-* files were loaded, run them all.
    if (!batch_entries.empty()) {
        int ok = 0;
        int fail = 0;
        for (const auto &entry : batch_entries) {
            easynmr::JobConfig bcfg = config;
            bcfg.input_value  = entry.input_value;
            bcfg.input_format = entry.input_format;
            if (!entry.name.empty()) bcfg.job_name = entry.name;
            easynmr::Pipeline bpipeline;
            const auto bres = bpipeline.run(bcfg);
            if (bres.status == "ok" || bres.status == "success") {
                std::cout << "[OK]   " << bcfg.job_name << " → " << bres.output_dir << "\n";
                ++ok;
            } else {
                std::cerr << "[FAIL] " << bcfg.job_name << ": " << bres.message << "\n";
                ++fail;
            }
        }
        std::cout << "Batch complete: " << ok << " ok, " << fail << " failed\n";
        return (fail > 0) ? 1 : 0;
    }

    if (config.input_value.empty()) {
        std::cerr << "Missing required --input\n";
        print_help();
        return 1;
    }

    if (config.input_format == easynmr::InputFormat::Unknown) {
        std::cerr << "Unsupported --input-format\n";
        return 1;
    }

    easynmr::Pipeline pipeline;
    const auto result = pipeline.run(config);

    std::cout << "Job: " << result.job_id << "\n"
              << "Status: " << result.status << "\n"
              << "Message: " << result.message << "\n"
              << "Output dir: " << result.output_dir << "\n"
              << "Spectrum CSV: " << result.spectrum_csv << "\n"
              << "Peaks CSV: " << result.peaks_csv << "\n"
              << "Assignments: " << result.assignments_json << "\n"
              << "Assignments CSV: " << result.assignments_csv << "\n"
              << "Spectra manifest CSV: " << result.spectra_manifest_csv << "\n"
              << "Structure SVG: " << result.structure_svg << "\n"
              << "Structure atoms CSV: " << result.structure_atoms_csv << "\n"
              << "Structure bonds CSV: " << result.structure_bonds_csv << "\n"
              << "Audit JSON: " << result.audit_json << "\n";
    if (!result.reaction_summary_json.empty()) {
        std::cout << "Reaction summary JSON: " << result.reaction_summary_json << "\n"
                  << "Product structure SVG: " << result.structure_product_svg << "\n";
    }

    if (!result.warnings.empty()) {
        std::cout << "Warnings:\n";
        for (const auto &warning : result.warnings) {
            std::cout << "  - " << warning << "\n";
        }
    }

    return result.status == "ok" ? 0 : 2;
}
