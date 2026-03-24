#include "core/pipeline.h"

#include <chrono>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace easynmr {
namespace {

std::string now_job_id() {
    const auto now = std::chrono::system_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
    return "job-" + std::to_string(micros);
}

std::string escape_json(const std::string &input) {
    std::string out;
    out.reserve(input.size());
    for (char c : input) {
        switch (c) {
        case '\\':
            out += "\\\\";
            break;
        case '"':
            out += "\\\"";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            out += c;
            break;
        }
    }
    return out;
}

std::string shell_quote(const std::string &value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
        if (c == '"' || c == '\\') {
            escaped.push_back('\\');
        }
        escaped.push_back(c);
    }
    return "\"" + escaped + "\"";
}

std::string extract_string_field(const std::string &json, const std::string &key) {
    std::regex re("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (std::regex_search(json, match, re) && match.size() > 1) {
        return match[1].str();
    }
    return {};
}

std::vector<std::string> extract_warnings(const std::string &json) {
    std::vector<std::string> warnings;
    const std::regex block_re("\"warnings\"\\s*:\\s*\\[(.*?)\\]");
    std::smatch block_match;
    if (!std::regex_search(json, block_match, block_re) || block_match.size() <= 1) {
        return warnings;
    }
    const std::string block = block_match[1].str();
    const std::regex item_re("\"([^\"]*)\"");
    for (auto it = std::sregex_iterator(block.begin(), block.end(), item_re); it != std::sregex_iterator(); ++it) {
        warnings.push_back((*it)[1].str());
    }
    return warnings;
}

double extract_number_field(const std::string &json, const std::string &key, double fallback = 0.0) {
    std::regex re("\"" + key + "\"\\s*:\\s*([0-9]+(?:\\.[0-9]+)?)");
    std::smatch match;
    if (std::regex_search(json, match, re) && match.size() > 1) {
        try {
            return std::stod(match[1].str());
        } catch (...) {
            return fallback;
        }
    }
    return fallback;
}

std::string read_file_if_exists(const std::filesystem::path &path) {
    std::ifstream in(path);
    if (!in) {
        return {};
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

PipelineProgress parse_progress_json(const std::string &json) {
    PipelineProgress progress;
    progress.stage = extract_string_field(json, "stage");
    progress.message = extract_string_field(json, "message");
    progress.fraction = extract_number_field(json, "fraction", 0.0);
    return progress;
}

} // namespace

JobOutputs Pipeline::run(const JobConfig &config) const {
    return run_impl(config, "predict", nullptr);
}

JobOutputs Pipeline::run(const JobConfig &config, const std::function<void(const PipelineProgress &)> &progress_cb) const {
    return run_impl(config, "predict", &progress_cb);
}

JobOutputs Pipeline::preview(const JobConfig &config) const {
    return run_impl(config, "preview", nullptr);
}

JobOutputs Pipeline::run_impl(
    const JobConfig &config,
    const std::string &mode,
    const std::function<void(const PipelineProgress &)> *progress_cb) const {
    JobOutputs outputs;
    outputs.job_id = now_job_id();
    outputs.status = "failed";

    const std::filesystem::path output_dir = std::filesystem::path(config.output_root) / outputs.job_id;
    std::filesystem::create_directories(output_dir);

    const std::filesystem::path request_path = output_dir / "request.json";
    const std::filesystem::path response_path = output_dir / "response.json";
    const std::filesystem::path progress_path = output_dir / "progress.json";

    {
        std::ofstream request_file(request_path);
        request_file << "{\n";
        request_file << "  \"job_id\": \"" << escape_json(outputs.job_id) << "\",\n";
        request_file << "  \"mode\": \"" << escape_json(mode) << "\",\n";
        request_file << "  \"job_name\": \"" << escape_json(config.job_name) << "\",\n";
        request_file << "  \"input\": {\n";
        request_file << "    \"format\": \"" << to_string(config.input_format) << "\",\n";
        request_file << "    \"value\": \"" << escape_json(config.input_value) << "\"\n";
        request_file << "  },\n";
        request_file << "  \"settings\": {\n";
        request_file << "    \"frequency_mhz\": " << config.frequency_mhz << ",\n";
        request_file << "    \"solvent\": \"" << escape_json(config.solvent) << "\",\n";
        request_file << "    \"nucleus\": \"" << escape_json(config.nucleus) << "\",\n";
        request_file << "    \"auto_protomer\": " << (config.auto_protomer ? "true" : "false") << ",\n";
        request_file << "    \"manual_lock\": " << (config.manual_lock ? "true" : "false") << ",\n";
        request_file << "    \"ph\": " << config.ph << ",\n";
        request_file << "    \"max_conformers\": " << config.max_conformers << ",\n";
        request_file << "    \"boltzmann_cutoff\": " << config.boltzmann_cutoff << ",\n";
        request_file << "    \"energy_window_kcal\": " << config.energy_window_kcal << ",\n";
        request_file << "    \"line_shape\": \"" << escape_json(config.line_shape) << "\",\n";
        request_file << "    \"fwhm_hz\": " << config.fwhm_hz << "\n";
        request_file << "  },\n";
        request_file << "  \"output_dir\": \"" << escape_json(output_dir.string()) << "\",\n";
        request_file << "  \"progress_json\": \"" << escape_json(progress_path.string()) << "\"\n";
        request_file << "}\n";
    }

    std::ostringstream cmd;
    cmd << shell_quote(config.python_executable) << ' '
        << shell_quote(config.backend_script) << " --request "
        << shell_quote(request_path.string()) << " --response "
        << shell_quote(response_path.string());

    if (progress_cb != nullptr) {
        PipelineProgress p;
        p.stage = "launch";
        p.message = "Launching backend process";
        p.fraction = 0.0;
        (*progress_cb)(p);
    }

    std::atomic<bool> finished = false;
    int rc = 0;
    std::thread backend_thread([&]() {
        rc = std::system(cmd.str().c_str());
        finished = true;
    });

    std::string last_progress_json;
    while (!finished.load()) {
        if (progress_cb != nullptr) {
            const std::string json = read_file_if_exists(progress_path);
            if (!json.empty() && json != last_progress_json) {
                last_progress_json = json;
                (*progress_cb)(parse_progress_json(json));
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    backend_thread.join();

    if (progress_cb != nullptr) {
        const std::string json = read_file_if_exists(progress_path);
        if (!json.empty() && json != last_progress_json) {
            (*progress_cb)(parse_progress_json(json));
        }
    }

    if (rc != 0) {
        outputs.message = "Backend execution failed with exit code " + std::to_string(rc);
        outputs.output_dir = output_dir.string();
        return outputs;
    }

    std::ifstream response_file(response_path);
    if (!response_file) {
        outputs.message = "Backend did not produce response.json";
        outputs.output_dir = output_dir.string();
        return outputs;
    }

    std::stringstream buffer;
    buffer << response_file.rdbuf();
    const std::string json = buffer.str();

    outputs.status = extract_string_field(json, "status");
    outputs.message = extract_string_field(json, "message");
    outputs.output_dir = extract_string_field(json, "output_dir");
    outputs.spectrum_csv = extract_string_field(json, "spectrum_csv");
    outputs.peaks_csv = extract_string_field(json, "peaks_csv");
    outputs.assignments_json = extract_string_field(json, "assignments_json");
    outputs.assignments_csv = extract_string_field(json, "assignments_csv");
    outputs.structure_svg = extract_string_field(json, "structure_svg");
    outputs.structure_atoms_csv = extract_string_field(json, "structure_atoms_csv");
    outputs.structure_bonds_csv = extract_string_field(json, "structure_bonds_csv");
    outputs.structure_xyz = extract_string_field(json, "structure_xyz");
    outputs.spectra_manifest_csv = extract_string_field(json, "spectra_manifest_csv");
    outputs.audit_json = extract_string_field(json, "audit_json");
    outputs.warnings = extract_warnings(json);

    if (outputs.status.empty()) {
        outputs.status = "failed";
    }
    if (outputs.output_dir.empty()) {
        outputs.output_dir = output_dir.string();
    }

    return outputs;
}

} // namespace easynmr
