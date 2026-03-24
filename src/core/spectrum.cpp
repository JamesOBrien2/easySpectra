#include "core/spectrum.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>

namespace easynmr {
namespace {

std::string trim_copy(const std::string &value) {
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])) != 0) {
        ++begin;
    }
    std::size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1])) != 0) {
        --end;
    }
    return value.substr(begin, end - begin);
}

bool try_parse_double(const std::string &text, double &out) {
    const std::string clean = trim_copy(text);
    if (clean.empty()) {
        return false;
    }
    try {
        std::size_t consumed = 0;
        const double parsed = std::stod(clean, &consumed);
        if (consumed != clean.size()) {
            return false;
        }
        out = parsed;
        return true;
    } catch (...) {
        return false;
    }
}

std::string to_lower_copy(const std::string &value) {
    std::string out = value;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

bool has_bruker_raw_layout(const std::filesystem::path &path) {
    if (!std::filesystem::is_directory(path)) {
        return false;
    }
    std::error_code ec;
    const bool has_acqus = std::filesystem::exists(path / "acqus", ec) && !ec;
    const bool has_fid = std::filesystem::exists(path / "fid", ec) && !ec;
    const bool has_ser = std::filesystem::exists(path / "ser", ec) && !ec;
    const bool has_pdata = std::filesystem::exists(path / "pdata", ec) && !ec;
    if (has_acqus && (has_fid || has_ser || has_pdata)) {
        return true;
    }

    // Some exported Bruker directories are nested one level below the selected path.
    for (const auto &entry : std::filesystem::directory_iterator(path, ec)) {
        if (ec || !entry.is_directory()) {
            continue;
        }
        const std::filesystem::path sub = entry.path();
        const bool sub_has_acqus = std::filesystem::exists(sub / "acqus", ec) && !ec;
        const bool sub_has_fid = std::filesystem::exists(sub / "fid", ec) && !ec;
        const bool sub_has_ser = std::filesystem::exists(sub / "ser", ec) && !ec;
        const bool sub_has_pdata = std::filesystem::exists(sub / "pdata", ec) && !ec;
        if (sub_has_acqus && (sub_has_fid || sub_has_ser || sub_has_pdata)) {
            return true;
        }
    }
    return false;
}

std::string guess_exp_format(const std::filesystem::path &path, const std::string &line_sample) {
    const std::string ext = to_lower_copy(path.extension().string());
    const std::string sample = to_lower_copy(line_sample);
    if (sample.find("mnova") != std::string::npos || sample.find("mestrenova") != std::string::npos
        || ext == ".mnova" || ext == ".mnv") {
        return "mnova_text_export";
    }
    if (sample.find("bruker") != std::string::npos || sample.find("topspin") != std::string::npos
        || sample.find("##title") != std::string::npos || sample.find("jcamp") != std::string::npos) {
        return "bruker_text_export";
    }
    if (ext == ".txt" || ext == ".csv" || ext == ".asc" || ext == ".dat" || ext == ".dx") {
        return "generic_text_export";
    }
    return "generic_text_export";
}

std::vector<double> parse_numeric_tokens(const std::string &line) {
    std::vector<double> nums;
    if (line.empty()) {
        return nums;
    }
    std::string normalized = line;

    // Common MNova/European export uses decimal comma with ';' separators.
    if (normalized.find(';') != std::string::npos && normalized.find('.') == std::string::npos) {
        std::replace(normalized.begin(), normalized.end(), ',', '.');
    }

    static const std::regex number_re(R"([+-]?(?:\d+(?:[.,]\d+)?|\.\d+)(?:[eE][+-]?\d+)?)");
    for (auto it = std::sregex_iterator(normalized.begin(), normalized.end(), number_re);
         it != std::sregex_iterator();
         ++it) {
        std::string token = (*it).str();
        std::replace(token.begin(), token.end(), ',', '.');
        double parsed = 0.0;
        if (try_parse_double(token, parsed)) {
            nums.push_back(parsed);
        }
    }
    return nums;
}

bool is_comment_or_header_line(const std::string &line) {
    const std::string clean = trim_copy(line);
    if (clean.empty()) {
        return true;
    }
    if (clean.rfind("#", 0) == 0 || clean.rfind("//", 0) == 0 || clean.rfind(";", 0) == 0) {
        return true;
    }
    if (clean.rfind("##", 0) == 0) {
        return true;
    }
    return false;
}

} // namespace

std::vector<SpectrumPoint> load_spectrum_csv(const std::string &path) {
    std::vector<SpectrumPoint> points;
    std::ifstream input(path);
    if (!input) {
        return points;
    }

    std::string line;
    bool header_skipped = false;
    while (std::getline(input, line)) {
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string ppm_str;
        std::string intensity_str;
        if (!std::getline(ss, ppm_str, ',')) {
            continue;
        }
        if (!std::getline(ss, intensity_str, ',')) {
            continue;
        }

        SpectrumPoint point;
        if (!try_parse_double(ppm_str, point.ppm)) {
            continue;
        }
        if (!try_parse_double(intensity_str, point.intensity)) {
            continue;
        }
        points.push_back(point);
    }

    return points;
}

ExperimentalSpectrumLoadResult load_experimental_spectrum(const std::string &path_text) {
    ExperimentalSpectrumLoadResult result;
    const std::filesystem::path path(path_text);
    if (path_text.empty()) {
        result.error_message = "No experimental spectrum path provided.";
        return result;
    }

    std::error_code ec;
    if (!std::filesystem::exists(path, ec) || ec) {
        result.error_message = "Experimental spectrum path does not exist.";
        return result;
    }

    if (has_bruker_raw_layout(path)) {
        result.error_message =
            "Bruker raw directories are not parsed directly yet. Export a 2-column ASCII/CSV from TopSpin or MNova.";
        result.detected_format = "bruker_raw_directory";
        return result;
    }

    if (!std::filesystem::is_regular_file(path, ec) || ec) {
        result.error_message = "Experimental spectrum path is not a readable file.";
        return result;
    }

    const std::string ext = to_lower_copy(path.extension().string());
    if (ext == ".mnova" || ext == ".mnv") {
        result.error_message =
            "MNova project files are not parsed directly. Export as text/CSV (x,y) and import that file.";
        result.detected_format = "mnova_project_file";
        return result;
    }

    std::ifstream input(path);
    if (!input) {
        result.error_message = "Could not open experimental spectrum file.";
        return result;
    }

    std::string first_nonempty_line;
    std::string line;
    while (std::getline(input, line)) {
        const std::string clean = trim_copy(line);
        if (!clean.empty()) {
            first_nonempty_line = clean;
            break;
        }
    }
    input.clear();
    input.seekg(0);

    result.detected_format = guess_exp_format(path, first_nonempty_line);

    std::vector<SpectrumPoint> parsed;
    while (std::getline(input, line)) {
        if (is_comment_or_header_line(line)) {
            continue;
        }
        const auto nums = parse_numeric_tokens(line);
        if (nums.size() < 2) {
            continue;
        }
        SpectrumPoint p;
        p.ppm = nums[0];
        p.intensity = nums[1];
        if (!std::isfinite(p.ppm) || !std::isfinite(p.intensity)) {
            continue;
        }
        parsed.push_back(p);
    }

    if (parsed.size() < 8) {
        result.error_message =
            "Experimental spectrum parse found too few numeric points. Expected a 2-column x,y text/CSV export.";
        return result;
    }

    std::sort(parsed.begin(), parsed.end(), [](const SpectrumPoint &a, const SpectrumPoint &b) {
        return a.ppm > b.ppm;
    });
    result.points = std::move(parsed);
    return result;
}

} // namespace easynmr
