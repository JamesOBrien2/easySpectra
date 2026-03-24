#include "core/spectrum.h"

#include <cctype>
#include <fstream>
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

} // namespace easynmr
