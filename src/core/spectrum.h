#pragma once

#include <string>
#include <vector>

namespace easynmr {

struct SpectrumPoint {
    double ppm = 0.0;
    double intensity = 0.0;
};

struct ExperimentalSpectrumLoadResult {
    std::vector<SpectrumPoint> points;
    std::string detected_format;
    std::string error_message;
};

std::vector<SpectrumPoint> load_spectrum_csv(const std::string &path);
ExperimentalSpectrumLoadResult load_experimental_spectrum(const std::string &path);

} // namespace easynmr
