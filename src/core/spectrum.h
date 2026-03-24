#pragma once

#include <string>
#include <vector>

namespace easynmr {

struct SpectrumPoint {
    double ppm = 0.0;
    double intensity = 0.0;
};

std::vector<SpectrumPoint> load_spectrum_csv(const std::string &path);

} // namespace easynmr
