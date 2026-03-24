#pragma once

#include <string>
#include <vector>

namespace easynmr {

struct SpectralProductFiles {
    std::string label;
    std::string spectrum_csv;
    std::string peaks_csv;
    std::string assignments_csv;
};

std::vector<SpectralProductFiles> load_spectral_products_manifest(const std::string &path);

} // namespace easynmr
