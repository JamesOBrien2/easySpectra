#include "core/spectrum.h"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: easynmr-expcheck <experimental-spectrum-file>\n";
        return 1;
    }
    const std::string path = argv[1];
    const auto loaded = easynmr::load_experimental_spectrum(path);
    if (!loaded.error_message.empty()) {
        std::cerr << "error: " << loaded.error_message << "\n";
        if (!loaded.detected_format.empty()) {
            std::cerr << "detected_format: " << loaded.detected_format << "\n";
        }
        return 2;
    }
    std::cout << "ok\n"
              << "format: " << loaded.detected_format << "\n"
              << "points: " << loaded.points.size() << "\n";
    if (!loaded.points.empty()) {
        std::cout << "x_max: " << loaded.points.front().ppm << "\n"
                  << "x_min: " << loaded.points.back().ppm << "\n";
    }
    return 0;
}
