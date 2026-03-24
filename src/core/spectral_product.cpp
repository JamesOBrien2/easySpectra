#include "core/spectral_product.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <utility>

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

std::vector<std::string> parse_csv_line(const std::string &line) {
    std::vector<std::string> fields;
    std::string cur;
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
            continue;
        }
        if (c == ',' && !in_quotes) {
            fields.push_back(cur);
            cur.clear();
            continue;
        }
        cur.push_back(c);
    }
    fields.push_back(cur);
    return fields;
}

} // namespace

std::vector<SpectralProductFiles> load_spectral_products_manifest(const std::string &path) {
    std::vector<SpectralProductFiles> out;
    std::ifstream input(path);
    if (!input) {
        return out;
    }

    std::string line;
    bool header = true;
    while (std::getline(input, line)) {
        if (header) {
            header = false;
            continue;
        }
        if (line.empty()) {
            continue;
        }

        const auto fields = parse_csv_line(line);
        if (fields.size() < 4) {
            continue;
        }

        SpectralProductFiles files;
        files.label = trim_copy(fields[0]);
        files.spectrum_csv = trim_copy(fields[1]);
        files.peaks_csv = trim_copy(fields[2]);
        files.assignments_csv = trim_copy(fields[3]);

        if (!files.label.empty() && !files.spectrum_csv.empty()) {
            out.push_back(std::move(files));
        }
    }

    return out;
}

} // namespace easynmr
