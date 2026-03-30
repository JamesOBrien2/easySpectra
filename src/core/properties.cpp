#include "core/properties.h"

#include <fstream>
#include <regex>
#include <sstream>
#include <string>

namespace easynmr {
namespace {

// ── minimal JSON field extractors ────────────────────────────────────────────

std::string extract_string_field(const std::string &json, const std::string &key) {
    std::regex re("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch m;
    if (std::regex_search(json, m, re) && m.size() > 1) {
        return m[1].str();
    }
    return {};
}

double extract_number_field(const std::string &json, const std::string &key, double fallback = 0.0) {
    std::regex re("\"" + key + "\"\\s*:\\s*(-?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)");
    std::smatch m;
    if (std::regex_search(json, m, re) && m.size() > 1) {
        try { return std::stod(m[1].str()); } catch (...) {}
    }
    return fallback;
}

int extract_int_field(const std::string &json, const std::string &key, int fallback = 0) {
    return static_cast<int>(extract_number_field(json, key, static_cast<double>(fallback)));
}

bool extract_bool_field(const std::string &json, const std::string &key, bool fallback = false) {
    std::regex re("\"" + key + "\"\\s*:\\s*(true|false)");
    std::smatch m;
    if (std::regex_search(json, m, re) && m.size() > 1) {
        return m[1].str() == "true";
    }
    return fallback;
}

// Extract the raw text of a JSON array or object value for a given key.
// Handles simple one-level nesting (finds matching bracket).
std::string extract_block(const std::string &json, const std::string &key, char open, char close) {
    const std::string search = "\"" + key + "\"";
    auto pos = json.find(search);
    if (pos == std::string::npos) return {};
    auto bracket = json.find(open, pos + search.size());
    if (bracket == std::string::npos) return {};
    int depth = 0;
    for (std::size_t i = bracket; i < json.size(); ++i) {
        if (json[i] == open) ++depth;
        else if (json[i] == close) {
            --depth;
            if (depth == 0) {
                return json.substr(bracket, i - bracket + 1);
            }
        }
    }
    return {};
}

// Split a JSON array of objects into individual object strings.
std::vector<std::string> split_json_objects(const std::string &array_json) {
    std::vector<std::string> objects;
    if (array_json.empty() || array_json.front() != '[') return objects;
    int depth = 0;
    std::size_t start = std::string::npos;
    for (std::size_t i = 0; i < array_json.size(); ++i) {
        char c = array_json[i];
        if (c == '{') {
            if (depth == 0) start = i;
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0 && start != std::string::npos) {
                objects.push_back(array_json.substr(start, i - start + 1));
                start = std::string::npos;
            }
        }
    }
    return objects;
}

// Extract string array values from a JSON array like ["a","b","c"].
std::vector<std::string> extract_string_array(const std::string &json, const std::string &key) {
    std::vector<std::string> result;
    std::string block = extract_block(json, key, '[', ']');
    if (block.empty()) return result;
    std::regex item_re("\"([^\"]*)\"");
    for (auto it = std::sregex_iterator(block.begin(), block.end(), item_re);
         it != std::sregex_iterator(); ++it) {
        result.push_back((*it)[1].str());
    }
    return result;
}

} // anonymous namespace

MolecularProperties load_properties_json(const std::string &path) {
    MolecularProperties props;
    std::ifstream in(path);
    if (!in) return props;
    std::stringstream buf;
    buf << in.rdbuf();
    const std::string json = buf.str();
    if (json.empty()) return props;

    // RDKit section
    std::string rdkit_block = extract_block(json, "rdkit", '{', '}');
    if (!rdkit_block.empty()) {
        props.formula       = extract_string_field(rdkit_block, "formula");
        props.mw            = extract_number_field(rdkit_block, "mw");
        props.exact_mw      = extract_number_field(rdkit_block, "exact_mw");
        props.logp          = extract_number_field(rdkit_block, "logp");
        props.tpsa          = extract_number_field(rdkit_block, "tpsa");
        props.hbd           = extract_int_field(rdkit_block, "hbd");
        props.hba           = extract_int_field(rdkit_block, "hba");
        props.rotbonds      = extract_int_field(rdkit_block, "rotbonds");
        props.ar_rings      = extract_int_field(rdkit_block, "ar_rings");
        props.heavy_atoms   = extract_int_field(rdkit_block, "heavy_atom_count");
        props.formal_charge = extract_int_field(rdkit_block, "formal_charge");
    }

    // Electronic section
    std::string elec_block = extract_block(json, "electronic", '{', '}');
    if (!elec_block.empty()) {
        props.homo_ev      = extract_number_field(elec_block, "homo_ev");
        props.lumo_ev      = extract_number_field(elec_block, "lumo_ev");
        props.gap_ev       = extract_number_field(elec_block, "gap_ev");
        props.dipole_debye = extract_number_field(elec_block, "dipole_debye");
        props.has_electronic = (props.homo_ev != 0.0 || props.lumo_ev != 0.0
                                || props.gap_ev != 0.0 || props.dipole_debye != 0.0);
    }

    // Fukui indices
    std::string fukui_array = extract_block(json, "fukui", '[', ']');
    if (!fukui_array.empty()) {
        for (const auto &obj : split_json_objects(fukui_array)) {
            FukuiAtom fa;
            fa.atom_index = extract_int_field(obj, "atom_idx");
            fa.element    = extract_string_field(obj, "element");
            fa.f_plus     = extract_number_field(obj, "f_plus");
            fa.f_minus    = extract_number_field(obj, "f_minus");
            fa.f_zero     = extract_number_field(obj, "f_zero");
            props.fukui.push_back(fa);
        }
        props.has_fukui = !props.fukui.empty();
    }

    // pKa groups
    std::string pka_array = extract_block(json, "pka_groups", '[', ']');
    if (!pka_array.empty()) {
        for (const auto &obj : split_json_objects(pka_array)) {
            PkaGroup pg;
            pg.atom_index  = extract_int_field(obj, "atom_idx");
            pg.pka_est     = extract_number_field(obj, "pka_est");
            pg.pka_low     = extract_number_field(obj, "pka_low");
            pg.pka_high    = extract_number_field(obj, "pka_high");
            pg.group_name  = extract_string_field(obj, "group_name");
            pg.site_type   = extract_string_field(obj, "site_type");
            props.pka_groups.push_back(pg);
        }
    }

    // Redox section
    std::string redox_block = extract_block(json, "redox", '{', '}');
    if (!redox_block.empty()) {
        props.e_ox_nhe  = extract_number_field(redox_block, "e_ox_nhe");
        props.e_red_nhe = extract_number_field(redox_block, "e_red_nhe");
        props.e_ox_fc   = extract_number_field(redox_block, "e_ox_fc");
        props.e_red_fc  = extract_number_field(redox_block, "e_red_fc");
        props.has_redox = (props.e_ox_nhe != 0.0 || props.e_red_nhe != 0.0);
    }

    // MS adducts + isotope pattern
    std::string ms_block = extract_block(json, "ms", '{', '}');
    if (!ms_block.empty()) {
        props.ms_monoisotopic = extract_number_field(ms_block, "monoisotopic_mass");

        std::string adducts_arr = extract_block(ms_block, "adducts", '[', ']');
        if (!adducts_arr.empty()) {
            for (const auto &obj : split_json_objects(adducts_arr)) {
                MsAdduct a;
                a.name   = extract_string_field(obj, "name");
                a.mz     = extract_number_field(obj, "mz");
                a.charge = extract_int_field(obj, "charge", 1);
                a.mode   = extract_string_field(obj, "mode");
                props.ms_adducts.push_back(a);
            }
        }

        std::string iso_arr = extract_block(ms_block, "isotope_pattern", '[', ']');
        if (!iso_arr.empty()) {
            for (const auto &obj : split_json_objects(iso_arr)) {
                IsotopePeak p;
                p.mass_offset        = extract_int_field(obj, "mass_offset");
                p.relative_abundance = extract_number_field(obj, "relative_abundance");
                props.isotope_peaks.push_back(p);
            }
        }

        props.has_ms = (props.ms_monoisotopic > 0.0);
    }

    // Warnings
    props.warnings = extract_string_array(json, "warnings");

    props.valid = true;
    return props;
}

} // namespace easynmr
