#include "core/batch_parser.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace easynmr {
namespace {

std::string trim(const std::string &s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

// Split a line by the given delimiter; trims each field.
std::vector<std::string> split_line(const std::string &line, char delim) {
    std::vector<std::string> fields;
    std::istringstream ss(line);
    std::string field;
    while (std::getline(ss, field, delim)) {
        fields.push_back(trim(field));
    }
    return fields;
}

// Wrap a string in single quotes for shell use, escaping any embedded single quotes.
std::string shell_quote(const std::string &s) {
    std::string r = "'";
    for (char c : s) {
        if (c == '\'') r += "'\\''";
        else r += c;
    }
    r += "'";
    return r;
}

} // namespace

std::vector<BatchEntry> parse_batch_sdf(const std::string &content) {
    std::vector<BatchEntry> entries;
    std::istringstream ss(content);
    std::string line;

    std::vector<std::string> block_lines;
    while (std::getline(ss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line == "$$$$") {
            if (!block_lines.empty()) {
                BatchEntry entry;
                // First line of a V2000/V3000 molfile is the molecule name (may be blank).
                entry.name = trim(block_lines.front());
                std::string mol_block;
                for (const auto &bl : block_lines) mol_block += bl + "\n";
                mol_block += "$$$$\n";
                entry.input_value = mol_block;
                entry.input_format = InputFormat::Mol;
                entries.push_back(std::move(entry));
                block_lines.clear();
            }
        } else {
            block_lines.push_back(line);
        }
    }
    // Handle file without trailing "$$$$"
    if (!block_lines.empty()) {
        BatchEntry entry;
        entry.name = trim(block_lines.front());
        std::string mol_block;
        for (const auto &bl : block_lines) mol_block += bl + "\n";
        entry.input_value = mol_block;
        entry.input_format = InputFormat::Mol;
        entries.push_back(std::move(entry));
    }
    return entries;
}

std::vector<BatchEntry> parse_batch_csv(const std::string &content) {
    std::vector<BatchEntry> entries;
    std::istringstream ss(content);
    std::string line;

    // Detect delimiter from the first non-blank line.
    char delim = ',';
    std::string first_line;
    while (std::getline(ss, first_line)) {
        if (!first_line.empty() && first_line.back() == '\r') first_line.pop_back();
        if (!trim(first_line).empty()) break;
    }
    if (first_line.find('\t') != std::string::npos) delim = '\t';

    // Detect header row and column indices.
    auto fields = split_line(first_line, delim);
    int smiles_col = -1;
    int name_col = -1;
    bool has_header = false;
    for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
        const std::string lf = to_lower(fields[i]);
        if (lf == "smiles" || lf == "smi" || lf == "structure" || lf == "canonical_smiles") {
            smiles_col = i;
            has_header = true;
        }
        if (lf == "name" || lf == "compound" || lf == "id" || lf == "compound_name" ||
            lf == "molecule_name" || lf == "mol_name") {
            name_col = i;
            has_header = true;
        }
    }

    if (!has_header) {
        // Treat first line as data: column 0 = SMILES, column 1 = name.
        smiles_col = 0;
        name_col = (fields.size() > 1) ? 1 : -1;
        if (smiles_col < static_cast<int>(fields.size())) {
            const std::string smi = fields[smiles_col];
            if (!smi.empty()) {
                BatchEntry entry;
                entry.input_value = smi;
                entry.input_format = InputFormat::Smiles;
                if (name_col >= 0 && name_col < static_cast<int>(fields.size())) {
                    entry.name = fields[name_col];
                }
                entries.push_back(std::move(entry));
            }
        }
    } else if (smiles_col < 0) {
        // Header found but no SMILES column — fall back to column 0.
        smiles_col = 0;
    }

    while (std::getline(ss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string tl = trim(line);
        if (tl.empty() || tl.front() == '#') continue;
        auto row = split_line(line, delim);
        if (smiles_col >= static_cast<int>(row.size())) continue;
        const std::string smi = row[smiles_col];
        if (smi.empty()) continue;
        BatchEntry entry;
        entry.input_value = smi;
        entry.input_format = InputFormat::Smiles;
        if (name_col >= 0 && name_col < static_cast<int>(row.size())) {
            entry.name = row[name_col];
        }
        entries.push_back(std::move(entry));
    }
    return entries;
}

std::vector<BatchEntry> parse_batch_smi(const std::string &content) {
    std::vector<BatchEntry> entries;
    std::istringstream ss(content);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string tl = trim(line);
        if (tl.empty() || tl.front() == '#') continue;
        const auto ws = tl.find_first_of(" \t");
        std::string smi = (ws == std::string::npos) ? tl : tl.substr(0, ws);
        std::string name;
        if (ws != std::string::npos) name = trim(tl.substr(ws + 1));
        if (smi.empty()) continue;
        BatchEntry entry;
        entry.input_value = smi;
        entry.input_format = InputFormat::Smiles;
        entry.name = name;
        entries.push_back(std::move(entry));
    }
    return entries;
}

std::vector<BatchEntry> parse_batch_xyz(const std::string &content) {
    std::vector<BatchEntry> entries;
    std::istringstream ss(content);
    std::string line;

    while (std::getline(ss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string t = trim(line);
        if (t.empty()) continue;

        // First line of a block must be the atom count.
        int n_atoms = 0;
        try {
            const std::size_t consumed_len = 0;
            n_atoms = std::stoi(t);
            (void)consumed_len;
        } catch (...) {
            continue;
        }
        if (n_atoms <= 0 || n_atoms > 100000) continue;

        std::string xyz_block = t + "\n";
        std::string name;

        // Comment / name line
        if (std::getline(ss, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            name = trim(line);
            xyz_block += line + "\n";
        } else {
            break;
        }

        // Coordinate lines
        bool valid = true;
        for (int i = 0; i < n_atoms; ++i) {
            if (!std::getline(ss, line)) {
                valid = false;
                break;
            }
            if (!line.empty() && line.back() == '\r') line.pop_back();
            xyz_block += line + "\n";
        }

        if (valid) {
            BatchEntry entry;
            entry.input_value = xyz_block;
            entry.input_format = InputFormat::Xyz;
            entry.name = name;
            entries.push_back(std::move(entry));
        }
    }
    return entries;
}

std::vector<BatchEntry> parse_batch_cdxml(const std::string &path,
                                           const std::string &python_executable) {
    if (python_executable.empty()) return {};

    // Minimal Python script: extracts SMILES from each molecule in a CDXML file.
    // Requires RDKit >= 2022.03 (Chem.MolsFromCDXMLFile).
    static const char *script_src =
        "import sys\n"
        "try:\n"
        "    from rdkit import Chem\n"
        "    mols = Chem.MolsFromCDXMLFile(sys.argv[1])\n"
        "    for mol in mols:\n"
        "        if mol is None: continue\n"
        "        try:\n"
        "            nm = mol.GetProp('_Name') if mol.HasProp('_Name') else ''\n"
        "            nm = nm.strip()\n"
        "            smi = Chem.MolToSmiles(mol)\n"
        "            if smi:\n"
        "                print(smi + ('\\t' + nm if nm else ''))\n"
        "        except Exception:\n"
        "            pass\n"
        "except Exception as e:\n"
        "    print(str(e), file=sys.stderr)\n"
        "    sys.exit(1)\n";

    // Write helper script to temp directory.
    const std::filesystem::path tmp_script =
        std::filesystem::temp_directory_path() / "easynmr_cdxml_helper.py";
    {
        std::ofstream out(tmp_script);
        if (!out.is_open()) return {};
        out << script_src;
    }

    const std::string cmd = shell_quote(python_executable) + " " +
                            shell_quote(tmp_script.string()) + " " +
                            shell_quote(path) + " 2>/dev/null";

    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::filesystem::remove(tmp_script);
        return {};
    }

    std::string output;
    char buf[512];
    while (fgets(buf, sizeof(buf), pipe) != nullptr) output += buf;
    pclose(pipe);
    std::filesystem::remove(tmp_script);

    if (output.empty()) return {};
    return parse_batch_smi(output);  // output is tab-separated SMILES\tname, handled by smi parser
}

std::vector<BatchEntry> parse_batch_file(const std::string &path,
                                          const std::string &python_executable) {
    const std::string ext = to_lower(std::filesystem::path(path).extension().string());

    if (ext == ".cdxml" || ext == ".cdx") {
        // Binary .cdx is not supported by RDKit; .cdxml (XML) is.
        // Both extensions are routed through the Python helper; .cdx will return empty.
        return parse_batch_cdxml(path, python_executable);
    }

    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) return {};
    const std::string content((std::istreambuf_iterator<char>(ifs)),
                               std::istreambuf_iterator<char>());

    if (ext == ".sdf") return parse_batch_sdf(content);
    if (ext == ".xyz") return parse_batch_xyz(content);
    if (ext == ".csv") return parse_batch_csv(content);

    if (ext == ".txt") {
        // Try CSV first (handles tab/comma-delimited structured files).
        auto csv_result = parse_batch_csv(content);
        if (!csv_result.empty()) return csv_result;
        return parse_batch_smi(content);
    }

    // .smi, .smiles, or anything else — SMI parser.
    return parse_batch_smi(content);
}

} // namespace easynmr
