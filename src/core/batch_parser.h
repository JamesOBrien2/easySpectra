#pragma once

#include "core/job.h"

#include <string>
#include <vector>

namespace easynmr {

struct BatchEntry {
    std::string input_value;
    InputFormat input_format = InputFormat::Smiles;
    std::string name;
};

// Parse a multi-molecule SDF file. Splits on "$$$$" lines.
std::vector<BatchEntry> parse_batch_sdf(const std::string &content);

// Parse a CSV or TXT with SMILES + optional name columns.
// Recognises headers: smiles/smi/structure (SMILES) and name/compound/id (name).
// Auto-detects comma or tab delimiter. Falls back to column 0 = SMILES, column 1 = name.
std::vector<BatchEntry> parse_batch_csv(const std::string &content);

// Parse a .smi file — one SMILES per line, optional name after first whitespace.
std::vector<BatchEntry> parse_batch_smi(const std::string &content);

// Parse a multi-XYZ file (concatenated XYZ blocks).
// Each block: first line = atom count, second line = molecule name/comment, then coordinates.
std::vector<BatchEntry> parse_batch_xyz(const std::string &content);

// Extract molecules from a ChemDraw XML (.cdxml) file via Python+RDKit.
// Requires RDKit >= 2022.03. Returns empty if python_executable is empty or unavailable.
std::vector<BatchEntry> parse_batch_cdxml(const std::string &path,
                                           const std::string &python_executable);

// Dispatch by file extension. python_executable is required only for .cdxml files.
// Supported extensions:
//   .sdf            — multi-molecule SDF
//   .csv / .txt     — SMILES + name columns (comma or tab delimited)
//   .smi / .smiles  — one SMILES per line
//   .xyz            — multi-XYZ blocks
//   .cdxml          — ChemDraw XML (requires python_executable + RDKit)
std::vector<BatchEntry> parse_batch_file(const std::string &path,
                                          const std::string &python_executable = {});

} // namespace easynmr
