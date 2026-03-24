#include "gui/app_window.h"

#include <FL/Fl.H>
#include <FL/Fl_Box.H>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

namespace easynmr {
namespace {

Fl_Color ui(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

struct PeakRow {
    int group_id = 0;
    double center_ppm = 0.0;
    std::string multiplicity;
    double integral = 0.0;
};

struct WorkflowStep {
    std::string label;
    std::string method;
};

std::string normalize_nucleus_label(const std::string &raw);
bool is_proton_nucleus(const std::string &nucleus);
std::string nucleus_symbol(const std::string &nucleus);

std::string truncate_text(const std::string &value, std::size_t max_len) {
    if (value.size() <= max_len) {
        return value;
    }
    return value.substr(0, max_len - 1) + "...";
}

std::string shell_quote(const std::string &value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
        if (c == '"' || c == '\\') {
            escaped.push_back('\\');
        }
        escaped.push_back(c);
    }
    return "\"" + escaped + "\"";
}

bool is_executable_file(const std::string &path) {
    if (path.empty()) {
        return false;
    }
    return access(path.c_str(), X_OK) == 0;
}

std::string guess_xyzedit_workdir(const std::string &editor_bin) {
    namespace fs = std::filesystem;

    std::vector<fs::path> candidates;
    const fs::path bin_path(editor_bin);
    if (!bin_path.empty()) {
        const fs::path bin_dir = bin_path.parent_path();
        if (!bin_dir.empty()) {
            candidates.push_back(bin_dir);
            if (!bin_dir.parent_path().empty()) {
                candidates.push_back(bin_dir.parent_path());
            }
        }
    }

    candidates.emplace_back("/Users/user/Desktop/PhD/GitHub/xyzedit");

    for (const auto &candidate : candidates) {
        if (candidate.empty()) {
            continue;
        }
        std::error_code ec;
        const fs::path data_dir = candidate / "data" / "fragments";
        if (fs::exists(data_dir, ec) && !ec) {
            return candidate.string();
        }
    }
    return {};
}

std::string find_xyzedit_binary() {
    const std::vector<std::string> fixed_candidates = {
        "/Users/user/Desktop/PhD/GitHub/xyzedit/build/xyzedit-gui",
        "/Users/user/Desktop/PhD/GitHub/xyzedit/xyzedit-gui",
        "/Users/user/Desktop/PhD/GitHub/xyzedit/build/xyzedit",
        "/Users/user/Desktop/PhD/GitHub/xyzedit/xyzedit",
    };
    for (const auto &candidate : fixed_candidates) {
        if (is_executable_file(candidate)) {
            return candidate;
        }
    }

    const char *path_env = std::getenv("PATH");
    if (path_env == nullptr) {
        return {};
    }

    std::stringstream ss(path_env);
    std::string dir;
    while (std::getline(ss, dir, ':')) {
        if (dir.empty()) {
            continue;
        }
        const std::string gui_candidate = dir + "/xyzedit-gui";
        if (is_executable_file(gui_candidate)) {
            return gui_candidate;
        }
        const std::string cli_candidate = dir + "/xyzedit";
        if (is_executable_file(cli_candidate)) {
            return cli_candidate;
        }
    }
    return {};
}

std::string format_label_for(const QueuedJob &job, bool active) {
    std::string state = "PENDING";
    if (job.status == "running") {
        state = "RUNNING";
    } else if (job.status == "done") {
        state = "DONE";
    } else if (job.status == "failed") {
        state = "FAILED";
    }

    std::ostringstream oss;
    if (active) {
        oss << ">> ";
    } else {
        oss << "   ";
    }
    oss << "[" << state << "] "
        << job.id << " | "
        << truncate_text(job.config.job_name, 10) << " | "
        << to_string(job.config.input_format) << " | "
        << job.config.nucleus;
    if (job.status == "running" && !job.message.empty()) {
        oss << " | " << truncate_text(job.message, 26);
    } else if (job.status == "failed" && !job.message.empty()) {
        oss << " | " << truncate_text(job.message, 14);
    }
    return oss.str();
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

bool try_parse_int(const std::string &text, int &out) {
    const std::string clean = trim_copy(text);
    if (clean.empty()) {
        return false;
    }
    try {
        std::size_t consumed = 0;
        const int parsed = std::stoi(clean, &consumed);
        if (consumed != clean.size()) {
            return false;
        }
        out = parsed;
        return true;
    } catch (...) {
        return false;
    }
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

std::vector<PeakRow> read_peak_rows(const std::string &path) {
    std::vector<PeakRow> rows;
    std::ifstream input(path);
    if (!input) {
        return rows;
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
        if (fields.size() < 6) {
            continue;
        }

        PeakRow row;
        if (!try_parse_int(fields[0], row.group_id)) {
            continue;
        }
        if (!try_parse_double(fields[1], row.center_ppm)) {
            continue;
        }
        row.multiplicity = trim_copy(fields[3]);
        if (!try_parse_double(fields[5], row.integral)) {
            continue;
        }
        rows.push_back(row);
    }

    return rows;
}

std::unordered_map<int, std::vector<int>> read_assignments(const std::string &path) {
    std::unordered_map<int, std::vector<int>> out;
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

        int group_id = 0;
        if (!try_parse_int(fields[0], group_id)) {
            continue;
        }
        std::vector<int> atoms;
        std::stringstream ss(fields[3]);
        std::string token;
        while (std::getline(ss, token, ',')) {
            int atom_idx = 0;
            if (try_parse_int(token, atom_idx)) {
                atoms.push_back(atom_idx);
            }
        }
        out[group_id] = std::move(atoms);
    }

    return out;
}

std::vector<StructureAtom> read_structure_atoms(const std::string &path) {
    std::vector<StructureAtom> atoms;
    std::ifstream input(path);
    if (!input) {
        return atoms;
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
        if (fields.size() < 5) {
            continue;
        }

        StructureAtom atom;
        if (!try_parse_int(fields[0], atom.atom_index)) {
            continue;
        }
        atom.element = trim_copy(fields[1]);
        if (!try_parse_double(fields[2], atom.x)) {
            continue;
        }
        if (!try_parse_double(fields[3], atom.y)) {
            continue;
        }

        std::stringstream hs(fields[4]);
        std::string token;
        while (std::getline(hs, token, ';')) {
            int h_idx = 0;
            if (try_parse_int(token, h_idx)) {
                atom.attached_hydrogens.push_back(h_idx);
            }
        }
        atoms.push_back(std::move(atom));
    }

    return atoms;
}

std::vector<StructureBond> read_structure_bonds(const std::string &path) {
    std::vector<StructureBond> bonds;
    std::ifstream input(path);
    if (!input) {
        return bonds;
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
        if (fields.size() < 3) {
            continue;
        }

        StructureBond bond;
        if (!try_parse_int(fields[0], bond.atom_a)) {
            continue;
        }
        if (!try_parse_int(fields[1], bond.atom_b)) {
            continue;
        }
        if (!try_parse_int(fields[2], bond.order)) {
            continue;
        }
        if (bond.order < 1) {
            bond.order = 1;
        }
        bonds.push_back(std::move(bond));
    }

    return bonds;
}

std::map<std::string, NucleusResultFiles> read_spectra_manifest(const std::string &path) {
    std::map<std::string, NucleusResultFiles> out;
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

        const std::string nucleus = normalize_nucleus_label(fields[0]);
        NucleusResultFiles files;
        files.spectrum_csv = trim_copy(fields[1]);
        files.peaks_csv = trim_copy(fields[2]);
        files.assignments_csv = trim_copy(fields[3]);
        if (!files.spectrum_csv.empty()) {
            out[nucleus] = std::move(files);
        }
    }

    return out;
}

std::vector<ReferencePeak> reference_peaks_for_solvent(const std::string &solvent_raw, const std::string &nucleus_raw) {
    std::string solvent = solvent_raw;
    std::transform(solvent.begin(), solvent.end(), solvent.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    std::string nucleus = nucleus_raw;
    std::transform(nucleus.begin(), nucleus.end(), nucleus.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    std::vector<ReferencePeak> peaks;
    if (nucleus == "13c") {
        peaks.push_back({0.00, "TMS 0.00"});
        if (solvent == "cdcl3" || solvent == "chcl3" || solvent == "chloroform") {
            peaks.push_back({77.16, "CDCl3 77.16"});
        } else if (solvent == "dmso") {
            peaks.push_back({39.52, "DMSO 39.52"});
        }
    } else if (nucleus == "19f") {
        peaks.push_back({0.00, "CFCl3 0.00"});
    } else {
        peaks.push_back({0.00, "TMS 0.00"});
        if (solvent == "cdcl3" || solvent == "chcl3" || solvent == "chloroform") {
            peaks.push_back({7.26, "CDCl3 7.26"});
            peaks.push_back({1.56, "H2O 1.56"});
        } else if (solvent == "dmso") {
            peaks.push_back({2.50, "DMSO 2.50"});
            peaks.push_back({3.33, "H2O 3.33"});
        } else if (solvent == "h2o" || solvent == "water") {
            peaks.push_back({4.79, "H2O 4.79"});
        }
    }

    return peaks;
}

const std::vector<WorkflowStep> &workflow_steps() {
    static const std::vector<WorkflowStep> steps = {
        {"Input validation", "Parse format and sanitize molecule"},
        {"2D atom mapping", "Build depiction and atom/hydrogen map"},
        {"Conformer generation", "ETKDG + MMFF pre-opt"},
        {"Conformer optimization", "xTB GFN2-xTB + ALPB (MMFF fallback)"},
        {"Boltzmann selection", "99% population / 6 kcal mol^-1 window"},
        {"NMR parameter estimation", "Target-nucleus shifts + couplings (first-order)"},
        {"Spectrum simulation", "Apply line shape and splitting model"},
        {"Outputs and audit", "Write spectra, assignments, metadata"},
    };
    return steps;
}

int workflow_stage_index(const std::string &stage_raw) {
    std::string stage = stage_raw;
    std::transform(stage.begin(), stage.end(), stage.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (stage == "launch" || stage == "initializing" || stage == "parse_input") {
        return 0;
    }
    if (stage == "structure_2d") {
        return 1;
    }
    if (stage == "conformer_generation") {
        return 2;
    }
    if (stage == "conformer_optimization") {
        return 3;
    }
    if (stage == "boltzmann_selection") {
        return 4;
    }
    if (stage == "nmr_parameter_estimation") {
        return 5;
    }
    if (stage == "spectrum_simulation") {
        return 6;
    }
    if (stage == "write_outputs" || stage == "done") {
        return 7;
    }
    return -1;
}

std::string normalize_nucleus_label(const std::string &raw) {
    std::string nucleus = raw;
    std::transform(nucleus.begin(), nucleus.end(), nucleus.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (nucleus == "13c" || nucleus == "c13" || nucleus == "c") {
        return "13C";
    }
    if (nucleus == "19f" || nucleus == "f19" || nucleus == "f") {
        return "19F";
    }
    return "1H";
}

bool is_proton_nucleus(const std::string &nucleus) {
    return normalize_nucleus_label(nucleus) == "1H";
}

std::string nucleus_symbol(const std::string &nucleus) {
    const std::string normalized = normalize_nucleus_label(nucleus);
    if (normalized == "13C") {
        return "C";
    }
    if (normalized == "19F") {
        return "F";
    }
    return "H";
}

} // namespace

AppWindow::AppWindow(int w, int h, const char *title)
    : Fl_Double_Window(w, h, title) {
    color(ui(237, 242, 248));
    begin();

    auto *top_bar = new Fl_Box(FL_FLAT_BOX, 0, 0, w, 36, "");
    top_bar->color(ui(49, 58, 72));
    auto *top_accent = new Fl_Box(FL_FLAT_BOX, 0, 34, w, 2, "");
    top_accent->color(ui(154, 184, 204));
    auto *top_title = new Fl_Box(12, 6, 420, 24, "EasyNMR  |  Local NMR Predictor");
    top_title->box(FL_NO_BOX);
    top_title->labelfont(FL_HELVETICA_BOLD);
    top_title->labelsize(12);
    top_title->labelcolor(ui(238, 245, 252));
    top_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    const int panel_x = 10;
    const int panel_y = 46;
    const int panel_w = 334;
    const int panel_h = h - panel_y - 12;

    auto *left_panel_bg = new Fl_Box(FL_THIN_UP_BOX, panel_x, panel_y, panel_w, panel_h, "");
    left_panel_bg->color(ui(248, 251, 255));
    auto *left_panel_accent = new Fl_Box(FL_FLAT_BOX, panel_x, panel_y, 4, panel_h, "");
    left_panel_accent->color(ui(154, 184, 204));

    auto *left_title = new Fl_Box(panel_x + 10, panel_y + 8, panel_w - 20, 22, "Molecule Setup");
    left_title->box(FL_NO_BOX);
    left_title->labelfont(FL_HELVETICA_BOLD);
    left_title->labelsize(11);
    left_title->labelcolor(ui(67, 77, 92));
    left_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    auto *job_name_label = new Fl_Box(panel_x + 10, panel_y + 34, 120, 16, "Job name");
    job_name_label->box(FL_NO_BOX);
    job_name_label->labelsize(12);
    job_name_label->labelcolor(ui(86, 97, 112));
    job_name_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    job_name_input_ = new Fl_Input(panel_x + 10, panel_y + 50, panel_w - 20, 26);
    job_name_input_->value("untitled");
    job_name_input_->box(FL_DOWN_BOX);
    job_name_input_->color(ui(255, 255, 255));

    auto *format_label = new Fl_Box(panel_x + 10, panel_y + 82, 120, 16, "Format");
    format_label->box(FL_NO_BOX);
    format_label->labelsize(12);
    format_label->labelcolor(ui(86, 97, 112));
    format_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    format_choice_ = new Fl_Choice(panel_x + 10, panel_y + 98, 156, 26);
    format_choice_->add("smiles");
    format_choice_->add("mol");
    format_choice_->add("sdf");
    format_choice_->add("xyz");
    format_choice_->value(0);
    format_choice_->box(FL_DOWN_BOX);
    format_choice_->color(ui(255, 255, 255));
    format_choice_->callback(on_preview_cb, this);

    auto *solvent_label = new Fl_Box(panel_x + 178, panel_y + 82, 120, 16, "Solvent");
    solvent_label->box(FL_NO_BOX);
    solvent_label->labelsize(12);
    solvent_label->labelcolor(ui(86, 97, 112));
    solvent_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    solvent_choice_ = new Fl_Choice(panel_x + 178, panel_y + 98, 146, 26);
    solvent_choice_->add("cdcl3");
    solvent_choice_->add("dmso");
    solvent_choice_->add("h2o");
    solvent_choice_->value(0);
    solvent_choice_->box(FL_DOWN_BOX);
    solvent_choice_->color(ui(255, 255, 255));
    solvent_choice_->callback(on_preview_cb, this);

    auto *line_shape_label = new Fl_Box(panel_x + 10, panel_y + 128, 120, 16, "Line shape");
    line_shape_label->box(FL_NO_BOX);
    line_shape_label->labelsize(12);
    line_shape_label->labelcolor(ui(86, 97, 112));
    line_shape_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    line_shape_choice_ = new Fl_Choice(panel_x + 10, panel_y + 144, 156, 26);
    line_shape_choice_->add("lorentzian");
    line_shape_choice_->add("gaussian");
    line_shape_choice_->add("voigt");
    line_shape_choice_->value(0);
    line_shape_choice_->box(FL_DOWN_BOX);
    line_shape_choice_->color(ui(255, 255, 255));

    auto *nucleus_label = new Fl_Box(panel_x + 178, panel_y + 128, 120, 16, "Nucleus");
    nucleus_label->box(FL_NO_BOX);
    nucleus_label->labelsize(12);
    nucleus_label->labelcolor(ui(86, 97, 112));
    nucleus_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    nucleus_choice_ = new Fl_Choice(panel_x + 178, panel_y + 144, 146, 26);
    nucleus_choice_->add("auto (all present)");
    nucleus_choice_->value(0);
    nucleus_choice_->box(FL_DOWN_BOX);
    nucleus_choice_->color(ui(255, 255, 255));
    nucleus_choice_->deactivate();

    auto *input_label = new Fl_Box(panel_x + 10, panel_y + 176, panel_w - 20, 16, "SMILES / structure text");
    input_label->box(FL_NO_BOX);
    input_label->labelsize(12);
    input_label->labelcolor(ui(86, 97, 112));
    input_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    input_box_ = new Fl_Multiline_Input(panel_x + 10, panel_y + 192, panel_w - 20, 64);
    input_box_->value("CCO");
    input_box_->box(FL_DOWN_BOX);
    input_box_->color(ui(255, 255, 255));
    input_box_->when(FL_WHEN_RELEASE | FL_WHEN_ENTER_KEY);
    input_box_->callback(on_preview_cb, this);

    auto *structure_title = new Fl_Box(panel_x + 10, panel_y + 264, panel_w - 214, 20, "Interactive Structure");
    structure_title->box(FL_NO_BOX);
    structure_title->labelsize(12);
    structure_title->labelcolor(ui(86, 97, 112));
    structure_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    edit_button_ = new Fl_Button(panel_x + panel_w - 196, panel_y + 262, 104, 24, "Edit xyzedit");
    edit_button_->callback(on_edit_structure_cb, this);
    edit_button_->box(FL_UP_BOX);
    edit_button_->color(ui(207, 218, 233));
    edit_button_->labelcolor(ui(63, 73, 86));
    edit_button_->labelfont(FL_HELVETICA_BOLD);
    edit_button_->labelsize(11);

    preview_button_ = new Fl_Button(panel_x + panel_w - 86, panel_y + 262, 76, 24, "Preview");
    preview_button_->callback(on_preview_cb, this);
    preview_button_->box(FL_UP_BOX);
    preview_button_->color(ui(198, 211, 228));
    preview_button_->labelcolor(ui(63, 73, 86));
    preview_button_->labelfont(FL_HELVETICA_BOLD);
    preview_button_->labelsize(11);

    structure_widget_ = new StructureWidget(panel_x + 10, panel_y + 288, panel_w - 20, 160, nullptr);
    structure_widget_->set_on_atom_selected([this](int atom_index, const std::vector<int> &hydrogens) {
        on_structure_atom_picked(atom_index, hydrogens);
    });

    queue_button_ = new Fl_Button(panel_x + 10, panel_y + 454, 100, 30, "Queue");
    queue_button_->callback(on_queue_job_cb, this);
    queue_button_->box(FL_UP_BOX);
    queue_button_->color(ui(218, 225, 236));
    queue_button_->labelcolor(ui(63, 73, 86));
    queue_button_->labelfont(FL_HELVETICA_BOLD);

    start_button_ = new Fl_Button(panel_x + 116, panel_y + 454, 110, 30, "Run Pending");
    start_button_->callback(on_start_queue_cb, this);
    start_button_->box(FL_UP_BOX);
    start_button_->color(ui(168, 205, 194));
    start_button_->labelcolor(ui(52, 66, 66));
    start_button_->labelfont(FL_HELVETICA_BOLD);

    cancel_button_ = new Fl_Button(panel_x + 232, panel_y + 454, 92, 30, "Cancel");
    cancel_button_->callback(on_cancel_cb, this);
    cancel_button_->box(FL_UP_BOX);
    cancel_button_->color(ui(223, 229, 238));
    cancel_button_->labelcolor(ui(88, 98, 112));

    auto *queue_title = new Fl_Box(panel_x + 10, panel_y + 492, panel_w - 20, 20, "Queue (up to 50 jobs)");
    queue_title->box(FL_NO_BOX);
    queue_title->labelsize(12);
    queue_title->labelcolor(ui(86, 97, 112));
    queue_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    queue_browser_ = new Fl_Hold_Browser(panel_x + 10, panel_y + 516, panel_w - 20, 140);
    queue_browser_->callback(on_select_job_cb, this);
    queue_browser_->box(FL_DOWN_BOX);
    queue_browser_->color(ui(255, 255, 255));
    queue_browser_->textsize(11);
    queue_browser_->selection_color(ui(224, 236, 248));

    run_selected_button_ = new Fl_Button(panel_x + 10, panel_y + 666, 130, 26, "Run Selected");
    run_selected_button_->callback(on_run_selected_cb, this);
    run_selected_button_->box(FL_UP_BOX);
    run_selected_button_->color(ui(186, 204, 229));
    run_selected_button_->labelcolor(ui(59, 73, 94));
    run_selected_button_->labelfont(FL_HELVETICA_BOLD);

    status_box_ = new Fl_Box(panel_x + 146, panel_y + 666, panel_w - 156, 26, "Idle");
    status_box_->box(FL_FLAT_BOX);
    status_box_->color(ui(229, 236, 246));
    status_box_->labelsize(12);
    status_box_->labelcolor(ui(67, 77, 92));
    status_box_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    const int right_x = panel_x + panel_w + 10;
    const int right_y = panel_y;
    const int right_w = w - right_x - 10;
    const int spectrum_h = h - right_y - 250;

    auto *spectrum_bg = new Fl_Box(FL_THIN_UP_BOX, right_x, right_y, right_w, spectrum_h, "");
    spectrum_bg->color(ui(250, 252, 255));

    auto *spectrum_title = new Fl_Box(right_x + 10, right_y + 8, 220, 20, "Spectrum");
    spectrum_title->box(FL_NO_BOX);
    spectrum_title->labelfont(FL_HELVETICA_BOLD);
    spectrum_title->labelsize(12);
    spectrum_title->labelcolor(ui(68, 78, 92));
    spectrum_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    spectrum_widget_ = new SpectrumWidget(right_x + 8, right_y + 30, right_w - 16, spectrum_h - 38, "Spectrum");
    spectrum_widget_->set_on_peak_selected([this](int group_id) { on_peak_picked(group_id); });

    spectrum_nucleus_choice_ = new Fl_Choice(right_x + 210, right_y + 6, 110, 24);
    spectrum_nucleus_choice_->box(FL_DOWN_BOX);
    spectrum_nucleus_choice_->color(ui(255, 255, 255));
    spectrum_nucleus_choice_->add("1H");
    spectrum_nucleus_choice_->value(0);
    spectrum_nucleus_choice_->callback(on_select_spectrum_nucleus_cb, this);

    reference_choice_ = new Fl_Choice(right_x + right_w - 250, right_y + 6, 240, 24);
    reference_choice_->box(FL_DOWN_BOX);
    reference_choice_->color(ui(255, 255, 255));
    reference_choice_->labelsize(11);
    reference_choice_->add("Refs: none");
    reference_choice_->value(0);
    reference_choice_->callback(on_select_reference_cb, this);

    auto *peak_title = new Fl_Box(right_x, right_y + spectrum_h + 8, right_w / 2 - 6, 20, "Peak Groups");
    peak_title->box(FL_NO_BOX);
    peak_title->labelsize(12);
    peak_title->labelcolor(ui(80, 90, 105));
    peak_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    peak_browser_ = new Fl_Hold_Browser(right_x, right_y + spectrum_h + 30, right_w / 2 - 6, 88);
    peak_browser_->callback(on_select_peak_cb, this);
    peak_browser_->box(FL_DOWN_BOX);
    peak_browser_->color(ui(255, 255, 255));
    peak_browser_->selection_color(ui(224, 236, 248));
    peak_browser_->textsize(11);

    auto *atom_title = new Fl_Box(right_x + right_w / 2 + 6, right_y + spectrum_h + 8, right_w / 2 - 6, 20, "Nucleus Assignments");
    atom_title->box(FL_NO_BOX);
    atom_title->labelsize(12);
    atom_title->labelcolor(ui(80, 90, 105));
    atom_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    atom_browser_ = new Fl_Hold_Browser(right_x + right_w / 2 + 6, right_y + spectrum_h + 30, right_w / 2 - 6, 88);
    atom_browser_->callback(on_select_atom_cb, this);
    atom_browser_->box(FL_DOWN_BOX);
    atom_browser_->color(ui(255, 255, 255));
    atom_browser_->selection_color(ui(224, 236, 248));
    atom_browser_->textsize(11);

    auto *workflow_title = new Fl_Box(right_x, right_y + spectrum_h + 126, right_w, 20, "Workflow Monitor");
    workflow_title->box(FL_NO_BOX);
    workflow_title->labelsize(12);
    workflow_title->labelcolor(ui(80, 90, 105));
    workflow_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    workflow_progress_widget_ = new WorkflowProgressWidget(right_x, right_y + spectrum_h + 148, right_w, 30, nullptr);
    {
        std::vector<std::string> step_labels;
        for (const auto &step : workflow_steps()) {
            step_labels.push_back(step.label);
        }
        workflow_progress_widget_->set_steps(std::move(step_labels));
        workflow_progress_widget_->set_progress_state(-1, -1, false, false, false, 0.0, "Workflow idle");
    }

    workflow_browser_ = new Fl_Hold_Browser(right_x, right_y + spectrum_h + 180, right_w, 52);
    workflow_browser_->box(FL_DOWN_BOX);
    workflow_browser_->color(ui(255, 255, 255));
    workflow_browser_->textsize(11);
    workflow_browser_->selection_color(ui(255, 255, 255));

    end();
    Fl::add_timeout(0.2, on_ui_tick_cb, this);
    preview_current_input(false);
    refresh_workflow_browser(nullptr);
}

AppWindow::~AppWindow() {
    cancel_requested_ = true;
    Fl::remove_timeout(on_debounced_preview_cb, this);
    Fl::remove_timeout(on_ui_tick_cb, this);
    if (preview_future_active_ && preview_future_.valid()) {
        preview_future_.wait();
        (void)preview_future_.get();
        preview_future_active_ = false;
    }
    if (worker_thread_.joinable()) {
        worker_thread_.join();
    }
}

void AppWindow::on_queue_job_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->queue_current_input();
}

void AppWindow::on_start_queue_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->start_queue();
}

void AppWindow::on_cancel_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->cancel_queue();
}

void AppWindow::on_run_selected_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->run_selected();
}

void AppWindow::on_preview_cb(Fl_Widget *widget, void *userdata) {
    auto *self = static_cast<AppWindow *>(userdata);
    if (widget == self->preview_button_) {
        if (self->preview_debounce_pending_) {
            Fl::remove_timeout(on_debounced_preview_cb, self);
            self->preview_debounce_pending_ = false;
            self->preview_debounce_show_status_ = false;
        }
        self->preview_current_input(true);
        return;
    }
    self->schedule_preview(0.55, false);
}

void AppWindow::on_edit_structure_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->edit_current_structure();
}

void AppWindow::on_select_job_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_job();
}

void AppWindow::on_select_peak_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_peak();
}

void AppWindow::on_select_atom_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_atom();
}

void AppWindow::on_select_reference_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_reference();
}

void AppWindow::on_select_spectrum_nucleus_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_spectrum_nucleus();
}

void AppWindow::on_debounced_preview_cb(void *userdata) {
    auto *self = static_cast<AppWindow *>(userdata);
    self->preview_debounce_pending_ = false;
    const bool show_status = self->preview_debounce_show_status_;
    self->preview_debounce_show_status_ = false;
    self->preview_current_input(show_status);
}

void AppWindow::on_worker_awake(void *userdata) {
    auto *self = static_cast<AppWindow *>(userdata);
    self->refresh_queue_browser();
    const int selected = self->queue_browser_ != nullptr ? self->queue_browser_->value() : 0;
    if (selected > 0) {
        bool should_reload_selected = false;
        {
            std::lock_guard<std::mutex> lock(self->jobs_mutex_);
            const std::size_t idx = static_cast<std::size_t>(selected - 1);
            if (idx < self->jobs_.size()) {
                const auto &job = self->jobs_[idx];
                should_reload_selected = (!self->worker_running_ && (job.status == "done" || job.status == "failed"));
            }
        }
        if (should_reload_selected) {
            self->load_selected_job_visuals();
        }
    }
    self->redraw();
}

void AppWindow::on_ui_tick_cb(void *userdata) {
    auto *self = static_cast<AppWindow *>(userdata);
    self->poll_preview_async();
    bool has_running_work = self->worker_running_.load();
    if (!has_running_work) {
        std::lock_guard<std::mutex> lock(self->jobs_mutex_);
        for (const auto &job : self->jobs_) {
            if (job.status == "running") {
                has_running_work = true;
                break;
            }
        }
    }
    if (has_running_work) {
        self->refresh_queue_browser();
        self->redraw();
    }
    Fl::repeat_timeout(0.2, on_ui_tick_cb, self);
}

void AppWindow::queue_current_input() {
    QueuedJob job;
    job.status = "pending";

    job.config.job_name = job_name_input_->value();
    job.config.input_value = input_box_->value();

    const char *fmt = format_choice_->text(format_choice_->value());
    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    const char *shape = line_shape_choice_->text(line_shape_choice_->value());

    job.config.input_format = input_format_from_string(fmt ? fmt : "smiles");
    job.config.solvent = solvent ? solvent : "cdcl3";
    job.config.nucleus = "auto";
    job.config.line_shape = shape ? shape : "lorentzian";

    if (job.config.input_value.empty()) {
        status_box_->label("Cannot queue: input is empty");
        return;
    }

    std::size_t queued_index = 0;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (jobs_.size() >= 50) {
            status_box_->label("Queue limit reached (50)");
            return;
        }
        job.id = "q" + std::to_string(jobs_.size() + 1);
        jobs_.push_back(job);
        queued_index = jobs_.size() - 1;
        if (worker_running_ && run_scope_ == RunScope::PendingQueue) {
            ++run_total_;
        }
    }

    status_box_->label("Job queued; generating structure preview");
    refresh_queue_browser();
    queue_browser_->value(static_cast<int>(queued_index + 1));
    request_preview_for_job_index(queued_index, true);
}

void AppWindow::start_queue() {
    if (worker_running_) {
        status_box_->label("Queue already running");
        return;
    }

    int pending_count = 0;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        for (const auto &job : jobs_) {
            if (job.status == "pending") {
                ++pending_count;
            }
        }
        run_scope_ = RunScope::PendingQueue;
        selected_run_index_ = static_cast<std::size_t>(-1);
        run_total_ = pending_count;
        run_done_ = 0;
        active_job_index_ = -1;
    }

    if (pending_count <= 0) {
        status_box_->label("No pending jobs to run");
        return;
    }

    cancel_requested_ = false;
    worker_running_ = true;
    status_box_->copy_label(("Run all pending: " + std::to_string(pending_count)).c_str());

    if (worker_thread_.joinable()) {
        worker_thread_.join();
    }

    worker_thread_ = std::thread([this]() { run_worker_loop(); });
}

void AppWindow::cancel_queue() {
    cancel_requested_ = true;
    status_box_->label("Cancel requested: stopping after active job");
}

void AppWindow::run_selected() {
    if (worker_running_) {
        status_box_->label("Queue already running");
        return;
    }

    const int selected = queue_browser_->value();
    if (selected <= 0) {
        status_box_->label("Select a queue job first");
        return;
    }

    std::size_t selected_index = static_cast<std::size_t>(selected - 1);
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (selected_index >= jobs_.size()) {
            return;
        }

        auto &job = jobs_[selected_index];
        if (job.status == "running") {
            status_box_->label("Selected job is already running");
            return;
        }
        // Re-run selected job from scratch even if it was done/failed.
        job.status = "pending";
        job.message.clear();
        job.output_dir.clear();
        job.spectrum_csv.clear();
        job.peaks_csv.clear();
        job.assignments_csv.clear();
        job.spectra_manifest_csv.clear();
        job.structure_svg.clear();
        job.structure_atoms_csv.clear();
        job.structure_bonds_csv.clear();
        job.structure_xyz.clear();
        job.progress_stage.clear();
        job.progress_message.clear();
        job.progress_fraction = 0.0;

        run_scope_ = RunScope::SelectedOnly;
        selected_run_index_ = selected_index;
        run_total_ = 1;
        run_done_ = 0;
        active_job_index_ = -1;
    }

    cancel_requested_ = false;
    worker_running_ = true;
    status_box_->copy_label(("Run selected: q" + std::to_string(selected)).c_str());

    if (worker_thread_.joinable()) {
        worker_thread_.join();
    }
    worker_thread_ = std::thread([this]() { run_worker_loop(); });

    refresh_queue_browser();
}

void AppWindow::schedule_preview(double delay_seconds, bool show_status) {
    if (worker_running_) {
        return;
    }

    preview_debounce_show_status_ = preview_debounce_show_status_ || show_status;
    if (preview_debounce_pending_) {
        Fl::remove_timeout(on_debounced_preview_cb, this);
    }
    preview_debounce_pending_ = true;
    Fl::add_timeout(delay_seconds, on_debounced_preview_cb, this);
}

void AppWindow::launch_preview_async(const JobConfig &cfg, bool show_status) {
    if (preview_future_active_) {
        preview_pending_cfg_ = cfg;
        preview_pending_ = true;
        preview_pending_show_status_ = preview_pending_show_status_ || show_status;
        if (show_status) {
            status_box_->label("Preview queued...");
        }
        return;
    }

    preview_inflight_cfg_ = cfg;
    preview_inflight_show_status_ = show_status;
    preview_future_ = std::async(std::launch::async, [this, cfg]() { return pipeline_.preview(cfg); });
    preview_future_active_ = true;
    if (show_status) {
        status_box_->label("Rendering preview...");
    }
}

void AppWindow::poll_preview_async() {
    if (!preview_future_active_ || !preview_future_.valid()) {
        return;
    }
    if (preview_future_.wait_for(std::chrono::milliseconds(0)) != std::future_status::ready) {
        return;
    }

    const JobOutputs preview = preview_future_.get();
    const JobConfig cfg = preview_inflight_cfg_;
    const bool show_status = preview_inflight_show_status_;
    preview_future_active_ = false;

    if (preview.status != "ok") {
        if (show_status) {
            status_box_->copy_label(("Preview unavailable: " + preview.message).c_str());
        }
    } else {
        auto structure_atoms = read_structure_atoms(preview.structure_atoms_csv);
        auto structure_bonds = read_structure_bonds(preview.structure_bonds_csv);
        if (structure_widget_ != nullptr) {
            structure_widget_->set_structure(std::move(structure_atoms), std::move(structure_bonds));
        }
        active_nucleus_ = normalize_nucleus_label(cfg.nucleus);
        apply_reference_peaks(cfg.solvent, active_nucleus_);
        spectrum_widget_->set_nucleus_label(active_nucleus_ + " NMR Spectrum");
        if (show_status) {
            status_box_->label("Preview updated from current input");
        }
    }

    if (preview_pending_) {
        const JobConfig pending_cfg = preview_pending_cfg_;
        const bool pending_show_status = preview_pending_show_status_;
        preview_pending_ = false;
        preview_pending_show_status_ = false;
        launch_preview_async(pending_cfg, pending_show_status);
    }
}

void AppWindow::preview_current_input(bool show_status) {
    if (worker_running_) {
        if (show_status) {
            status_box_->label("Preview is paused while calculations are running");
        }
        return;
    }

    JobConfig cfg;
    cfg.job_name = job_name_input_->value();
    cfg.input_value = input_box_->value();

    const char *fmt = format_choice_->text(format_choice_->value());
    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    cfg.input_format = input_format_from_string(fmt ? fmt : "smiles");
    cfg.solvent = solvent ? solvent : "cdcl3";
    cfg.nucleus = "auto";
    cfg.need_editable_xyz = false;

    if (cfg.input_value.empty() || cfg.input_format == InputFormat::Unknown) {
        return;
    }
    launch_preview_async(cfg, show_status);
}

void AppWindow::edit_current_structure() {
    if (worker_running_) {
        status_box_->label("Editor is disabled while calculations are running");
        return;
    }

    JobConfig cfg;
    cfg.job_name = job_name_input_->value();
    cfg.input_value = input_box_->value();
    const char *fmt = format_choice_->text(format_choice_->value());
    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    cfg.input_format = input_format_from_string(fmt ? fmt : "smiles");
    cfg.solvent = solvent ? solvent : "cdcl3";
    cfg.nucleus = "auto";
    cfg.need_editable_xyz = true;

    if (cfg.input_value.empty() || cfg.input_format == InputFormat::Unknown) {
        status_box_->label("Cannot open editor: invalid or empty input");
        return;
    }
    if (preview_debounce_pending_) {
        Fl::remove_timeout(on_debounced_preview_cb, this);
        preview_debounce_pending_ = false;
        preview_debounce_show_status_ = false;
    }

    status_box_->label("Preparing editable structure...");
    const JobOutputs preview = pipeline_.preview(cfg);
    if (preview.status != "ok") {
        status_box_->copy_label(("Editor prep failed: " + preview.message).c_str());
        return;
    }

    std::filesystem::path xyz_path = preview.structure_xyz;
    if (xyz_path.is_relative()) {
        xyz_path = std::filesystem::absolute(xyz_path);
    }
    if (xyz_path.empty() || !std::filesystem::exists(xyz_path)) {
        status_box_->label("Editor prep failed: missing editable XYZ");
        return;
    }

    auto structure_atoms = read_structure_atoms(preview.structure_atoms_csv);
    auto structure_bonds = read_structure_bonds(preview.structure_bonds_csv);
    if (structure_widget_ != nullptr) {
        structure_widget_->set_structure(std::move(structure_atoms), std::move(structure_bonds));
    }
    active_nucleus_ = normalize_nucleus_label(cfg.nucleus);
    apply_reference_peaks(cfg.solvent, active_nucleus_);
    spectrum_widget_->set_nucleus_label(active_nucleus_ + " NMR Spectrum");

    const std::string editor_bin = find_xyzedit_binary();
    if (editor_bin.empty()) {
        status_box_->label("xyzedit not found (install/build xyzedit-gui or xyzedit)");
        return;
    }

    const std::string editor_workdir = guess_xyzedit_workdir(editor_bin);
    status_box_->label("Opening xyzedit... close editor to import changes");
    redraw();

    std::string command;
    if (!editor_workdir.empty()) {
        command = "cd " + shell_quote(editor_workdir) + " && "
                + shell_quote(editor_bin) + " " + shell_quote(xyz_path.string());
    } else {
        command = shell_quote(editor_bin) + " " + shell_quote(xyz_path.string());
    }
    const int rc = std::system(command.c_str());
    if (rc != 0) {
        status_box_->copy_label(("xyzedit exited with code " + std::to_string(rc)).c_str());
        return;
    }

    std::ifstream edited(xyz_path);
    if (!edited) {
        status_box_->label("Edited XYZ could not be read");
        return;
    }
    std::stringstream xyz_buffer;
    xyz_buffer << edited.rdbuf();
    const std::string xyz_text = xyz_buffer.str();
    if (xyz_text.empty()) {
        status_box_->label("Edited XYZ was empty; keeping existing input");
        return;
    }

    input_box_->value(xyz_text.c_str());
    for (int i = 0; i < format_choice_->size(); ++i) {
        const char *option = format_choice_->text(i);
        if (option != nullptr && std::string(option) == "xyz") {
            format_choice_->value(i);
            break;
        }
    }

    preview_current_input(true);
    status_box_->label("Imported edited XYZ from xyzedit");
}

void AppWindow::request_preview_for_job_index(std::size_t index, bool load_visuals_if_selected) {
    JobConfig cfg;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (index >= jobs_.size()) {
            return;
        }
        cfg = jobs_[index].config;
    }

    const JobOutputs preview = pipeline_.preview(cfg);
    if (preview.status != "ok") {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (index < jobs_.size() && jobs_[index].message.empty()) {
            jobs_[index].message = "Preview failed: " + preview.message;
        }
        return;
    }

    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (index >= jobs_.size()) {
            return;
        }
        auto &job = jobs_[index];
        job.structure_svg = preview.structure_svg;
        job.structure_atoms_csv = preview.structure_atoms_csv;
        job.structure_bonds_csv = preview.structure_bonds_csv;
        job.structure_xyz = preview.structure_xyz;
    }

    refresh_queue_browser();
    if (load_visuals_if_selected) {
        const int selected = queue_browser_->value();
        if (selected > 0 && static_cast<std::size_t>(selected - 1) == index) {
            load_selected_job_visuals();
        }
    }
}

void AppWindow::on_select_job() {
    const int selected = queue_browser_->value();
    if (selected > 0) {
        const std::size_t index = static_cast<std::size_t>(selected - 1);
        bool needs_preview = false;
        {
            std::lock_guard<std::mutex> lock(jobs_mutex_);
            if (index < jobs_.size()) {
                const auto &job = jobs_[index];
                needs_preview = job.structure_atoms_csv.empty() || job.structure_bonds_csv.empty();
            }
        }
        if (needs_preview) {
            request_preview_for_job_index(index, false);
        }
    }
    load_selected_job_visuals();
}

void AppWindow::on_select_peak() {
    const int row = peak_browser_->value();
    if (row <= 0) {
        return;
    }

    auto it = peak_row_to_group_.find(row);
    if (it == peak_row_to_group_.end()) {
        return;
    }
    on_peak_picked(it->second);
}

void AppWindow::on_select_atom() {
    const int row = atom_browser_->value();
    if (row <= 0) {
        return;
    }

    auto atom_it = atom_row_to_atom_.find(row);
    if (atom_it == atom_row_to_atom_.end()) {
        return;
    }
    const int atom_idx = atom_it->second;

    std::vector<int> matched_groups;
    for (const auto &entry : group_to_atoms_) {
        const auto &atoms = entry.second;
        if (std::find(atoms.begin(), atoms.end(), atom_idx) != atoms.end()) {
            matched_groups.push_back(entry.first);
        }
    }

    std::sort(matched_groups.begin(), matched_groups.end());
    spectrum_widget_->set_selected_groups(matched_groups);
    highlight_hydrogen_rows({atom_idx});
    if (structure_widget_ != nullptr) {
        if (is_proton_nucleus(active_nucleus_)) {
            structure_widget_->set_selected_atom(-1);
            structure_widget_->set_highlight_hydrogens({atom_idx});
        } else {
            structure_widget_->set_highlight_hydrogens({});
            structure_widget_->set_selected_atom(atom_idx);
        }
    }

    if (!matched_groups.empty()) {
        const int group_id = matched_groups.front();
        auto pr = group_to_peak_row_.find(group_id);
        if (pr != group_to_peak_row_.end()) {
            peak_browser_->value(pr->second);
        }
    }

    std::ostringstream status;
    status << "Selected " << nucleus_symbol(active_nucleus_) << atom_idx;
    if (!matched_groups.empty()) {
        status << " -> groups ";
        for (std::size_t i = 0; i < matched_groups.size(); ++i) {
            status << matched_groups[i];
            if (i + 1 < matched_groups.size()) {
                status << ",";
            }
        }
    } else {
        status << " (no assigned groups)";
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_select_reference() {
    if (reference_choice_ == nullptr) {
        return;
    }
    const int selected = reference_choice_->value();
    if (selected <= 0) {
        spectrum_widget_->set_highlighted_reference(-1);
    } else {
        spectrum_widget_->set_highlighted_reference(selected - 1);
    }
}

void AppWindow::on_select_spectrum_nucleus() {
    if (spectrum_nucleus_choice_ == nullptr) {
        return;
    }
    const int selected = spectrum_nucleus_choice_->value();
    if (selected < 0 || selected >= spectrum_nucleus_choice_->size()) {
        return;
    }
    const char *label = spectrum_nucleus_choice_->text(selected);
    if (label == nullptr) {
        return;
    }
    const int selected_job = queue_browser_->value();
    if (selected_job <= 0) {
        return;
    }

    QueuedJob job;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        const std::size_t index = static_cast<std::size_t>(selected_job - 1);
        if (index >= jobs_.size()) {
            return;
        }
        job = jobs_[index];
    }

    active_nucleus_ = normalize_nucleus_label(label);
    apply_nucleus_visuals(job, active_nucleus_);
}

void AppWindow::on_peak_picked(int group_id) {
    if (group_id <= 0) {
        return;
    }

    spectrum_widget_->set_selected_group(group_id);
    if (structure_widget_ != nullptr) {
        structure_widget_->set_selected_atom(-1);
    }

    auto pr = group_to_peak_row_.find(group_id);
    if (pr != group_to_peak_row_.end()) {
        peak_browser_->value(pr->second);
    }

    auto ga = group_to_atoms_.find(group_id);
    if (ga != group_to_atoms_.end() && !ga->second.empty()) {
        highlight_hydrogen_rows(ga->second);
        if (structure_widget_ != nullptr) {
            if (is_proton_nucleus(active_nucleus_)) {
                structure_widget_->set_selected_atom(-1);
                structure_widget_->set_highlight_hydrogens(ga->second);
            } else {
                structure_widget_->set_highlight_hydrogens({});
                structure_widget_->set_selected_atom(ga->second.front());
            }
        }
        const int atom_idx = ga->second.front();
        auto ar = atom_to_row_.find(atom_idx);
        if (ar != atom_to_row_.end()) {
            atom_browser_->value(ar->second);
        }
    } else {
        highlight_hydrogen_rows({});
        if (structure_widget_ != nullptr) {
            structure_widget_->set_highlight_hydrogens({});
        }
    }

    std::ostringstream status;
    status << "Selected group " << group_id;
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_structure_atom_picked(int atom_index, const std::vector<int> &attached_hydrogens) {
    if (atom_index <= 0) {
        return;
    }

    std::vector<int> nucleus_targets = attached_hydrogens;
    if (!is_proton_nucleus(active_nucleus_)) {
        nucleus_targets = {atom_index};
    }

    if (structure_widget_ != nullptr) {
        structure_widget_->set_selected_atom(atom_index);
        if (is_proton_nucleus(active_nucleus_)) {
            structure_widget_->set_highlight_hydrogens(attached_hydrogens);
        } else {
            structure_widget_->set_highlight_hydrogens({});
        }
    }

    highlight_hydrogen_rows(nucleus_targets);

    std::vector<int> matched_groups;
    for (const auto &entry : group_to_atoms_) {
        bool intersects = false;
        for (int target_atom : nucleus_targets) {
            if (std::find(entry.second.begin(), entry.second.end(), target_atom) != entry.second.end()) {
                intersects = true;
                break;
            }
        }
        if (intersects) {
            matched_groups.push_back(entry.first);
        }
    }
    std::sort(matched_groups.begin(), matched_groups.end());

    spectrum_widget_->set_selected_groups(matched_groups);

    if (!matched_groups.empty()) {
        const int primary_group = matched_groups.front();
        auto pr = group_to_peak_row_.find(primary_group);
        if (pr != group_to_peak_row_.end()) {
            peak_browser_->value(pr->second);
        }
    }

    if (!nucleus_targets.empty()) {
        const int first_h = nucleus_targets.front();
        auto ar = atom_to_row_.find(first_h);
        if (ar != atom_to_row_.end()) {
            atom_browser_->value(ar->second);
        }
    }

    std::ostringstream status;
    status << "Selected atom " << atom_index;
    if (!nucleus_targets.empty()) {
        status << " | " << nucleus_symbol(active_nucleus_) << ":";
        for (std::size_t i = 0; i < nucleus_targets.size(); ++i) {
            status << nucleus_targets[i];
            if (i + 1 < nucleus_targets.size()) {
                status << ",";
            }
        }
    }
    if (!matched_groups.empty()) {
        status << " | groups:";
        for (std::size_t i = 0; i < matched_groups.size(); ++i) {
            status << matched_groups[i];
            if (i + 1 < matched_groups.size()) {
                status << ",";
            }
        }
    } else {
        status << " (no assigned " << active_nucleus_ << " groups)";
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::apply_nucleus_visuals(const QueuedJob &job, const std::string &nucleus) {
    active_nucleus_ = normalize_nucleus_label(nucleus);

    NucleusResultFiles files;
    auto it = active_nucleus_results_.find(active_nucleus_);
    if (it != active_nucleus_results_.end()) {
        files = it->second;
    } else {
        files.spectrum_csv = job.spectrum_csv;
        files.peaks_csv = job.peaks_csv;
        files.assignments_csv = job.assignments_csv;
    }

    if (!files.spectrum_csv.empty()) {
        if (!spectrum_widget_->load_from_csv(files.spectrum_csv)) {
            status_box_->label("Selected job has no readable spectrum CSV");
        }
    } else {
        spectrum_widget_->set_points({});
        spectrum_widget_->set_peak_markers({});
        spectrum_widget_->set_selected_groups({});
    }

    const auto peak_rows = read_peak_rows(files.peaks_csv);
    group_to_atoms_ = read_assignments(files.assignments_csv);
    apply_reference_peaks(job.config.solvent, active_nucleus_);
    spectrum_widget_->set_nucleus_label(active_nucleus_ + " NMR Spectrum");

    peak_browser_->clear();
    atom_browser_->clear();
    group_to_peak_row_.clear();
    peak_row_to_group_.clear();
    atom_to_row_.clear();
    atom_row_to_atom_.clear();
    atom_row_base_label_.clear();

    std::vector<PeakMarker> markers;
    for (const auto &row : peak_rows) {
        std::ostringstream label;
        label << "G" << row.group_id << " | " << row.center_ppm << " ppm | " << row.multiplicity << " | int " << row.integral;
        peak_browser_->add(label.str().c_str());
        const int browser_row = peak_browser_->size();
        group_to_peak_row_[row.group_id] = browser_row;
        peak_row_to_group_[browser_row] = row.group_id;

        PeakMarker marker;
        marker.group_id = row.group_id;
        marker.center_ppm = row.center_ppm;
        marker.label = row.multiplicity;
        markers.push_back(marker);
    }
    spectrum_widget_->set_peak_markers(std::move(markers));

    std::vector<int> atoms;
    for (const auto &entry : group_to_atoms_) {
        atoms.insert(atoms.end(), entry.second.begin(), entry.second.end());
    }
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    const std::string atom_prefix = nucleus_symbol(active_nucleus_);

    for (int atom_idx : atoms) {
        std::vector<int> groups;
        for (const auto &entry : group_to_atoms_) {
            if (std::find(entry.second.begin(), entry.second.end(), atom_idx) != entry.second.end()) {
                groups.push_back(entry.first);
            }
        }
        std::sort(groups.begin(), groups.end());

        std::ostringstream label;
        label << atom_prefix << atom_idx << " | groups: ";
        for (std::size_t i = 0; i < groups.size(); ++i) {
            label << groups[i];
            if (i + 1 < groups.size()) {
                label << ",";
            }
        }
        atom_browser_->add(label.str().c_str());
        const int row = atom_browser_->size();
        atom_to_row_[atom_idx] = row;
        atom_row_to_atom_[row] = atom_idx;
        atom_row_base_label_[row] = label.str();
    }

    if (!peak_rows.empty()) {
        on_peak_picked(peak_rows.front().group_id);
    } else {
        highlight_hydrogen_rows({});
        if (structure_widget_ != nullptr) {
            structure_widget_->clear_highlight();
        }
    }
}

void AppWindow::load_selected_job_visuals() {
    const int selected = queue_browser_->value();
    if (selected <= 0) {
        return;
    }

    QueuedJob job;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        const std::size_t index = static_cast<std::size_t>(selected - 1);
        if (index >= jobs_.size()) {
            return;
        }
        job = jobs_[index];
    }

    if (nucleus_choice_ != nullptr && nucleus_choice_->size() > 0) {
        nucleus_choice_->value(0);
    }

    auto structure_atoms = read_structure_atoms(job.structure_atoms_csv);
    auto structure_bonds = read_structure_bonds(job.structure_bonds_csv);
    if (structure_widget_ != nullptr) {
        structure_widget_->set_structure(std::move(structure_atoms), std::move(structure_bonds));
    }

    active_nucleus_results_.clear();
    if (!job.spectra_manifest_csv.empty()) {
        active_nucleus_results_ = read_spectra_manifest(job.spectra_manifest_csv);
    }
    if (active_nucleus_results_.empty() && !job.spectrum_csv.empty()) {
        const std::string fallback_nucleus = normalize_nucleus_label(job.config.nucleus);
        active_nucleus_results_[fallback_nucleus] = {
            job.spectrum_csv,
            job.peaks_csv,
            job.assignments_csv,
        };
    }

    if (spectrum_nucleus_choice_ != nullptr) {
        spectrum_nucleus_choice_->clear();

        std::vector<std::string> ordered_labels;
        const std::vector<std::string> preferred = {"1H", "13C", "19F"};
        for (const auto &label : preferred) {
            auto entry = active_nucleus_results_.find(label);
            if (entry != active_nucleus_results_.end() && !entry->second.spectrum_csv.empty()) {
                ordered_labels.push_back(label);
            }
        }
        for (const auto &entry : active_nucleus_results_) {
            if (std::find(ordered_labels.begin(), ordered_labels.end(), entry.first) == ordered_labels.end()
                && !entry.second.spectrum_csv.empty()) {
                ordered_labels.push_back(entry.first);
            }
        }
        if (ordered_labels.empty()) {
            ordered_labels.push_back("1H");
        }

        for (const auto &label : ordered_labels) {
            spectrum_nucleus_choice_->add(label.c_str());
        }

        std::string desired_nucleus = normalize_nucleus_label(active_nucleus_);
        if (std::find(ordered_labels.begin(), ordered_labels.end(), desired_nucleus) == ordered_labels.end()) {
            if (std::find(ordered_labels.begin(), ordered_labels.end(), "1H") != ordered_labels.end()) {
                desired_nucleus = "1H";
            } else {
                desired_nucleus = ordered_labels.front();
            }
        }

        int desired_index = 0;
        for (int i = 0; i < spectrum_nucleus_choice_->size(); ++i) {
            const char *text = spectrum_nucleus_choice_->text(i);
            if (text != nullptr && desired_nucleus == text) {
                desired_index = i;
                break;
            }
        }
        spectrum_nucleus_choice_->value(desired_index);
        if (ordered_labels.size() <= 1) {
            spectrum_nucleus_choice_->deactivate();
        } else {
            spectrum_nucleus_choice_->activate();
        }
        active_nucleus_ = desired_nucleus;
    }

    apply_nucleus_visuals(job, active_nucleus_);
}

void AppWindow::highlight_hydrogen_rows(const std::vector<int> &highlighted_hydrogens) {
    std::unordered_set<int> selected(highlighted_hydrogens.begin(), highlighted_hydrogens.end());
    for (const auto &entry : atom_row_to_atom_) {
        const int row = entry.first;
        const int atom_idx = entry.second;
        auto base_it = atom_row_base_label_.find(row);
        if (base_it == atom_row_base_label_.end()) {
            continue;
        }
        std::string label = base_it->second;
        if (selected.count(atom_idx) > 0) {
            label = "* " + label;
        }
        atom_browser_->text(row, label.c_str());
    }
}

void AppWindow::apply_reference_peaks(const std::string &solvent, const std::string &nucleus) {
    const auto refs = reference_peaks_for_solvent(solvent, nucleus);
    spectrum_widget_->set_reference_peaks(refs);
    spectrum_widget_->set_highlighted_reference(-1);

    if (reference_choice_ == nullptr) {
        return;
    }
    reference_choice_->clear();
    reference_choice_->add("Refs: none");
    for (const auto &ref : refs) {
        std::string item = "Ref: " + ref.label;
        reference_choice_->add(item.c_str());
    }
    reference_choice_->value(0);
}

void AppWindow::refresh_workflow_browser(const QueuedJob *job) {
    if (workflow_browser_ == nullptr) {
        return;
    }

    workflow_browser_->clear();

    if (job == nullptr) {
        if (workflow_progress_widget_ != nullptr) {
            workflow_progress_widget_->set_progress_state(-1, -1, false, false, false, 0.0, "Workflow idle");
        }
        workflow_browser_->add("No active job");
        workflow_browser_->add("Queue a job, then click Run Pending / Run Selected.");
        return;
    }

    std::ostringstream header;
    header << "Job " << job->id << " | " << truncate_text(job->config.job_name, 22)
           << " | " << to_string(job->config.input_format) << " | " << job->config.nucleus;
    workflow_browser_->add(header.str().c_str());

    std::ostringstream current;
    if (!job->progress_message.empty()) {
        current << "Current: " << truncate_text(job->progress_message, 92);
    } else if (job->status == "done") {
        current << "Current: Prediction complete";
    } else if (job->status == "failed") {
        current << "Current: Failed - " << truncate_text(job->message, 72);
    } else if (job->status == "running") {
        current << "Current: Running...";
    } else {
        current << "Current: Waiting to run";
    }
    workflow_browser_->add(current.str().c_str());

    const auto &steps = workflow_steps();
    int active_idx = workflow_stage_index(job->progress_stage);
    if (active_idx < 0) {
        if (job->status == "done") {
            active_idx = static_cast<int>(steps.size() - 1);
        } else if (job->status == "running") {
            active_idx = 0;
        }
    }

    if (workflow_progress_widget_ != nullptr) {
        if (job->status == "done") {
            workflow_progress_widget_->set_progress_state(
                static_cast<int>(steps.size() - 1),
                -1,
                false,
                true,
                false,
                1.0,
                "Complete: " + truncate_text(job->config.job_name, 40));
        } else if (job->status == "failed") {
            workflow_progress_widget_->set_progress_state(
                active_idx >= 0 ? active_idx : 0,
                active_idx >= 0 ? active_idx : 0,
                false,
                false,
                true,
                job->progress_fraction,
                "Failed: " + truncate_text(job->message, 64));
        } else if (job->status == "running") {
            workflow_progress_widget_->set_progress_state(
                active_idx,
                -1,
                true,
                false,
                false,
                job->progress_fraction,
                job->progress_message.empty() ? "Running workflow" : job->progress_message);
        } else {
            workflow_progress_widget_->set_progress_state(
                -1,
                -1,
                false,
                false,
                false,
                0.0,
                "Queued: waiting to run");
        }
    }

    for (std::size_t i = 0; i < steps.size(); ++i) {
        std::string icon = "[ ]";
        if (job->status == "done") {
            icon = "[x]";
        } else if (job->status == "failed") {
            if (active_idx >= 0 && static_cast<int>(i) < active_idx) {
                icon = "[x]";
            } else if (active_idx >= 0 && static_cast<int>(i) == active_idx) {
                icon = "[!]";
            }
        } else if (job->status == "running") {
            if (active_idx >= 0 && static_cast<int>(i) < active_idx) {
                icon = "[x]";
            } else if (active_idx >= 0 && static_cast<int>(i) == active_idx) {
                icon = "[>]";
            }
        }

        std::ostringstream line;
        line << icon << " " << steps[i].label << " - " << steps[i].method;
        workflow_browser_->add(line.str().c_str());
    }

    if (job->status == "running" && job->progress_fraction > 0.0) {
        const int pct = static_cast<int>(job->progress_fraction * 100.0 + 0.5);
        std::ostringstream frac;
        frac << "Overall progress: " << pct << "%";
        workflow_browser_->add(frac.str().c_str());
    }
}

void AppWindow::run_worker_loop() {
    const std::size_t no_index = static_cast<std::size_t>(-1);

    while (true) {
        if (cancel_requested_) {
            break;
        }

        std::size_t next_index = no_index;
        {
            std::lock_guard<std::mutex> lock(jobs_mutex_);
            if (run_scope_ == RunScope::SelectedOnly) {
                if (run_done_ == 0 && selected_run_index_ < jobs_.size() && jobs_[selected_run_index_].status == "pending") {
                    next_index = selected_run_index_;
                    jobs_[next_index].status = "running";
                    const std::string nucleus_mode = (normalize_nucleus_label(jobs_[next_index].config.nucleus) == "1H"
                                                      && jobs_[next_index].config.nucleus != "auto")
                                                         ? "1H model"
                                                         : "multi-nucleus (1H/13C/19F as present)";
                    jobs_[next_index].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> "
                                              + nucleus_mode;
                    jobs_[next_index].progress_stage = "launch";
                    jobs_[next_index].progress_message = "Launching backend process";
                    jobs_[next_index].progress_fraction = 0.0;
                }
            } else {
                for (std::size_t i = 0; i < jobs_.size(); ++i) {
                    if (jobs_[i].status == "pending") {
                        next_index = i;
                        jobs_[i].status = "running";
                        const std::string nucleus_mode = (normalize_nucleus_label(jobs_[i].config.nucleus) == "1H"
                                                          && jobs_[i].config.nucleus != "auto")
                                                             ? "1H model"
                                                             : "multi-nucleus (1H/13C/19F as present)";
                        jobs_[i].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> "
                                         + nucleus_mode;
                        jobs_[i].progress_stage = "launch";
                        jobs_[i].progress_message = "Launching backend process";
                        jobs_[i].progress_fraction = 0.0;
                        break;
                    }
                }
            }
            if (next_index != no_index) {
                active_job_index_ = static_cast<int>(next_index);
            }
        }

        if (next_index == no_index) {
            break;
        }

        Fl::awake(on_worker_awake, this);

        JobConfig config;
        {
            std::lock_guard<std::mutex> lock(jobs_mutex_);
            config = jobs_[next_index].config;
        }

        const JobOutputs result = pipeline_.run(config, [this, next_index](const PipelineProgress &progress) {
            {
                std::lock_guard<std::mutex> lock(jobs_mutex_);
                if (next_index < jobs_.size() && jobs_[next_index].status == "running") {
                    std::ostringstream msg;
                    if (!progress.message.empty()) {
                        msg << progress.message;
                    } else if (!progress.stage.empty()) {
                        msg << progress.stage;
                    } else {
                        msg << "running";
                    }
                    if (progress.fraction > 0.0) {
                        const int pct = static_cast<int>(progress.fraction * 100.0 + 0.5);
                        msg << " (" << pct << "%)";
                    }
                    jobs_[next_index].message = msg.str();
                    jobs_[next_index].progress_stage = progress.stage;
                    jobs_[next_index].progress_message = msg.str();
                    jobs_[next_index].progress_fraction = progress.fraction;
                }
            }
            Fl::awake(on_worker_awake, this);
        });

        {
            std::lock_guard<std::mutex> lock(jobs_mutex_);
            if (next_index < jobs_.size()) {
                auto &job = jobs_[next_index];
                job.status = result.status == "ok" ? "done" : "failed";
                job.output_dir = result.output_dir;
                job.spectrum_csv = result.spectrum_csv;
                job.peaks_csv = result.peaks_csv;
                job.assignments_csv = result.assignments_csv;
                job.spectra_manifest_csv = result.spectra_manifest_csv;
                job.structure_svg = result.structure_svg;
                job.structure_atoms_csv = result.structure_atoms_csv;
                job.structure_bonds_csv = result.structure_bonds_csv;
                job.structure_xyz = result.structure_xyz;
                if (result.status == "ok") {
                    job.message = "done";
                    job.progress_stage = "done";
                    job.progress_message = "Prediction complete";
                    job.progress_fraction = 1.0;
                } else {
                    job.message = result.message;
                    if (job.progress_stage.empty()) {
                        job.progress_stage = "failed";
                    }
                    if (job.progress_message.empty()) {
                        job.progress_message = result.message;
                    }
                }
            }
            ++run_done_;
            active_job_index_ = -1;
        }

        Fl::awake(on_worker_awake, this);

        if (run_scope_ == RunScope::SelectedOnly) {
            break;
        }
    }

    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        active_job_index_ = -1;
    }
    worker_running_ = false;
    Fl::awake(on_worker_awake, this);
}

void AppWindow::refresh_queue_browser() {
    std::lock_guard<std::mutex> lock(jobs_mutex_);

    const int selected = queue_browser_->value();
    queue_browser_->clear();

    for (std::size_t i = 0; i < jobs_.size(); ++i) {
        const bool active = (static_cast<int>(i) == active_job_index_);
        queue_browser_->add(format_label_for(jobs_[i], active).c_str());
    }

    if (selected > 0 && selected <= queue_browser_->size()) {
        queue_browser_->value(selected);
    }

    int pending = 0;
    int running = 0;
    int done = 0;
    int failed = 0;
    for (const auto &job : jobs_) {
        if (job.status == "pending") {
            ++pending;
        } else if (job.status == "running") {
            ++running;
        } else if (job.status == "done") {
            ++done;
        } else if (job.status == "failed") {
            ++failed;
        }
    }

    std::ostringstream status;
    if (worker_running_) {
        status << "RUN ";
        status << (run_scope_ == RunScope::SelectedOnly ? "sel" : "all");
        status << " " << run_done_ << "/" << std::max(run_total_, 1);
        if (active_job_index_ >= 0 && active_job_index_ < static_cast<int>(jobs_.size())) {
            status << " q" << (active_job_index_ + 1);
            const auto &active = jobs_[static_cast<std::size_t>(active_job_index_)];
            if (!active.message.empty()) {
                status << " " << truncate_text(active.message, 46);
            }
        }
        status << " p" << pending << " d" << done << " f" << failed;
    } else {
        status << "p" << pending << " r" << running << " d" << done << " f" << failed;
    }
    status_box_->copy_label(status.str().c_str());

    const QueuedJob *workflow_job = nullptr;
    if (active_job_index_ >= 0 && active_job_index_ < static_cast<int>(jobs_.size())) {
        workflow_job = &jobs_[static_cast<std::size_t>(active_job_index_)];
    }
    if (workflow_job == nullptr && selected > 0 && selected <= queue_browser_->size()) {
        const std::size_t idx = static_cast<std::size_t>(selected - 1);
        if (idx < jobs_.size()) {
            workflow_job = &jobs_[idx];
        }
    }
    if (workflow_job == nullptr && !jobs_.empty()) {
        workflow_job = &jobs_.back();
    }
    refresh_workflow_browser(workflow_job);
}

} // namespace easynmr
