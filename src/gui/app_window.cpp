#include "gui/app_window.h"

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Image_Surface.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <FL/Fl_PNG_Image.H>
#include <FL/Fl_Pixmap.H>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
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
    double j_hz = 0.0;
    double integral = 0.0;
};

struct WorkflowStep {
    std::string label;
    std::string method;
};

struct ExampleCase {
    std::string difficulty;
    std::string case_name;
    std::string workflow;
    std::string target_product;
    std::string nucleus_arg;
    std::string smiles;
    std::string computed_reference_csv;
    std::string experimental_overlay_file;
    std::string experimental_format_hint;
    std::string notes;
    std::string csv_row;
};

std::string normalize_nucleus_label(const std::string &raw);
bool is_proton_nucleus(const std::string &nucleus);
std::string nucleus_symbol(const std::string &nucleus);
bool is_known_nmr_nucleus(const std::string &label);
std::string spectrum_title_for_label(const std::string &label);
std::string workflow_display_name(WorkflowKind kind);
std::string make_overlay_label_from_path(const std::string &path);
std::string make_unique_overlay_key(
    const std::map<std::string, std::vector<SpectrumPoint>> &overlays,
    const std::string &base_key);
InputFormat detect_input_format_auto(const std::string &value);
std::vector<std::string> split_smiles_inputs(const std::string &value);
std::string trim_copy(const std::string &value);
std::vector<std::string> parse_csv_line(const std::string &line);

std::string uppercase_copy(const std::string &value) {
    std::string out = trim_copy(value);
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return out;
}

std::string find_examples_case_csv() {
    namespace fs = std::filesystem;

    std::error_code ec;
    fs::path probe = fs::current_path(ec);
    if (ec) {
        return {};
    }

    const fs::path relative("tests/spectra_comparison_cases.csv");
    for (int depth = 0; depth < 8; ++depth) {
        const fs::path candidate = probe / relative;
        ec.clear();
        if (fs::exists(candidate, ec) && !ec) {
            return candidate.lexically_normal().string();
        }
        if (!probe.has_parent_path()) {
            break;
        }
        const fs::path parent = probe.parent_path();
        if (parent == probe) {
            break;
        }
        probe = parent;
    }
    return {};
}

std::string resolve_example_asset_path(const std::string &value) {
    namespace fs = std::filesystem;
    const std::string clean = trim_copy(value);
    if (clean.empty()) {
        return {};
    }

    fs::path p(clean);
    if (p.is_absolute()) {
        return p.lexically_normal().string();
    }

    std::error_code ec;
    fs::path cwd = fs::current_path(ec);
    if (!ec) {
        fs::path direct = cwd / p;
        ec.clear();
        if (fs::exists(direct, ec) && !ec) {
            return direct.lexically_normal().string();
        }
    }

    const std::string cases_csv = find_examples_case_csv();
    if (!cases_csv.empty()) {
        const fs::path repo_root = fs::path(cases_csv).parent_path().parent_path();
        fs::path candidate = repo_root / p;
        ec.clear();
        if (fs::exists(candidate, ec) && !ec) {
            return candidate.lexically_normal().string();
        }
    }

    return p.lexically_normal().string();
}

std::string lowercase_copy(const std::string &value) {
    std::string out = trim_copy(value);
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

std::vector<std::string> split_tokens_alnum(const std::string &value) {
    std::vector<std::string> tokens;
    std::string token;
    for (char c : value) {
        if (std::isalnum(static_cast<unsigned char>(c)) != 0) {
            token.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        } else if (!token.empty()) {
            tokens.push_back(token);
            token.clear();
        }
    }
    if (!token.empty()) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string title_case_tokens(const std::vector<std::string> &tokens) {
    if (tokens.empty()) {
        return {};
    }
    std::ostringstream oss;
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        std::string token = tokens[i];
        if (!token.empty()) {
            token[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
        }
        oss << token;
        if (i + 1 < tokens.size()) {
            oss << " ";
        }
    }
    return oss.str();
}

std::string infer_example_bundle_name(const ExampleCase &example) {
    const std::unordered_set<std::string> ignore = {
        "easy",
        "medium",
        "hard",
        "nmr",
        "cd",
        "auto",
        "h",
        "c",
        "f",
        "p",
        "1h",
        "13c",
        "19f",
        "31p",
    };

    const auto all_tokens = split_tokens_alnum(example.case_name);
    std::vector<std::string> kept;
    for (const auto &token : all_tokens) {
        if (ignore.count(token) > 0) {
            continue;
        }
        kept.push_back(token);
    }
    if (kept.empty() && !example.smiles.empty()) {
        return trim_copy(example.smiles);
    }
    return title_case_tokens(kept);
}

std::string pick_default_product(const std::set<std::string> &products) {
    if (products.count("1H") > 0) {
        return "1H";
    }
    if (products.count("13C") > 0) {
        return "13C";
    }
    if (products.count("19F") > 0) {
        return "19F";
    }
    if (products.count("31P") > 0) {
        return "31P";
    }
    if (products.count("CD") > 0) {
        return "CD";
    }
    if (!products.empty()) {
        return *products.begin();
    }
    return "1H";
}

std::string sanitize_filename_component(const std::string &value) {
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        if (std::isalnum(static_cast<unsigned char>(c)) != 0) {
            out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        } else if (c == '_' || c == '-') {
            out.push_back(c);
        } else {
            out.push_back('_');
        }
    }
    if (out.empty()) {
        out = "example";
    }
    return out;
}

std::string csv_escape(const std::string &value) {
    if (value.find(',') == std::string::npos && value.find('"') == std::string::npos) {
        return value;
    }
    std::string escaped;
    escaped.reserve(value.size() + 4);
    escaped.push_back('"');
    for (char c : value) {
        if (c == '"') {
            escaped.push_back('"');
        }
        escaped.push_back(c);
    }
    escaped.push_back('"');
    return escaped;
}

std::string write_example_manifest(
    const std::string &bundle_name,
    const std::map<std::string, std::string> &computed_by_product) {
    namespace fs = std::filesystem;
    if (computed_by_product.empty()) {
        return {};
    }

    std::error_code ec;
    fs::path out_dir = fs::temp_directory_path(ec);
    if (ec) {
        return {};
    }
    out_dir /= "easynmr_example_manifests";
    fs::create_directories(out_dir, ec);
    if (ec) {
        return {};
    }

    const auto now = std::chrono::system_clock::now();
    const auto stamp = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    const fs::path path = out_dir / (sanitize_filename_component(bundle_name) + "_" + std::to_string(stamp) + ".csv");

    std::ofstream out(path);
    if (!out) {
        return {};
    }
    out << "label,spectrum_csv,peaks_csv,assignments_csv\n";
    for (const auto &entry : computed_by_product) {
        out << csv_escape(entry.first) << ","
            << csv_escape(entry.second) << ",,\n";
    }
    if (!out) {
        return {};
    }
    return path.string();
}

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

bool write_ppm_rgb(
    const std::string &path,
    const unsigned char *rgb,
    int width,
    int height,
    std::string *error_message) {
    if (rgb == nullptr || width <= 0 || height <= 0) {
        if (error_message != nullptr) {
            *error_message = "Invalid image buffer for export.";
        }
        return false;
    }

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error_message != nullptr) {
            *error_message = "Could not open output file for writing.";
        }
        return false;
    }

    out << "P6\n" << width << " " << height << "\n255\n";
    out.write(reinterpret_cast<const char *>(rgb), static_cast<std::streamsize>(width * height * 3));
    if (!out) {
        if (error_message != nullptr) {
            *error_message = "Failed while writing image data.";
        }
        return false;
    }
    return true;
}

void trim_white_border(
    std::vector<unsigned char> *rgb,
    int *width,
    int *height) {
    if (rgb == nullptr || width == nullptr || height == nullptr) {
        return;
    }
    if (*width <= 0 || *height <= 0) {
        return;
    }
    const int src_w = *width;
    const int src_h = *height;
    if (static_cast<int>(rgb->size()) < src_w * src_h * 3) {
        return;
    }

    auto is_background = [&](int x, int y) {
        const std::size_t idx = static_cast<std::size_t>((y * src_w + x) * 3);
        const unsigned char r = (*rgb)[idx + 0];
        const unsigned char g = (*rgb)[idx + 1];
        const unsigned char b = (*rgb)[idx + 2];
        // Only treat near-pure white as padding so we keep the light plot colors.
        return r >= 254 && g >= 254 && b >= 254;
    };

    int min_x = src_w;
    int min_y = src_h;
    int max_x = -1;
    int max_y = -1;
    for (int y = 0; y < src_h; ++y) {
        for (int x = 0; x < src_w; ++x) {
            if (is_background(x, y)) {
                continue;
            }
            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
        }
    }

    if (max_x < min_x || max_y < min_y) {
        return;
    }

    const int kPad = 2;
    min_x = std::max(0, min_x - kPad);
    min_y = std::max(0, min_y - kPad);
    max_x = std::min(src_w - 1, max_x + kPad);
    max_y = std::min(src_h - 1, max_y + kPad);

    const int dst_w = max_x - min_x + 1;
    const int dst_h = max_y - min_y + 1;
    if (dst_w >= src_w && dst_h >= src_h) {
        return;
    }

    std::vector<unsigned char> cropped;
    cropped.resize(static_cast<std::size_t>(dst_w * dst_h * 3));
    for (int y = 0; y < dst_h; ++y) {
        const unsigned char *src_row = rgb->data() + static_cast<std::size_t>(((y + min_y) * src_w + min_x) * 3);
        unsigned char *dst_row = cropped.data() + static_cast<std::size_t>(y * dst_w * 3);
        std::copy(src_row, src_row + static_cast<std::size_t>(dst_w * 3), dst_row);
    }

    *rgb = std::move(cropped);
    *width = dst_w;
    *height = dst_h;
}

bool capture_widget_snapshot_rgb(
    Fl_Widget *widget,
    std::vector<unsigned char> *rgb,
    int *out_w,
    int *out_h,
    std::string *error_message) {
    if (rgb == nullptr || out_w == nullptr || out_h == nullptr) {
        if (error_message != nullptr) {
            *error_message = "Internal export buffer is not available.";
        }
        return false;
    }
    if (widget == nullptr || widget->w() <= 0 || widget->h() <= 0) {
        if (error_message != nullptr) {
            *error_message = "Spectrum widget is not ready for export.";
        }
        return false;
    }

    Fl_Image_Surface surface(widget->w(), widget->h());
    Fl_Surface_Device::push_current(&surface);
    fl_color(FL_WHITE);
    fl_rectf(0, 0, widget->w(), widget->h());
    surface.draw(widget, -widget->x(), -widget->y());
    Fl_Surface_Device::pop_current();

    Fl_RGB_Image *img = surface.image();
    if (img == nullptr) {
        if (error_message != nullptr) {
            *error_message = "Could not capture spectrum widget image.";
        }
        return false;
    }

    const char *const *raw_data = img->data();
    if (raw_data == nullptr || img->count() < 1 || raw_data[0] == nullptr) {
        img->release();
        if (error_message != nullptr) {
            *error_message = "Captured image buffer is empty.";
        }
        return false;
    }

    const int width = img->data_w();
    const int height = img->data_h();
    const int depth = img->d();
    const int line_stride = (img->ld() > 0) ? img->ld() : width * depth;
    if (width <= 0 || height <= 0 || depth < 3 || line_stride < width * depth) {
        img->release();
        if (error_message != nullptr) {
            *error_message = "Unsupported exported image format.";
        }
        return false;
    }

    const unsigned char *src = reinterpret_cast<const unsigned char *>(raw_data[0]);
    rgb->resize(static_cast<std::size_t>(width * height * 3));
    for (int y = 0; y < height; ++y) {
        const unsigned char *src_row = src + static_cast<std::size_t>(y * line_stride);
        unsigned char *dst_row = rgb->data() + static_cast<std::size_t>(y * width * 3);
        for (int x = 0; x < width; ++x) {
            const int src_offset = x * depth;
            const int dst_offset = x * 3;
            dst_row[dst_offset + 0] = src_row[src_offset + 0];
            dst_row[dst_offset + 1] = src_row[src_offset + 1];
            dst_row[dst_offset + 2] = src_row[src_offset + 2];
        }
    }

    img->release();
    *out_w = width;
    *out_h = height;
    trim_white_border(rgb, out_w, out_h);
    return true;
}

bool export_widget_snapshot_ppm(
    Fl_Widget *widget,
    const std::string &output_path,
    std::string *error_message) {
    std::vector<unsigned char> rgb;
    int out_w = 0;
    int out_h = 0;
    if (!capture_widget_snapshot_rgb(widget, &rgb, &out_w, &out_h, error_message)) {
        return false;
    }
    return write_ppm_rgb(output_path, rgb.data(), out_w, out_h, error_message);
}

bool export_widget_snapshot_png(
    Fl_Widget *widget,
    const std::string &output_path,
    std::string *error_message) {
    std::vector<unsigned char> rgb;
    int out_w = 0;
    int out_h = 0;
    if (!capture_widget_snapshot_rgb(widget, &rgb, &out_w, &out_h, error_message)) {
        return false;
    }
    std::error_code ec;
    std::filesystem::remove(output_path, ec);
    (void)fl_write_png(output_path.c_str(), rgb.data(), out_w, out_h, 3, out_w * 3);
    ec.clear();
    if (!std::filesystem::exists(output_path, ec) || ec) {
        if (error_message != nullptr) {
            *error_message = "Could not write PNG output.";
        }
        return false;
    }
    const auto written_size = std::filesystem::file_size(output_path, ec);
    if (ec || written_size <= 0) {
        if (error_message != nullptr) {
            *error_message = "PNG output file is empty.";
        }
        return false;
    }
    return true;
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

const char *kRoleBlueDotXpm[] = {
    "9 9 2 1",
    "  c None",
    ". c #6CA6E8",
    "         ",
    "         ",
    "   ...   ",
    "  .....  ",
    "  .....  ",
    "  .....  ",
    "   ...   ",
    "         ",
    "         ",
};

const char *kRoleGreenDotXpm[] = {
    "9 9 2 1",
    "  c None",
    ". c #D7A97E",
    "         ",
    "         ",
    "   ...   ",
    "  .....  ",
    "  .....  ",
    "  .....  ",
    "   ...   ",
    "         ",
    "         ",
};

Fl_Pixmap kRoleBlueDot(kRoleBlueDotXpm);
Fl_Pixmap kRoleGreenDot(kRoleGreenDotXpm);

std::string format_label_for(
    const QueuedJob &job,
    bool selected,
    bool is_primary_role,
    bool is_compare_role) {
    std::string state_icon = "WAIT";
    if (job.status == "running") {
        state_icon = "RUN";
    } else if (job.status == "done") {
        state_icon = "OK";
    } else if (job.status == "failed") {
        state_icon = "FAIL";
    }

    std::string role = ".";
    if (is_primary_role) {
        role = "A";
    } else if (is_compare_role) {
        role = "B";
    }

    std::string run_or_select = " ";
    if (selected) {
        run_or_select = ">";
    } else if (job.status == "running") {
        run_or_select = "*";
    }

    std::string mode = workflow_display_name(job.config.workflow_kind);
    if (job.config.workflow_kind == WorkflowKind::Nmr) {
        mode += ":" + normalize_nucleus_label(job.config.nucleus);
    }

    std::string tail;
    if (job.status == "running" && !job.message.empty()) {
        tail = truncate_text(job.message, 22);
    } else if (job.status == "failed" && !job.message.empty()) {
        tail = truncate_text(job.message, 22);
    } else {
        tail = to_string(job.config.input_format);
    }

    std::ostringstream oss;
    oss << run_or_select
        << "\t"
        << role
        << "\t"
        << state_icon
        << "\t"
        << truncate_text(job.config.job_name, 20)
        << "\t"
        << job.id
        << "\t"
        << mode
        << "\t"
        << tail;
    return oss.str();
}

std::vector<ExampleCase> load_example_cases(const std::string &path) {
    std::vector<ExampleCase> out;
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
        if (trim_copy(line).empty()) {
            continue;
        }
        const auto fields = parse_csv_line(line);
        if (fields.size() < 10) {
            continue;
        }
        ExampleCase c;
        c.difficulty = trim_copy(fields[0]);
        c.case_name = trim_copy(fields[1]);
        c.workflow = trim_copy(fields[2]);
        c.target_product = trim_copy(fields[3]);
        c.nucleus_arg = trim_copy(fields[4]);
        c.smiles = trim_copy(fields[5]);
        c.computed_reference_csv = trim_copy(fields[6]);
        c.experimental_overlay_file = trim_copy(fields[7]);
        c.experimental_format_hint = trim_copy(fields[8]);
        c.notes = trim_copy(fields[9]);
        c.csv_row = line;
        if (!c.case_name.empty() && !c.computed_reference_csv.empty()) {
            out.push_back(std::move(c));
        }
    }
    return out;
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

std::vector<std::string> split_whitespace_tokens(const std::string &line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

bool looks_like_xyz_block(const std::vector<std::string> &lines) {
    if (lines.size() < 3) {
        return false;
    }

    int atom_count = 0;
    if (!try_parse_int(lines[0], atom_count) || atom_count <= 0) {
        return false;
    }
    if (lines.size() < static_cast<std::size_t>(atom_count + 2)) {
        return false;
    }

    for (int i = 0; i < atom_count; ++i) {
        const std::size_t line_idx = static_cast<std::size_t>(i + 2);
        if (line_idx >= lines.size()) {
            return false;
        }
        const auto tokens = split_whitespace_tokens(lines[line_idx]);
        if (tokens.size() < 4) {
            return false;
        }
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!try_parse_double(tokens[1], x) || !try_parse_double(tokens[2], y) || !try_parse_double(tokens[3], z)) {
            return false;
        }
    }

    return true;
}

InputFormat detect_input_format_auto(const std::string &value) {
    const std::string clean = trim_copy(value);
    if (clean.empty()) {
        return InputFormat::Unknown;
    }

    std::string lowered = clean;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (lowered.find("$$$$") != std::string::npos) {
        return InputFormat::Sdf;
    }
    if (lowered.find("v2000") != std::string::npos
        || lowered.find("v3000") != std::string::npos
        || lowered.find("m  end") != std::string::npos) {
        return InputFormat::Mol;
    }

    std::vector<std::string> nonempty_lines;
    {
        std::istringstream iss(clean);
        std::string line;
        while (std::getline(iss, line)) {
            const std::string trimmed_line = trim_copy(line);
            if (!trimmed_line.empty()) {
                nonempty_lines.push_back(trimmed_line);
            }
        }
    }

    if (looks_like_xyz_block(nonempty_lines)) {
        return InputFormat::Xyz;
    }

    return InputFormat::Smiles;
}

std::vector<std::string> split_smiles_inputs(const std::string &value) {
    std::vector<std::string> entries;
    std::string current;
    int bracket_depth = 0;
    int brace_depth = 0;

    const auto flush_current = [&]() {
        const std::string trimmed = trim_copy(current);
        if (!trimmed.empty()) {
            entries.push_back(trimmed);
        }
        current.clear();
    };

    for (char c : value) {
        if (c == '[') {
            ++bracket_depth;
        } else if (c == ']' && bracket_depth > 0) {
            --bracket_depth;
        } else if (c == '{') {
            ++brace_depth;
        } else if (c == '}' && brace_depth > 0) {
            --brace_depth;
        }

        const bool is_separator =
            (c == ',' || c == '\n' || c == ';' || c == '\r') && bracket_depth == 0 && brace_depth == 0;
        if (is_separator) {
            flush_current();
            continue;
        }
        current.push_back(c);
    }
    flush_current();

    return entries;
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
        try_parse_double(fields[4], row.j_hz);
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
        if (fields.size() >= 4) {
            bond.stereo_style = trim_copy(fields[3]);
            std::transform(
                bond.stereo_style.begin(),
                bond.stereo_style.end(),
                bond.stereo_style.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (bond.stereo_style != "wedge" && bond.stereo_style != "dash") {
                bond.stereo_style = "none";
            }
        }
        if (fields.size() >= 5) {
            int stereo_from_atom = 0;
            if (try_parse_int(fields[4], stereo_from_atom)) {
                bond.stereo_from_atom = stereo_from_atom;
            }
        }
        bonds.push_back(std::move(bond));
    }

    return bonds;
}

std::map<std::string, SpectralProductFiles> read_spectra_manifest(const std::string &path) {
    std::map<std::string, SpectralProductFiles> out;
    auto manifest = load_spectral_products_manifest(path);
    for (auto &files : manifest) {
        const std::string normalized = normalize_nucleus_label(files.label);
        if (!normalized.empty()) {
            files.label = normalized;
        }
        out[files.label] = std::move(files);
    }
    return out;
}

std::map<std::string, SpectralProductFiles> spectral_products_for_job(const QueuedJob &job) {
    std::map<std::string, SpectralProductFiles> products;
    if (!job.spectra_manifest_csv.empty()) {
        products = read_spectra_manifest(job.spectra_manifest_csv);
    }
    if (!products.empty()) {
        return products;
    }

    if (job.spectrum_csv.empty()) {
        return products;
    }

    const std::string fallback_nucleus =
        (job.config.workflow_kind == WorkflowKind::Cd) ? "CD" : normalize_nucleus_label(job.config.nucleus);
    SpectralProductFiles fallback;
    fallback.label = fallback_nucleus.empty() ? "1H" : fallback_nucleus;
    fallback.spectrum_csv = job.spectrum_csv;
    fallback.peaks_csv = job.peaks_csv;
    fallback.assignments_csv = job.assignments_csv;
    products[fallback.label] = std::move(fallback);
    return products;
}

int browser_line_at_mouse(const Fl_Hold_Browser *browser, int mouse_y) {
    if (browser == nullptr || browser->size() <= 0) {
        return 0;
    }
    const int bx = browser->x();
    const int by = browser->y();
    const int bh = browser->h();
    if (mouse_y < by || mouse_y > by + bh) {
        return 0;
    }

    const int row_height = std::max(14, static_cast<int>(browser->textsize()) + 6);
    if (row_height <= 0) {
        return 0;
    }
    const int y_in_list = mouse_y - by + browser->vposition();
    const int line = y_in_list / row_height + 1;
    if (line < 1 || line > browser->size()) {
        return 0;
    }
    return line;
}

std::vector<ReferencePeak> reference_peaks_for_solvent(const std::string &solvent_raw, const std::string &nucleus_raw) {
    std::string solvent = solvent_raw;
    std::transform(solvent.begin(), solvent.end(), solvent.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    std::string nucleus = nucleus_raw;
    std::transform(nucleus.begin(), nucleus.end(), nucleus.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (nucleus == "h1" || nucleus == "h") {
        nucleus = "1h";
    } else if (nucleus == "c13" || nucleus == "c") {
        nucleus = "13c";
    } else if (nucleus == "f19" || nucleus == "f") {
        nucleus = "19f";
    } else if (nucleus == "p31" || nucleus == "p") {
        nucleus = "31p";
    }

    if (nucleus != "1h" && nucleus != "13c" && nucleus != "19f" && nucleus != "31p") {
        return {};
    }

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
    } else if (nucleus == "31p") {
        peaks.push_back({0.00, "H3PO4 0.00"});
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
        {"Spectral parameter estimation", "Predict NMR or CD observables from selected conformers"},
        {"Spectrum simulation", "Apply product-specific broadening model"},
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
    if (stage == "nmr_parameter_estimation" || stage == "cd_parameter_estimation") {
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
    std::string nucleus = trim_copy(raw);
    if (nucleus.empty()) {
        return {};
    }
    std::transform(nucleus.begin(), nucleus.end(), nucleus.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (nucleus == "auto" || nucleus == "all") {
        return "1H";
    }
    if (nucleus == "1h" || nucleus == "h1" || nucleus == "h") {
        return "1H";
    }
    if (nucleus == "13c" || nucleus == "c13" || nucleus == "c") {
        return "13C";
    }
    if (nucleus == "19f" || nucleus == "f19" || nucleus == "f") {
        return "19F";
    }
    if (nucleus == "31p" || nucleus == "p31" || nucleus == "p") {
        return "31P";
    }
    return trim_copy(raw);
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
    if (normalized == "31P") {
        return "P";
    }
    if (normalized == "1H") {
        return "H";
    }
    return "Nuc";
}

bool is_known_nmr_nucleus(const std::string &label) {
    const std::string normalized = normalize_nucleus_label(label);
    return normalized == "1H" || normalized == "13C" || normalized == "19F" || normalized == "31P";
}

std::string spectrum_title_for_label(const std::string &label) {
    const std::string normalized = normalize_nucleus_label(label);
    if (is_known_nmr_nucleus(normalized)) {
        return normalized + " NMR Spectrum";
    }
    if (!normalized.empty()) {
        return normalized + " Spectrum";
    }
    return "Spectrum";
}

std::string workflow_display_name(WorkflowKind kind) {
    if (kind == WorkflowKind::All) {
        return "ALL";
    }
    if (kind == WorkflowKind::Cd) {
        return "CD";
    }
    if (kind == WorkflowKind::Compare) {
        return "COMPARE";
    }
    return "NMR";
}

std::string make_overlay_label_from_path(const std::string &path_text) {
    namespace fs = std::filesystem;
    std::string label;
    try {
        const fs::path path(path_text);
        label = path.stem().string();
    } catch (...) {
        label.clear();
    }
    label = trim_copy(label);
    if (label.empty()) {
        label = "experimental";
    }
    return label;
}

std::string make_unique_overlay_key(
    const std::map<std::string, std::vector<SpectrumPoint>> &overlays,
    const std::string &base_key) {
    std::string candidate = trim_copy(base_key);
    if (candidate.empty()) {
        candidate = "experimental";
    }
    if (overlays.find(candidate) == overlays.end()) {
        return candidate;
    }
    for (int idx = 2; idx <= 9999; ++idx) {
        const std::string with_suffix = candidate + " (" + std::to_string(idx) + ")";
        if (overlays.find(with_suffix) == overlays.end()) {
            return with_suffix;
        }
    }
    return candidate + " (copy)";
}

std::string experimental_import_failure_status(const ExperimentalSpectrumLoadResult &loaded) {
    std::string hint = " Use a 2-column text/CSV file (x ppm, y intensity).";
    if (loaded.detected_format == "bruker_raw_directory") {
        hint = " Select an exported ASCII/CSV file from TopSpin or MNova, not the raw Bruker folder.";
    } else if (loaded.detected_format == "mnova_project_file") {
        hint = " Export from MNova as text/CSV (x,y), then import that exported file.";
    }
    return "Experimental import failed: " + loaded.error_message + hint;
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
    auto *top_title = new Fl_Box(12, 6, 420, 24, "easySpectra  |  Local Spectra Predictor");
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

    auto *workflow_label = new Fl_Box(panel_x + 10, panel_y + 82, 100, 16, "Workflow");
    workflow_label->box(FL_NO_BOX);
    workflow_label->labelsize(12);
    workflow_label->labelcolor(ui(86, 97, 112));
    workflow_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    workflow_choice_ = new Fl_Choice(panel_x + 10, panel_y + 98, 100, 26);
    workflow_choice_->add("all");
    workflow_choice_->add("nmr");
    workflow_choice_->add("cd");
    workflow_choice_->add("compare");
    workflow_choice_->value(0);
    workflow_choice_->box(FL_DOWN_BOX);
    workflow_choice_->color(ui(255, 255, 255));
    workflow_choice_->callback([](Fl_Widget *w, void *ud) {
        auto *self = static_cast<AppWindow *>(ud);
        const Fl_Choice *ch = static_cast<Fl_Choice *>(w);
        const char *sel = ch->text(ch->value());
        const bool is_compare = sel != nullptr && std::string(sel) == "compare";
        if (self->compare_input_label_ != nullptr) {
            if (is_compare) {
                self->compare_input_label_->show();
            } else {
                self->compare_input_label_->hide();
            }
        }
        if (self->compare_input_box_ != nullptr) {
            if (is_compare) {
                self->compare_input_box_->show();
            } else {
                self->compare_input_box_->hide();
            }
        }
        on_preview_cb(w, ud);
    }, this);

    auto *solvent_label = new Fl_Box(panel_x + 116, panel_y + 82, 208, 16, "Solvent");
    solvent_label->box(FL_NO_BOX);
    solvent_label->labelsize(12);
    solvent_label->labelcolor(ui(86, 97, 112));
    solvent_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    solvent_choice_ = new Fl_Choice(panel_x + 116, panel_y + 98, 208, 26);
    solvent_choice_->add("cdcl3");
    solvent_choice_->add("dmso");
    solvent_choice_->add("h2o");
    solvent_choice_->value(0);
    solvent_choice_->box(FL_DOWN_BOX);
    solvent_choice_->color(ui(255, 255, 255));
    solvent_choice_->callback(on_preview_cb, this);

    auto *input_label = new Fl_Box(panel_x + 10, panel_y + 128, panel_w - 20, 16, "SMILES / structure text");
    input_label->box(FL_NO_BOX);
    input_label->labelsize(12);
    input_label->labelcolor(ui(86, 97, 112));
    input_label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    input_box_ = new ColoredInputEditor(panel_x + 10, panel_y + 144, panel_w - 20, 64);
    input_box_->value("CCO");
    input_box_->set_syntax_mode(InputSyntaxMode::SmilesLike);
    input_box_->when(FL_WHEN_CHANGED | FL_WHEN_ENTER_KEY);
    input_box_->callback(on_preview_cb, this);

    compare_input_label_ = new Fl_Box(panel_x + 10, panel_y + 212, panel_w - 20, 16, "Product SMILES");
    compare_input_label_->box(FL_NO_BOX);
    compare_input_label_->labelsize(12);
    compare_input_label_->labelcolor(ui(86, 97, 112));
    compare_input_label_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    compare_input_label_->hide();

    compare_input_box_ = new ColoredInputEditor(panel_x + 10, panel_y + 228, panel_w - 20, 48);
    compare_input_box_->value("");
    compare_input_box_->set_syntax_mode(InputSyntaxMode::SmilesLike);
    compare_input_box_->when(FL_WHEN_CHANGED | FL_WHEN_ENTER_KEY);
    compare_input_box_->callback(on_preview_cb, this);
    compare_input_box_->hide();

    auto *structure_title = new Fl_Box(panel_x + 10, panel_y + 216, panel_w - 214, 20, "Interactive Structure");
    structure_title->box(FL_NO_BOX);
    structure_title->labelsize(12);
    structure_title->labelcolor(ui(86, 97, 112));
    structure_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    edit_button_ = new Fl_Button(panel_x + panel_w - 196, panel_y + 214, 104, 24, "Edit xyzedit");
    edit_button_->callback(on_edit_structure_cb, this);
    edit_button_->box(FL_UP_BOX);
    edit_button_->color(ui(207, 218, 233));
    edit_button_->labelcolor(ui(63, 73, 86));
    edit_button_->labelfont(FL_HELVETICA_BOLD);
    edit_button_->labelsize(11);

    preview_button_ = new Fl_Button(panel_x + panel_w - 86, panel_y + 214, 76, 24, "Preview");
    preview_button_->callback(on_preview_cb, this);
    preview_button_->box(FL_UP_BOX);
    preview_button_->color(ui(198, 211, 228));
    preview_button_->labelcolor(ui(63, 73, 86));
    preview_button_->labelfont(FL_HELVETICA_BOLD);
    preview_button_->labelsize(11);

    structure_area_x_ = panel_x + 10;
    structure_area_y_ = panel_y + 240;
    structure_area_w_ = panel_w - 20;
    structure_area_h_ = 160;

    structure_current_label_ = new Fl_Box(structure_area_x_, structure_area_y_, structure_area_w_ / 2, 16, "");
    structure_current_label_->box(FL_NO_BOX);
    structure_current_label_->labelsize(11);
    structure_current_label_->labelcolor(ui(86, 112, 146));
    structure_current_label_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    structure_current_label_->hide();

    structure_compare_label_ =
        new Fl_Box(structure_area_x_ + structure_area_w_ / 2 + 4, structure_area_y_, structure_area_w_ / 2 - 4, 16, "");
    structure_compare_label_->box(FL_NO_BOX);
    structure_compare_label_->labelsize(11);
    structure_compare_label_->labelcolor(ui(164, 114, 86));
    structure_compare_label_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    structure_compare_label_->hide();

    structure_widget_ =
        new StructureWidget(structure_area_x_, structure_area_y_, structure_area_w_, structure_area_h_, nullptr);
    structure_widget_->set_highlight_palette(ui(201, 223, 245), ui(85, 121, 167), ui(222, 235, 251));
    structure_widget_->set_on_atom_selected([this](int atom_index, const std::vector<int> &hydrogens) {
        on_structure_atom_picked(atom_index, hydrogens);
    });

    compare_structure_widget_ =
        new StructureWidget(structure_area_x_, structure_area_y_, structure_area_w_, structure_area_h_, nullptr);
    compare_structure_widget_->set_highlight_palette(ui(244, 216, 198), ui(176, 114, 86), ui(252, 232, 220));
    compare_structure_widget_->set_on_atom_selected([this](int atom_index, const std::vector<int> &hydrogens) {
        on_compare_structure_atom_picked(atom_index, hydrogens);
    });
    compare_structure_widget_->set_empty_message("No comparison structure");
    compare_structure_widget_->hide();

    queue_button_ = new Fl_Button(panel_x + 10, panel_y + 406, 100, 30, "Queue");
    queue_button_->callback(on_queue_job_cb, this);
    queue_button_->box(FL_UP_BOX);
    queue_button_->color(ui(218, 225, 236));
    queue_button_->labelcolor(ui(63, 73, 86));
    queue_button_->labelfont(FL_HELVETICA_BOLD);

    start_button_ = new Fl_Button(panel_x + 116, panel_y + 406, 110, 30, "Run Pending");
    start_button_->callback(on_start_queue_cb, this);
    start_button_->box(FL_UP_BOX);
    start_button_->color(ui(168, 205, 194));
    start_button_->labelcolor(ui(52, 66, 66));
    start_button_->labelfont(FL_HELVETICA_BOLD);

    cancel_button_ = new Fl_Button(panel_x + 232, panel_y + 406, 92, 30, "Cancel");
    cancel_button_->callback(on_cancel_cb, this);
    cancel_button_->box(FL_UP_BOX);
    cancel_button_->color(ui(223, 229, 238));
    cancel_button_->labelcolor(ui(88, 98, 112));

    auto *queue_title = new Fl_Box(panel_x + 10, panel_y + 444, panel_w - 20, 20, "Queue (A=current, B=compare)");
    queue_title->box(FL_NO_BOX);
    queue_title->labelsize(12);
    queue_title->labelcolor(ui(86, 97, 112));
    queue_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    queue_browser_ = new Fl_Hold_Browser(panel_x + 10, panel_y + 468, panel_w - 20, 122);
    queue_browser_->callback(on_select_job_cb, this);
    queue_browser_->box(FL_DOWN_BOX);
    queue_browser_->color(ui(255, 255, 255));
    queue_browser_->textsize(12);
    queue_browser_->textfont(FL_HELVETICA_BOLD);
    queue_browser_->selection_color(ui(224, 236, 248));
    static int queue_col_widths[] = {16, 24, 40, 112, 36, 72, 0};
    queue_browser_->column_char('\t');
    queue_browser_->column_widths(queue_col_widths);

    auto *example_title = new Fl_Box(panel_x + 10, panel_y + 594, panel_w - 20, 18, "Examples");
    example_title->box(FL_NO_BOX);
    example_title->labelsize(12);
    example_title->labelcolor(ui(82, 96, 118));
    example_title->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    example_choice_ = new Fl_Choice(panel_x + 10, panel_y + 612, panel_w - 20, 24);
    example_choice_->box(FL_DOWN_BOX);
    example_choice_->color(ui(255, 255, 255));
    example_choice_->labelsize(11);
    example_choice_->add("Examples: none");
    example_choice_->value(0);
    example_choice_->callback(on_select_example_cb, this);

    load_example_calc_button_ = new Fl_Button(panel_x + 10, panel_y + 640, 102, 24, "Load Calc");
    load_example_calc_button_->callback(on_load_example_calc_cb, this);
    load_example_calc_button_->box(FL_UP_BOX);
    load_example_calc_button_->color(ui(206, 220, 238));
    load_example_calc_button_->labelcolor(ui(60, 74, 95));
    load_example_calc_button_->labelsize(11);

    load_example_bundle_button_ = new Fl_Button(panel_x + 116, panel_y + 640, 102, 24, "Calc+Exp");
    load_example_bundle_button_->callback(on_load_example_bundle_cb, this);
    load_example_bundle_button_->box(FL_UP_BOX);
    load_example_bundle_button_->color(ui(190, 220, 207));
    load_example_bundle_button_->labelcolor(ui(55, 82, 69));
    load_example_bundle_button_->labelsize(11);

    run_selected_button_ = new Fl_Button(panel_x + 224, panel_y + 640, 100, 24, "Run Selected");
    run_selected_button_->callback(on_run_selected_cb, this);
    run_selected_button_->box(FL_UP_BOX);
    run_selected_button_->color(ui(186, 204, 229));
    run_selected_button_->labelcolor(ui(59, 73, 94));
    run_selected_button_->labelfont(FL_HELVETICA_BOLD);
    run_selected_button_->labelsize(11);

    status_box_ = new Fl_Box(panel_x + 10, panel_y + 668, panel_w - 20, 26, "Idle");
    status_box_->box(FL_FLAT_BOX);
    status_box_->color(ui(229, 236, 246));
    status_box_->labelsize(11);
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
    spectrum_widget_->set_on_manual_shift_pair(
        [this](int primary_group_id, int compare_group_id) { on_manual_shift_pair(primary_group_id, compare_group_id); });

    spectrum_nucleus_choice_ = new Fl_Choice(right_x + 210, right_y + 6, 110, 24);
    spectrum_nucleus_choice_->box(FL_DOWN_BOX);
    spectrum_nucleus_choice_->color(ui(255, 255, 255));
    spectrum_nucleus_choice_->add("1H");
    spectrum_nucleus_choice_->value(0);
    spectrum_nucleus_choice_->callback(on_select_spectrum_nucleus_cb, this);

    load_experimental_button_ = new Fl_Button(right_x + 326, right_y + 6, 84, 24, "Load Exp");
    load_experimental_button_->callback(on_load_experimental_cb, this);
    load_experimental_button_->box(FL_UP_BOX);
    load_experimental_button_->color(ui(208, 221, 236));
    load_experimental_button_->labelcolor(ui(66, 78, 95));
    load_experimental_button_->labelsize(11);

    clear_experimental_button_ = new Fl_Button(right_x + 414, right_y + 6, 84, 24, "Clear Exp");
    clear_experimental_button_->callback(on_clear_experimental_cb, this);
    clear_experimental_button_->box(FL_UP_BOX);
    clear_experimental_button_->color(ui(226, 232, 241));
    clear_experimental_button_->labelcolor(ui(87, 98, 113));
    clear_experimental_button_->labelsize(11);

    const int top_gap = 6;
    const int toolbar_left = right_x + 504;
    const int toolbar_right = right_x + right_w - 10;
    const int toolbar_available = std::max(0, toolbar_right - toolbar_left);

    int export_w = 78;
    int experimental_w = 150;
    int reference_w = 150;
    const int min_export_w = 64;
    const int min_experimental_w = 110;
    const int min_reference_w = 104;

    int needed = export_w + experimental_w + reference_w + (2 * top_gap);
    if (needed > toolbar_available) {
        int overflow = needed - toolbar_available;
        const int ref_reduction = std::min(overflow, reference_w - min_reference_w);
        reference_w -= ref_reduction;
        overflow -= ref_reduction;

        const int exp_reduction = std::min(overflow, experimental_w - min_experimental_w);
        experimental_w -= exp_reduction;
        overflow -= exp_reduction;

        const int export_reduction = std::min(overflow, export_w - min_export_w);
        export_w -= export_reduction;
        overflow -= export_reduction;

        if (overflow > 0) {
            reference_w = std::max(min_reference_w, reference_w - overflow);
        }
    }

    const int reference_x = toolbar_right - reference_w;
    const int experimental_x = reference_x - top_gap - experimental_w;
    const int export_x = experimental_x - top_gap - export_w;

    export_spectrum_button_ = new Fl_Button(export_x, right_y + 6, export_w, 24, "Export");
    export_spectrum_button_->callback(on_export_spectrum_cb, this);
    export_spectrum_button_->box(FL_UP_BOX);
    export_spectrum_button_->color(ui(201, 216, 234));
    export_spectrum_button_->labelcolor(ui(60, 72, 89));
    export_spectrum_button_->labelsize(11);

    experimental_choice_ = new Fl_Choice(experimental_x, right_y + 6, experimental_w, 24);
    experimental_choice_->box(FL_DOWN_BOX);
    experimental_choice_->color(ui(255, 255, 255));
    experimental_choice_->labelsize(11);
    experimental_choice_->add("Exp: none");
    experimental_choice_->value(0);
    experimental_choice_->callback(on_select_experimental_cb, this);

    reference_choice_ = new Fl_Choice(reference_x, right_y + 6, reference_w, 24);
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

    workflow_info_bg_ = new Fl_Box(FL_DOWN_BOX, right_x, right_y + spectrum_h + 180, right_w, 52, "");
    workflow_info_bg_->color(ui(255, 255, 255));

    workflow_info_line1_ = new Fl_Box(right_x + 8, right_y + spectrum_h + 182, right_w - 16, 16, "");
    workflow_info_line1_->box(FL_NO_BOX);
    workflow_info_line1_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    workflow_info_line1_->labelsize(12);
    workflow_info_line1_->labelfont(FL_HELVETICA_BOLD);
    workflow_info_line1_->labelcolor(ui(84, 95, 112));

    workflow_info_line2_ = new Fl_Box(right_x + 8, right_y + spectrum_h + 198, right_w - 16, 16, "");
    workflow_info_line2_->box(FL_NO_BOX);
    workflow_info_line2_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    workflow_info_line2_->labelsize(12);
    workflow_info_line2_->labelcolor(ui(98, 110, 128));

    workflow_info_line3_ = new Fl_Box(right_x + 8, right_y + spectrum_h + 214, right_w - 16, 16, "");
    workflow_info_line3_->box(FL_NO_BOX);
    workflow_info_line3_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    workflow_info_line3_->labelsize(11);
    workflow_info_line3_->labelcolor(ui(120, 132, 149));

    end();
    refresh_example_choice();
    refresh_experimental_choice();
    apply_active_experimental_overlay();
    update_structure_compare_layout(false);
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

int AppWindow::handle(int event) {
    if (event == FL_PUSH && Fl::event_button() == FL_RIGHT_MOUSE) {
        if (maybe_show_queue_context_menu()) {
            return 1;
        }
    }

    if (event == FL_SHORTCUT) {
        const int state = Fl::event_state();
        const bool command_or_ctrl = ((state & FL_COMMAND) != 0) || ((state & FL_CTRL) != 0);
        if (command_or_ctrl) {
            const int key = Fl::event_key();
            if (key == 'e' || key == 'E') {
                on_export_spectrum();
                return 1;
            }
            if (key == 'l' || key == 'L') {
                on_load_experimental();
                return 1;
            }
            if (key == 'r' || key == 'R') {
                start_queue();
                return 1;
            }
            if (key == 'p' || key == 'P') {
                preview_current_input(true);
                return 1;
            }
            if (key == FL_Enter || key == FL_KP_Enter) {
                queue_current_input();
                return 1;
            }
        }
    }
    return Fl_Double_Window::handle(event);
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
    self->refresh_input_syntax_mode();
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

void AppWindow::on_select_experimental_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_experimental();
}

void AppWindow::on_load_experimental_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_load_experimental();
}

void AppWindow::on_clear_experimental_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_clear_experimental();
}

void AppWindow::on_export_spectrum_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_export_spectrum();
}

void AppWindow::on_select_example_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_select_example();
}

void AppWindow::on_load_example_calc_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_load_example(false);
}

void AppWindow::on_load_example_bundle_cb(Fl_Widget *, void *userdata) {
    static_cast<AppWindow *>(userdata)->on_load_example(true);
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

void AppWindow::refresh_input_syntax_mode() {
    if (input_box_ == nullptr) {
        return;
    }

    const std::string raw_input = input_box_->value() != nullptr ? input_box_->value() : "";
    const InputFormat detected = detect_input_format_auto(raw_input);
    if (detected == InputFormat::Xyz) {
        input_box_->set_syntax_mode(InputSyntaxMode::XyzLike);
    } else {
        input_box_->set_syntax_mode(InputSyntaxMode::SmilesLike);
    }
}

void AppWindow::queue_current_input() {
    refresh_input_syntax_mode();
    std::string raw_input;
    if (input_box_ != nullptr) {
        const char *raw_ptr = input_box_->value();
        if (raw_ptr != nullptr) {
            raw_input = raw_ptr;
        }
    }
    const InputFormat detected_format = detect_input_format_auto(raw_input);
    if (detected_format == InputFormat::Unknown) {
        status_box_->label("Cannot queue: input is empty");
        return;
    }
    std::vector<std::string> inputs;
    if (detected_format == InputFormat::Smiles) {
        inputs = split_smiles_inputs(raw_input);
    } else {
        const std::string clean = trim_copy(raw_input);
        if (!clean.empty()) {
            inputs.push_back(raw_input);
        }
    }
    if (inputs.empty()) {
        status_box_->label("Cannot queue: input is empty");
        return;
    }

    const std::string base_name = [&]() {
        const std::string clean_name = trim_copy(job_name_input_ != nullptr ? job_name_input_->value() : "");
        if (!clean_name.empty()) {
            return clean_name;
        }
        return std::string("untitled");
    }();

    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    const char *workflow = workflow_choice_ != nullptr ? workflow_choice_->text(workflow_choice_->value()) : "all";
    WorkflowKind workflow_kind = workflow_kind_from_string(workflow ? workflow : "all");
    if (workflow_kind == WorkflowKind::Unknown) {
        workflow_kind = WorkflowKind::All;
    }

    std::size_t queued_index = 0;
    std::size_t queued_count = 0;
    bool queued_any = false;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (jobs_.size() + inputs.size() > 50) {
            status_box_->label("Queue limit reached (50)");
            return;
        }

        for (std::size_t i = 0; i < inputs.size(); ++i) {
            QueuedJob job;
            job.status = "pending";
            job.config.job_name = (inputs.size() > 1) ? (base_name + "_" + std::to_string(i + 1)) : base_name;
            job.config.input_value = inputs[i];
            job.config.workflow_kind = workflow_kind;
            job.config.input_format = detected_format;
            job.config.solvent = solvent ? solvent : "cdcl3";
            job.config.nucleus = "auto";
            job.config.line_shape = "lorentzian";
            if (workflow_kind == WorkflowKind::Compare && compare_input_box_ != nullptr) {
                const char *cmp_ptr = compare_input_box_->value();
                job.config.compare_input_value = cmp_ptr != nullptr ? cmp_ptr : "";
                job.config.compare_input_format = InputFormat::Smiles;
            }
            job.id = "q" + std::to_string(jobs_.size() + 1);
            jobs_.push_back(std::move(job));
            if (!queued_any) {
                queued_index = jobs_.size() - 1;
                queued_any = true;
            }
            ++queued_count;
        }
        if (worker_running_ && run_scope_ == RunScope::PendingQueue && queued_count > 0) {
            run_total_ += static_cast<int>(queued_count);
        }
    }

    if (queued_count == 1) {
        status_box_->label("Job queued; generating structure preview");
    } else {
        status_box_->copy_label(("Queued " + std::to_string(queued_count) + " jobs from SMILES list").c_str());
    }
    refresh_queue_browser();
    if (queued_any) {
        queue_browser_->value(static_cast<int>(queued_index + 1));
        request_preview_for_job_index(queued_index, true);
    }
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
        active_nucleus_ = (cfg.workflow_kind == WorkflowKind::Cd) ? "CD" : normalize_nucleus_label(cfg.nucleus);
        apply_reference_peaks(cfg.solvent, active_nucleus_);
        spectrum_widget_->set_nucleus_label(spectrum_title_for_label(active_nucleus_));
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
    refresh_input_syntax_mode();

    if (worker_running_) {
        if (show_status) {
            status_box_->label("Preview is paused while calculations are running");
        }
        return;
    }

    JobConfig cfg;
    cfg.job_name = job_name_input_->value();
    std::string raw_input;
    if (input_box_ != nullptr) {
        const char *raw_ptr = input_box_->value();
        if (raw_ptr != nullptr) {
            raw_input = raw_ptr;
        }
    }
    const InputFormat detected_format = detect_input_format_auto(raw_input);
    if (detected_format == InputFormat::Unknown) {
        return;
    }
    if (detected_format == InputFormat::Smiles) {
        const auto smiles_entries = split_smiles_inputs(raw_input);
        if (smiles_entries.size() > 1) {
            if (structure_widget_ != nullptr) {
                structure_widget_->set_empty_message("Multiple SMILES detected");
                structure_widget_->set_structure({}, {});
            }
            if (show_status) {
                status_box_->label("Multiple SMILES detected: queue to run all entries");
            }
            return;
        }
        if (!smiles_entries.empty()) {
            cfg.input_value = smiles_entries.front();
        }
    } else {
        cfg.input_value = raw_input;
    }
    if (structure_widget_ != nullptr) {
        structure_widget_->set_empty_message("No 2D structure");
    }

    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    const char *workflow = workflow_choice_ != nullptr ? workflow_choice_->text(workflow_choice_->value()) : "all";
    cfg.workflow_kind = workflow_kind_from_string(workflow ? workflow : "all");
    if (cfg.workflow_kind == WorkflowKind::Unknown) {
        cfg.workflow_kind = WorkflowKind::All;
    }
    cfg.input_format = detected_format;
    cfg.solvent = solvent ? solvent : "cdcl3";
    cfg.nucleus = "auto";
    cfg.need_editable_xyz = false;

    if (trim_copy(cfg.input_value).empty() || cfg.input_format == InputFormat::Unknown) {
        return;
    }
    launch_preview_async(cfg, show_status);
}

void AppWindow::edit_current_structure() {
    refresh_input_syntax_mode();

    if (worker_running_) {
        status_box_->label("Editor is disabled while calculations are running");
        return;
    }

    JobConfig cfg;
    cfg.job_name = job_name_input_->value();
    std::string raw_input;
    if (input_box_ != nullptr) {
        const char *raw_ptr = input_box_->value();
        if (raw_ptr != nullptr) {
            raw_input = raw_ptr;
        }
    }
    const InputFormat detected_format = detect_input_format_auto(raw_input);
    if (detected_format == InputFormat::Smiles) {
        const auto smiles_entries = split_smiles_inputs(raw_input);
        if (smiles_entries.size() > 1) {
            if (structure_widget_ != nullptr) {
                structure_widget_->set_empty_message("Multiple SMILES detected");
                structure_widget_->set_structure({}, {});
            }
            status_box_->label("Edit xyzedit supports one molecule at a time");
            return;
        }
        if (!smiles_entries.empty()) {
            cfg.input_value = smiles_entries.front();
        }
    } else {
        cfg.input_value = raw_input;
    }
    if (structure_widget_ != nullptr) {
        structure_widget_->set_empty_message("No 2D structure");
    }

    const char *solvent = solvent_choice_->text(solvent_choice_->value());
    const char *workflow = workflow_choice_ != nullptr ? workflow_choice_->text(workflow_choice_->value()) : "all";
    cfg.workflow_kind = workflow_kind_from_string(workflow ? workflow : "all");
    if (cfg.workflow_kind == WorkflowKind::Unknown) {
        cfg.workflow_kind = WorkflowKind::All;
    }
    cfg.input_format = detected_format;
    cfg.solvent = solvent ? solvent : "cdcl3";
    cfg.nucleus = "auto";
    cfg.need_editable_xyz = true;

    if (trim_copy(cfg.input_value).empty() || cfg.input_format == InputFormat::Unknown) {
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
    active_nucleus_ = (cfg.workflow_kind == WorkflowKind::Cd) ? "CD" : normalize_nucleus_label(cfg.nucleus);
    apply_reference_peaks(cfg.solvent, active_nucleus_);
    spectrum_widget_->set_nucleus_label(spectrum_title_for_label(active_nucleus_));

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
    refresh_input_syntax_mode();

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
    selected_job_index_ = (selected > 0) ? (selected - 1) : -1;
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
    maybe_apply_example_overlay_for_active_selection();
    refresh_queue_browser();
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
        highlight_comparison_for_primary_group(group_id);
    } else if (compare_structure_widget_ != nullptr) {
        compare_structure_widget_->clear_highlight();
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
    maybe_apply_example_overlay_for_active_selection();
}

void AppWindow::on_select_experimental() {
    if (experimental_choice_ == nullptr) {
        return;
    }
    const int selected = experimental_choice_->value();
    if (selected < 0 || selected >= static_cast<int>(experimental_choice_keys_.size())) {
        return;
    }
    active_experimental_overlay_key_ = experimental_choice_keys_[static_cast<std::size_t>(selected)];
    apply_active_experimental_overlay();

    if (active_experimental_overlay_key_.empty()) {
        status_box_->label("Experimental overlay: none");
        return;
    }

    const auto points_it = experimental_overlays_.find(active_experimental_overlay_key_);
    const auto format_it = experimental_overlay_formats_.find(active_experimental_overlay_key_);
    if (points_it == experimental_overlays_.end()) {
        status_box_->label("Experimental overlay selection unavailable");
        return;
    }

    std::ostringstream status;
    status << "Experimental overlay: " << active_experimental_overlay_key_
           << " (" << points_it->second.size() << " points";
    if (format_it != experimental_overlay_formats_.end() && !format_it->second.empty()) {
        status << ", " << format_it->second;
    }
    status << ")";
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_select_example() {
    if (example_choice_ == nullptr) {
        return;
    }
    const int selected = example_choice_->value();
    if (selected <= 0 || selected >= static_cast<int>(example_choice_case_indices_.size())) {
        return;
    }
    const int bundle_idx = example_choice_case_indices_[static_cast<std::size_t>(selected)];
    if (bundle_idx < 0 || bundle_idx >= static_cast<int>(example_bundle_row_indices_.size())) {
        return;
    }

    std::set<std::string> products;
    std::string workflow;
    for (int row_idx : example_bundle_row_indices_[static_cast<std::size_t>(bundle_idx)]) {
        if (row_idx < 0 || row_idx >= static_cast<int>(example_case_rows_.size())) {
            continue;
        }
        const auto fields = parse_csv_line(example_case_rows_[static_cast<std::size_t>(row_idx)]);
        if (fields.size() < 10) {
            continue;
        }
        products.insert(normalize_nucleus_label(fields[3]));
        if (workflow.empty()) {
            workflow = trim_copy(fields[2]);
        }
    }

    std::ostringstream products_text;
    bool first = true;
    for (const auto &p : products) {
        if (!first) {
            products_text << ", ";
        }
        products_text << p;
        first = false;
    }

    const std::string bundle_name = example_bundle_names_[static_cast<std::size_t>(bundle_idx)];
    std::ostringstream status;
    status << "Example selected: " << bundle_name;
    if (!workflow.empty()) {
        status << " | " << uppercase_copy(workflow);
    }
    if (!products.empty()) {
        status << " | products: " << products_text.str();
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_load_example(bool with_experimental) {
    if (example_choice_ == nullptr) {
        status_box_->label("No example menu available");
        return;
    }
    const int selected = example_choice_->value();
    if (selected <= 0 || selected >= static_cast<int>(example_choice_case_indices_.size())) {
        status_box_->label("Select an example first");
        return;
    }
    const int bundle_idx = example_choice_case_indices_[static_cast<std::size_t>(selected)];
    if (bundle_idx < 0 || bundle_idx >= static_cast<int>(example_bundle_row_indices_.size())) {
        status_box_->label("Selected example is unavailable");
        return;
    }
    const auto &row_indices = example_bundle_row_indices_[static_cast<std::size_t>(bundle_idx)];
    if (row_indices.empty()) {
        status_box_->label("Selected example has no spectra rows");
        return;
    }

    std::map<std::string, std::string> computed_by_product;
    std::map<std::string, std::string> experimental_by_product;
    std::set<std::string> products;
    std::string smiles;
    bool has_nmr = false;
    bool has_cd = false;

    namespace fs = std::filesystem;
    std::error_code ec;
    for (int row_idx : row_indices) {
        if (row_idx < 0 || row_idx >= static_cast<int>(example_case_rows_.size())) {
            continue;
        }
        const auto fields = parse_csv_line(example_case_rows_[static_cast<std::size_t>(row_idx)]);
        if (fields.size() < 10) {
            continue;
        }
        if (smiles.empty()) {
            smiles = trim_copy(fields[5]);
        }

        const std::string workflow = lowercase_copy(fields[2]);
        if (workflow == "cd") {
            has_cd = true;
        } else {
            has_nmr = true;
        }

        const std::string product = normalize_nucleus_label(fields[3]);
        if (product.empty()) {
            continue;
        }
        const std::string computed_ref = resolve_example_asset_path(fields[6]);
        ec.clear();
        if (computed_ref.empty() || !fs::exists(computed_ref, ec) || ec) {
            status_box_->copy_label(("Example missing computed spectrum file: " + trim_copy(fields[6])).c_str());
            return;
        }
        if (computed_by_product.find(product) == computed_by_product.end()) {
            computed_by_product[product] = computed_ref;
            products.insert(product);
        }

        const std::string exp_file = resolve_example_asset_path(fields[7]);
        if (!exp_file.empty() && experimental_by_product.find(product) == experimental_by_product.end()) {
            experimental_by_product[product] = exp_file;
        }
    }

    if (computed_by_product.empty()) {
        status_box_->label("Selected example has no valid computed spectra");
        return;
    }

    const std::string default_product = pick_default_product(products);
    const auto default_it = computed_by_product.find(default_product);
    if (default_it == computed_by_product.end()) {
        status_box_->label("Selected example has no default product");
        return;
    }

    const std::string manifest_path = write_example_manifest(
        example_bundle_names_[static_cast<std::size_t>(bundle_idx)],
        computed_by_product);
    if (manifest_path.empty()) {
        status_box_->label("Failed to build example spectra manifest");
        return;
    }

    QueuedJob job;
    job.status = "done";
    job.config.job_name = example_bundle_names_[static_cast<std::size_t>(bundle_idx)].empty()
                              ? "example"
                              : example_bundle_names_[static_cast<std::size_t>(bundle_idx)];
    job.config.input_value = smiles;
    job.config.input_format = InputFormat::Smiles;
    if (has_nmr && has_cd) {
        job.config.workflow_kind = WorkflowKind::All;
    } else if (has_cd) {
        job.config.workflow_kind = WorkflowKind::Cd;
    } else {
        job.config.workflow_kind = WorkflowKind::Nmr;
    }
    job.config.solvent = solvent_choice_ != nullptr ? solvent_choice_->text(solvent_choice_->value()) : "cdcl3";
    job.config.nucleus = default_product;
    job.id.clear();
    job.message = "example bundle loaded";
    job.spectrum_csv = default_it->second;
    job.spectra_manifest_csv = manifest_path;

    int row_number = 0;
    int job_index = -1;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        job.id = "q" + std::to_string(jobs_.size() + 1);
        jobs_.push_back(job);
        row_number = static_cast<int>(jobs_.size());
        job_index = row_number - 1;
    }

    if (job_index >= 0) {
        std::map<std::string, std::string> overlay_key_by_product;
        std::string overlay_warning;
        if (with_experimental) {
            for (const auto &entry : experimental_by_product) {
                const std::string &product = entry.first;
                const std::string &exp_file = entry.second;
                if (exp_file.empty()) {
                    continue;
                }
                ec.clear();
                if (!fs::exists(exp_file, ec) || ec) {
                    if (overlay_warning.empty()) {
                        overlay_warning = "Missing overlay: " + exp_file;
                    }
                    continue;
                }
                const auto loaded = load_experimental_spectrum(exp_file);
                if (!loaded.error_message.empty()) {
                    if (overlay_warning.empty()) {
                        overlay_warning = experimental_import_failure_status(loaded);
                    }
                    continue;
                }
                std::string base_label = example_bundle_names_[static_cast<std::size_t>(bundle_idx)];
                if (base_label.empty()) {
                    base_label = "example";
                }
                base_label += " " + product;
                const std::string overlay_key = make_unique_overlay_key(experimental_overlays_, base_label);
                experimental_overlays_[overlay_key] = loaded.points;
                experimental_overlay_paths_[overlay_key] = exp_file;
                experimental_overlay_formats_[overlay_key] = loaded.detected_format;
                overlay_key_by_product[product] = overlay_key;
            }
        }
        if (!overlay_key_by_product.empty()) {
            example_job_overlay_keys_[job_index] = overlay_key_by_product;
        } else {
            example_job_overlay_keys_.erase(job_index);
        }
        if (!overlay_warning.empty()) {
            status_box_->copy_label(overlay_warning.c_str());
        }
    }

    refresh_queue_browser();
    queue_browser_->value(row_number);
    on_select_job();
    maybe_apply_example_overlay_for_active_selection();

    std::ostringstream status;
    status << "Loaded example bundle: " << job.config.job_name << " | products ";
    bool first = true;
    for (const auto &product : products) {
        if (!first) {
            status << ",";
        }
        status << product;
        first = false;
    }
    if (with_experimental) {
        const auto overlay_it = example_job_overlay_keys_.find(job_index);
        const std::size_t count = (overlay_it != example_job_overlay_keys_.end()) ? overlay_it->second.size() : 0;
        status << " | overlays " << count;
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_load_experimental() {
    Fl_Native_File_Chooser chooser;
    chooser.title("Load Experimental Spectrum");
    chooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
    chooser.filter("Text/CSV\t*.csv\nText\t*.txt\nASCII\t*.asc\nData\t*.dat\nAll\t*");
    const int rc = chooser.show();
    if (rc != 0) {
        return;
    }

    const char *filename = chooser.filename();
    if (filename == nullptr || std::string(filename).empty()) {
        return;
    }

    const auto loaded = load_experimental_spectrum(filename);
    if (!loaded.error_message.empty()) {
        status_box_->copy_label(experimental_import_failure_status(loaded).c_str());
        return;
    }

    const std::string base_label = make_overlay_label_from_path(filename);
    const std::string overlay_key = make_unique_overlay_key(experimental_overlays_, base_label);
    experimental_overlays_[overlay_key] = loaded.points;
    experimental_overlay_paths_[overlay_key] = filename;
    experimental_overlay_formats_[overlay_key] = loaded.detected_format;
    active_experimental_overlay_key_ = overlay_key;
    refresh_experimental_choice();
    apply_active_experimental_overlay();

    std::ostringstream status;
    status << "Loaded experimental spectrum: " << overlay_key << " (" << loaded.points.size() << " points";
    if (!loaded.detected_format.empty()) {
        status << ", " << loaded.detected_format;
    }
    status << ")";
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_clear_experimental() {
    if (experimental_overlays_.empty()) {
        active_experimental_overlay_key_.clear();
        refresh_experimental_choice();
        apply_active_experimental_overlay();
        status_box_->label("No experimental overlays loaded");
        return;
    }

    if (active_experimental_overlay_key_.empty()) {
        experimental_overlays_.clear();
        experimental_overlay_paths_.clear();
        experimental_overlay_formats_.clear();
        refresh_experimental_choice();
        apply_active_experimental_overlay();
        status_box_->label("Cleared all experimental overlays");
        return;
    }

    const std::string removed = active_experimental_overlay_key_;
    experimental_overlays_.erase(removed);
    experimental_overlay_paths_.erase(removed);
    experimental_overlay_formats_.erase(removed);
    if (!experimental_overlays_.empty()) {
        active_experimental_overlay_key_ = experimental_overlays_.begin()->first;
    } else {
        active_experimental_overlay_key_.clear();
    }
    refresh_experimental_choice();
    apply_active_experimental_overlay();
    status_box_->copy_label(("Removed experimental overlay: " + removed).c_str());
}

void AppWindow::on_export_spectrum() {
    if (spectrum_widget_ == nullptr) {
        status_box_->label("No spectrum widget available to export");
        return;
    }

    Fl_Native_File_Chooser chooser;
    chooser.title("Export Spectrum Plot");
    chooser.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
    chooser.filter("PNG image\t*.png\nPPM image\t*.ppm\nAll\t*");
    chooser.options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    const int rc = chooser.show();
    if (rc != 0) {
        return;
    }

    const char *filename = chooser.filename();
    if (filename == nullptr || std::string(filename).empty()) {
        return;
    }

    std::string path = filename;
    namespace fs = std::filesystem;
    fs::path out_path(path);
    std::string ext = out_path.has_extension() ? out_path.extension().string() : std::string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ext.empty()) {
        out_path += ".png";
        ext = ".png";
    }

    std::string error_message;
    if (ext == ".png") {
        if (!export_widget_snapshot_png(spectrum_widget_, out_path.string(), &error_message)) {
            if (error_message.empty()) {
                error_message = "Unknown export failure.";
            }
            status_box_->copy_label(("Export failed: " + error_message).c_str());
            return;
        }
    } else {
        if (!export_widget_snapshot_ppm(spectrum_widget_, out_path.string(), &error_message)) {
            if (error_message.empty()) {
                error_message = "Unknown export failure.";
            }
            status_box_->copy_label(("Export failed: " + error_message).c_str());
            return;
        }
    }

    status_box_->copy_label(("Exported spectrum plot: " + out_path.string()).c_str());
}

bool AppWindow::maybe_show_queue_context_menu() {
    if (queue_browser_ == nullptr) {
        return false;
    }
    const int mx = Fl::event_x();
    const int my = Fl::event_y();
    if (mx < queue_browser_->x() || mx > queue_browser_->x() + queue_browser_->w()
        || my < queue_browser_->y() || my > queue_browser_->y() + queue_browser_->h()) {
        return false;
    }

    const int current_line = queue_browser_->value();
    int target_line = browser_line_at_mouse(queue_browser_, my);
    if (target_line <= 0) {
        target_line = current_line;
    }
    if (target_line <= 0) {
        return true;
    }

    const int current_index = (current_line > 0) ? (current_line - 1) : -1;
    const int target_index = target_line - 1;

    QueuedJob current_job;
    QueuedJob target_job;
    bool has_current = false;
    bool has_target = false;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (current_index >= 0 && current_index < static_cast<int>(jobs_.size())) {
            current_job = jobs_[static_cast<std::size_t>(current_index)];
            has_current = true;
        }
        if (target_index >= 0 && target_index < static_cast<int>(jobs_.size())) {
            target_job = jobs_[static_cast<std::size_t>(target_index)];
            has_target = true;
        }
    }

    const std::string current_name = has_current ? truncate_text(current_job.config.job_name, 20) : "selected job";
    const std::string compare_label = "Compare to " + current_name;
    const std::string clear_label = "Clear comparison";
    const bool compare_enabled =
        has_current && has_target && current_index != target_index
        && current_job.status == "done" && target_job.status == "done";
    const bool clear_enabled = comparison_job_index_ >= 0;

    Fl_Menu_Item menu[] = {
        {compare_label.c_str(), 0, nullptr, nullptr, compare_enabled ? 0 : FL_MENU_INACTIVE},
        {clear_label.c_str(), 0, nullptr, nullptr, clear_enabled ? 0 : FL_MENU_INACTIVE},
        {nullptr}
    };

    const Fl_Menu_Item *picked = menu->popup(mx, my);
    if (picked == nullptr) {
        return true;
    }

    if (picked == &menu[0] && compare_enabled) {
        comparison_job_index_ = target_index;
        manual_shift_pairs_.clear();
        if (spectrum_widget_ != nullptr) {
            spectrum_widget_->clear_manual_shift_pairs();
        }
        load_selected_job_visuals();
        refresh_queue_browser();
        status_box_->copy_label(
            ("Comparing " + target_job.config.job_name + " to " + current_job.config.job_name).c_str());
        return true;
    }

    if (picked == &menu[1] && clear_enabled) {
        clear_comparison_state(true);
        refresh_queue_browser();
        status_box_->label("Cleared calculated comparison");
        return true;
    }

    return true;
}

void AppWindow::clear_comparison_state(bool clear_selected_job) {
    if (clear_selected_job) {
        comparison_job_index_ = -1;
    }
    comparison_peak_rows_.clear();
    manual_shift_pairs_.clear();
    comparison_group_to_ppm_.clear();
    comparison_group_to_atoms_.clear();
    if (spectrum_widget_ != nullptr) {
        spectrum_widget_->clear_comparison_points();
    }
    if (compare_structure_widget_ != nullptr) {
        compare_structure_widget_->set_structure({}, {});
        compare_structure_widget_->clear_highlight();
    }
    update_structure_compare_layout(false);
}

void AppWindow::update_structure_compare_layout(bool compare_mode) {
    if (structure_widget_ == nullptr || compare_structure_widget_ == nullptr
        || structure_current_label_ == nullptr || structure_compare_label_ == nullptr) {
        return;
    }

    if (!compare_mode) {
        structure_current_label_->hide();
        structure_compare_label_->hide();
        compare_structure_widget_->hide();
        structure_widget_->resize(structure_area_x_, structure_area_y_, structure_area_w_, structure_area_h_);
        return;
    }

    const int label_h = 16;
    const int gap = 8;
    const int split_w = (structure_area_w_ - gap) / 2;
    const int widget_y = structure_area_y_ + label_h + 2;
    const int widget_h = std::max(56, structure_area_h_ - label_h - 2);

    structure_current_label_->resize(structure_area_x_, structure_area_y_, split_w, label_h);
    structure_compare_label_->resize(structure_area_x_ + split_w + gap, structure_area_y_, split_w, label_h);
    structure_current_label_->show();
    structure_compare_label_->show();

    structure_widget_->resize(structure_area_x_, widget_y, split_w, widget_h);
    compare_structure_widget_->resize(structure_area_x_ + split_w + gap, widget_y, split_w, widget_h);
    compare_structure_widget_->show();
}

void AppWindow::highlight_comparison_for_primary_group(int primary_group_id) {
    if (compare_structure_widget_ == nullptr || comparison_group_to_atoms_.empty()) {
        return;
    }

    int compare_group_id = -1;
    for (const auto &pair : manual_shift_pairs_) {
        if (pair.primary_group_id == primary_group_id) {
            compare_group_id = pair.compare_group_id;
            break;
        }
    }
    if (compare_group_id <= 0) {
        compare_structure_widget_->clear_highlight();
        return;
    }

    auto atoms_it = comparison_group_to_atoms_.find(compare_group_id);
    if (atoms_it == comparison_group_to_atoms_.end() || atoms_it->second.empty()) {
        compare_structure_widget_->clear_highlight();
        return;
    }

    if (is_proton_nucleus(active_nucleus_)) {
        compare_structure_widget_->set_selected_atom(-1);
        compare_structure_widget_->set_highlight_hydrogens(atoms_it->second);
    } else {
        compare_structure_widget_->set_highlight_hydrogens({});
        compare_structure_widget_->set_selected_atom(atoms_it->second.front());
    }
}

void AppWindow::apply_comparison_visuals(const QueuedJob &active_job) {
    if (spectrum_widget_ == nullptr) {
        return;
    }

    if (comparison_job_index_ < 0 || selected_job_index_ < 0 || comparison_job_index_ == selected_job_index_) {
        clear_comparison_state(false);
        return;
    }

    QueuedJob compare_job;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        if (comparison_job_index_ < 0 || comparison_job_index_ >= static_cast<int>(jobs_.size())) {
            clear_comparison_state(true);
            return;
        }
        compare_job = jobs_[static_cast<std::size_t>(comparison_job_index_)];
    }

    if (compare_job.status != "done") {
        clear_comparison_state(false);
        status_box_->copy_label(("Comparison job not ready: " + compare_job.config.job_name).c_str());
        return;
    }

    const auto compare_products = spectral_products_for_job(compare_job);
    auto it = compare_products.find(active_nucleus_);
    if (it == compare_products.end() || it->second.spectrum_csv.empty()) {
        clear_comparison_state(false);
        status_box_->copy_label(
            ("Comparison job has no " + active_nucleus_ + " spectrum: " + compare_job.config.job_name).c_str());
        return;
    }

    const auto compare_points = load_spectrum_csv(it->second.spectrum_csv);
    if (compare_points.empty()) {
        clear_comparison_state(false);
        status_box_->copy_label(
            ("Comparison spectrum unreadable for " + active_nucleus_ + ": " + compare_job.config.job_name).c_str());
        return;
    }

    comparison_peak_rows_.clear();
    comparison_group_to_ppm_.clear();
    comparison_group_to_atoms_.clear();
    const auto compare_peak_rows = read_peak_rows(it->second.peaks_csv);
    for (const auto &row : compare_peak_rows) {
        ComparisonPeakMarker marker;
        marker.group_id = row.group_id;
        marker.center_ppm = row.center_ppm;
        marker.label = row.multiplicity;
        comparison_peak_rows_.push_back(marker);
        comparison_group_to_ppm_[row.group_id] = row.center_ppm;
    }
    comparison_group_to_atoms_ = read_assignments(it->second.assignments_csv);

    manual_shift_pairs_.clear();
    spectrum_widget_->set_comparison_points(compare_points);
    spectrum_widget_->set_comparison_peak_markers(comparison_peak_rows_);
    spectrum_widget_->set_manual_shift_pairs(manual_shift_pairs_);

    if (compare_structure_widget_ != nullptr) {
        const auto compare_atoms = read_structure_atoms(compare_job.structure_atoms_csv);
        const auto compare_bonds = read_structure_bonds(compare_job.structure_bonds_csv);
        compare_structure_widget_->set_structure(compare_atoms, compare_bonds);
        compare_structure_widget_->set_empty_message("No comparison structure");
    }
    if (structure_current_label_ != nullptr) {
        structure_current_label_->copy_label(("A: " + truncate_text(active_job.config.job_name, 14)).c_str());
    }
    if (structure_compare_label_ != nullptr) {
        structure_compare_label_->copy_label(("B: " + truncate_text(compare_job.config.job_name, 14)).c_str());
    }
    update_structure_compare_layout(true);

    std::ostringstream status;
    status << "Comparing " << compare_job.config.job_name << " vs " << active_job.config.job_name
           << " (" << active_nucleus_ << ")";
    if (!comparison_peak_rows_.empty()) {
        status << " | " << comparison_peak_rows_.size() << " compare peaks";
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_manual_shift_pair(int primary_group_id, int compare_group_id) {
    auto primary_it = primary_group_to_ppm_.find(primary_group_id);
    auto compare_it = comparison_group_to_ppm_.find(compare_group_id);
    if (primary_it == primary_group_to_ppm_.end() || compare_it == comparison_group_to_ppm_.end()) {
        status_box_->label("Shift pair unavailable for selected peaks");
        return;
    }

    const double primary_ppm = primary_it->second;
    const double compare_ppm = compare_it->second;
    const double delta_ppm = compare_ppm - primary_ppm;

    manual_shift_pairs_.erase(
        std::remove_if(
            manual_shift_pairs_.begin(),
            manual_shift_pairs_.end(),
            [primary_group_id, compare_group_id](const ManualShiftPair &pair) {
                return pair.primary_group_id == primary_group_id || pair.compare_group_id == compare_group_id;
            }),
        manual_shift_pairs_.end());

    ManualShiftPair pair;
    pair.primary_group_id = primary_group_id;
    pair.compare_group_id = compare_group_id;
    pair.primary_ppm = primary_ppm;
    pair.compare_ppm = compare_ppm;
    pair.delta_ppm = delta_ppm;
    manual_shift_pairs_.push_back(pair);
    spectrum_widget_->set_manual_shift_pairs(manual_shift_pairs_);
    on_peak_picked(primary_group_id);

    std::ostringstream status;
    status << "Paired G" << primary_group_id << " with cmp G" << compare_group_id
           << " | delta " << std::showpos << std::fixed << std::setprecision(2) << delta_ppm
           << " ppm " << (delta_ppm > 0.0 ? "(downfield)" : "(upfield)");
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::on_peak_picked(int group_id) {
    if (group_id <= 0) {
        spectrum_widget_->set_selected_groups({});
        peak_browser_->value(0);
        atom_browser_->value(0);
        highlight_hydrogen_rows({});
        if (structure_widget_ != nullptr) {
            structure_widget_->set_selected_atom(-1);
            structure_widget_->set_highlight_hydrogens({});
        }
        if (compare_structure_widget_ != nullptr) {
            compare_structure_widget_->clear_highlight();
        }
        status_box_->label("Cleared peak highlight");
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
    highlight_comparison_for_primary_group(group_id);

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
        highlight_comparison_for_primary_group(primary_group);
    } else if (compare_structure_widget_ != nullptr) {
        compare_structure_widget_->clear_highlight();
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

void AppWindow::on_compare_structure_atom_picked(int atom_index, const std::vector<int> &attached_hydrogens) {
    if (atom_index <= 0 || compare_structure_widget_ == nullptr) {
        return;
    }

    std::vector<int> nucleus_targets = attached_hydrogens;
    if (!is_proton_nucleus(active_nucleus_)) {
        nucleus_targets = {atom_index};
    }

    compare_structure_widget_->set_selected_atom(atom_index);
    if (is_proton_nucleus(active_nucleus_)) {
        compare_structure_widget_->set_highlight_hydrogens(nucleus_targets);
    } else {
        compare_structure_widget_->set_highlight_hydrogens({});
    }

    std::vector<int> matched_compare_groups;
    for (const auto &entry : comparison_group_to_atoms_) {
        bool intersects = false;
        for (int target_atom : nucleus_targets) {
            if (std::find(entry.second.begin(), entry.second.end(), target_atom) != entry.second.end()) {
                intersects = true;
                break;
            }
        }
        if (intersects) {
            matched_compare_groups.push_back(entry.first);
        }
    }
    std::sort(matched_compare_groups.begin(), matched_compare_groups.end());

    std::vector<int> matched_primary_groups;
    for (int compare_group_id : matched_compare_groups) {
        for (const auto &pair : manual_shift_pairs_) {
            if (pair.compare_group_id == compare_group_id) {
                matched_primary_groups.push_back(pair.primary_group_id);
            }
        }
    }
    std::sort(matched_primary_groups.begin(), matched_primary_groups.end());
    matched_primary_groups.erase(
        std::unique(matched_primary_groups.begin(), matched_primary_groups.end()),
        matched_primary_groups.end());
    if (!matched_primary_groups.empty()) {
        spectrum_widget_->set_selected_groups(matched_primary_groups);
    }

    std::ostringstream status;
    status << "Selected compare atom " << atom_index;
    if (!matched_compare_groups.empty()) {
        status << " | cmp groups:";
        for (std::size_t i = 0; i < matched_compare_groups.size(); ++i) {
            status << matched_compare_groups[i];
            if (i + 1 < matched_compare_groups.size()) {
                status << ",";
            }
        }
    }
    if (!matched_primary_groups.empty()) {
        status << " | mapped A groups:";
        for (std::size_t i = 0; i < matched_primary_groups.size(); ++i) {
            status << matched_primary_groups[i];
            if (i + 1 < matched_primary_groups.size()) {
                status << ",";
            }
        }
    }
    status_box_->copy_label(status.str().c_str());
}

void AppWindow::apply_nucleus_visuals(const QueuedJob &job, const std::string &nucleus) {
    active_nucleus_ = normalize_nucleus_label(nucleus);
    if (active_nucleus_.empty()) {
        active_nucleus_ = "1H";
    }
    spectrum_widget_->set_render_settings(job.config.line_shape, job.config.fwhm_hz, job.config.frequency_mhz);
    apply_active_experimental_overlay();

    SpectralProductFiles files;
    auto it = active_spectral_products_.find(active_nucleus_);
    if (it != active_spectral_products_.end()) {
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
    spectrum_widget_->set_nucleus_label(spectrum_title_for_label(active_nucleus_));

    peak_browser_->clear();
    atom_browser_->clear();
    group_to_peak_row_.clear();
    peak_row_to_group_.clear();
    atom_to_row_.clear();
    atom_row_to_atom_.clear();
    atom_row_base_label_.clear();
    primary_group_to_ppm_.clear();

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
        marker.multiplicity = row.multiplicity;
        marker.j_hz = row.j_hz;
        markers.push_back(marker);
        primary_group_to_ppm_[row.group_id] = row.center_ppm;
    }
    spectrum_widget_->set_peak_markers(std::move(markers));

    // For compare workflow jobs: auto-overlay the product spectrum when a
    // "(Reactant)" label is active, so both traces are shown side-by-side.
    if (job.config.workflow_kind == WorkflowKind::Compare && comparison_job_index_ < 0) {
        static const std::string kReactantSuffix = " (Reactant)";
        const bool is_reactant_label = (active_nucleus_.size() > kReactantSuffix.size()
            && active_nucleus_.substr(active_nucleus_.size() - kReactantSuffix.size()) == kReactantSuffix);
        if (is_reactant_label) {
            const std::string base_nucleus = active_nucleus_.substr(0, active_nucleus_.size() - kReactantSuffix.size());
            const std::string product_label = base_nucleus + " (Product)";
            auto product_it = active_spectral_products_.find(product_label);
            if (product_it != active_spectral_products_.end() && !product_it->second.spectrum_csv.empty()) {
                const auto product_points = load_spectrum_csv(product_it->second.spectrum_csv);
                spectrum_widget_->set_comparison_points(product_points);
                const auto product_peak_rows = read_peak_rows(product_it->second.peaks_csv);
                comparison_peak_rows_.clear();
                comparison_group_to_ppm_.clear();
                comparison_group_to_atoms_ = read_assignments(product_it->second.assignments_csv);
                std::vector<ComparisonPeakMarker> cmp_markers;
                for (const auto &row : product_peak_rows) {
                    ComparisonPeakMarker m;
                    m.group_id = row.group_id;
                    m.center_ppm = row.center_ppm;
                    m.label = row.multiplicity;
                    cmp_markers.push_back(m);
                    comparison_peak_rows_.push_back(m);
                    comparison_group_to_ppm_[row.group_id] = row.center_ppm;
                }
                spectrum_widget_->set_comparison_peak_markers(std::move(cmp_markers));
            } else {
                spectrum_widget_->clear_comparison_points();
            }
        } else {
            spectrum_widget_->clear_comparison_points();
        }
    } else {
        apply_comparison_visuals(job);
    }

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
        selected_job_index_ = -1;
        return;
    }
    selected_job_index_ = selected - 1;

    QueuedJob job;
    {
        std::lock_guard<std::mutex> lock(jobs_mutex_);
        const std::size_t index = static_cast<std::size_t>(selected - 1);
        if (index >= jobs_.size()) {
            return;
        }
        job = jobs_[index];
    }

    if (workflow_choice_ != nullptr) {
        int idx = 0;
        for (int i = 0; i < workflow_choice_->size(); ++i) {
            const char *label = workflow_choice_->text(i);
            if (label != nullptr && workflow_kind_from_string(label) == job.config.workflow_kind) {
                idx = i;
                break;
            }
        }
        workflow_choice_->value(idx);
    }

    auto structure_atoms = read_structure_atoms(job.structure_atoms_csv);
    auto structure_bonds = read_structure_bonds(job.structure_bonds_csv);
    if (structure_widget_ != nullptr) {
        structure_widget_->set_structure(std::move(structure_atoms), std::move(structure_bonds));
    }

    if (job.config.workflow_kind == WorkflowKind::Compare && !job.structure_atoms_product_csv.empty()) {
        const auto product_atoms = read_structure_atoms(job.structure_atoms_product_csv);
        const auto product_bonds = read_structure_bonds(job.structure_bonds_product_csv);
        if (compare_structure_widget_ != nullptr) {
            compare_structure_widget_->set_structure(product_atoms, product_bonds);
        }
        if (structure_current_label_ != nullptr) {
            structure_current_label_->copy_label("A: Reactant");
        }
        if (structure_compare_label_ != nullptr) {
            structure_compare_label_->copy_label("B: Product");
        }
        update_structure_compare_layout(true);
    } else if (comparison_job_index_ < 0) {
        if (compare_structure_widget_ != nullptr) {
            compare_structure_widget_->set_structure({}, {});
        }
        update_structure_compare_layout(false);
    }

    active_spectral_products_ = spectral_products_for_job(job);

    if (spectrum_nucleus_choice_ != nullptr) {
        spectrum_nucleus_choice_->clear();

        std::vector<std::string> ordered_labels;
        const std::vector<std::string> preferred = {"1H", "13C", "19F", "31P"};
        for (const auto &label : preferred) {
            auto entry = active_spectral_products_.find(label);
            if (entry != active_spectral_products_.end() && !entry->second.spectrum_csv.empty()) {
                ordered_labels.push_back(label);
            }
        }
        for (const auto &entry : active_spectral_products_) {
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

void AppWindow::refresh_example_choice() {
    if (example_choice_ == nullptr) {
        return;
    }

    example_choice_->clear();
    example_case_rows_.clear();
    example_bundle_names_.clear();
    example_bundle_row_indices_.clear();
    example_choice_case_indices_.clear();

    example_choice_->add("Examples: choose molecule");
    example_choice_case_indices_.push_back(-1);

    const std::string cases_csv = find_examples_case_csv();
    const auto cases = load_example_cases(cases_csv);

    struct ExampleBundleRow {
        std::string key;
        std::string display_name;
        std::vector<int> row_indices;
    };

    std::vector<ExampleBundleRow> bundles;
    std::unordered_map<std::string, int> bundle_key_to_idx;

    for (const auto &example : cases) {
        example_case_rows_.push_back(example.csv_row);
        const int row_index = static_cast<int>(example_case_rows_.size()) - 1;

        std::string bundle_key = trim_copy(example.smiles);
        if (bundle_key.empty()) {
            bundle_key = lowercase_copy(example.case_name);
        }
        if (bundle_key.empty()) {
            continue;
        }

        auto it = bundle_key_to_idx.find(bundle_key);
        if (it == bundle_key_to_idx.end()) {
            ExampleBundleRow bundle;
            bundle.key = bundle_key;
            bundle.display_name = infer_example_bundle_name(example);
            if (bundle.display_name.empty()) {
                bundle.display_name = "Example";
            }
            bundle.row_indices.push_back(row_index);
            const int new_idx = static_cast<int>(bundles.size());
            bundles.push_back(std::move(bundle));
            bundle_key_to_idx[bundle_key] = new_idx;
        } else {
            bundles[static_cast<std::size_t>(it->second)].row_indices.push_back(row_index);
        }
    }

    std::sort(bundles.begin(), bundles.end(), [](const ExampleBundleRow &a, const ExampleBundleRow &b) {
        return lowercase_copy(a.display_name) < lowercase_copy(b.display_name);
    });

    for (const auto &bundle : bundles) {
        example_bundle_names_.push_back(bundle.display_name);
        example_bundle_row_indices_.push_back(bundle.row_indices);
        const int bundle_index = static_cast<int>(example_bundle_names_.size()) - 1;
        example_choice_case_indices_.push_back(bundle_index);
        example_choice_->add(bundle.display_name.c_str());
    }

    example_choice_->value(0);

    if (example_bundle_names_.empty()) {
        if (!cases_csv.empty()) {
            status_box_->copy_label(("No examples found in " + cases_csv).c_str());
        } else {
            status_box_->label("No examples file found (tests/spectra_comparison_cases.csv)");
        }
        if (load_example_calc_button_ != nullptr) {
            load_example_calc_button_->deactivate();
        }
        if (load_example_bundle_button_ != nullptr) {
            load_example_bundle_button_->deactivate();
        }
        return;
    }

    if (load_example_calc_button_ != nullptr) {
        load_example_calc_button_->activate();
    }
    if (load_example_bundle_button_ != nullptr) {
        load_example_bundle_button_->activate();
    }
}

void AppWindow::refresh_experimental_choice() {
    if (experimental_choice_ == nullptr) {
        return;
    }

    experimental_choice_->clear();
    experimental_choice_keys_.clear();

    experimental_choice_->add("Exp: none");
    experimental_choice_keys_.push_back("");

    for (const auto &entry : experimental_overlays_) {
        std::string item = "Exp: " + entry.first;
        experimental_choice_->add(item.c_str());
        experimental_choice_keys_.push_back(entry.first);
    }

    int selected_index = 0;
    if (!active_experimental_overlay_key_.empty()) {
        for (int i = 1; i < static_cast<int>(experimental_choice_keys_.size()); ++i) {
            if (experimental_choice_keys_[static_cast<std::size_t>(i)] == active_experimental_overlay_key_) {
                selected_index = i;
                break;
            }
        }
    }

    if (selected_index == 0 && !experimental_overlays_.empty() && !active_experimental_overlay_key_.empty()) {
        active_experimental_overlay_key_.clear();
    }
    experimental_choice_->value(selected_index);
    experimental_choice_->activate();
}

void AppWindow::maybe_apply_example_overlay_for_active_selection() {
    if (selected_job_index_ < 0) {
        return;
    }
    auto by_job = example_job_overlay_keys_.find(selected_job_index_);
    if (by_job == example_job_overlay_keys_.end()) {
        return;
    }

    const std::string nucleus = normalize_nucleus_label(active_nucleus_);
    auto by_nucleus = by_job->second.find(nucleus);
    if (by_nucleus == by_job->second.end()) {
        active_experimental_overlay_key_.clear();
        refresh_experimental_choice();
        apply_active_experimental_overlay();
        return;
    }

    if (experimental_overlays_.find(by_nucleus->second) == experimental_overlays_.end()) {
        active_experimental_overlay_key_.clear();
        refresh_experimental_choice();
        apply_active_experimental_overlay();
        return;
    }

    active_experimental_overlay_key_ = by_nucleus->second;
    refresh_experimental_choice();
    apply_active_experimental_overlay();
}

void AppWindow::apply_active_experimental_overlay() {
    if (spectrum_widget_ == nullptr) {
        return;
    }
    if (active_experimental_overlay_key_.empty()) {
        spectrum_widget_->clear_experimental_points();
        return;
    }

    const auto it = experimental_overlays_.find(active_experimental_overlay_key_);
    if (it == experimental_overlays_.end()) {
        spectrum_widget_->clear_experimental_points();
        return;
    }
    spectrum_widget_->set_experimental_points(it->second);
}

void AppWindow::refresh_workflow_browser(const QueuedJob *job) {
    if (workflow_info_line1_ == nullptr || workflow_info_line2_ == nullptr || workflow_info_line3_ == nullptr) {
        return;
    }

    if (job == nullptr) {
        if (workflow_progress_widget_ != nullptr) {
            workflow_progress_widget_->set_progress_state(-1, -1, false, false, false, 0.0, "Workflow idle");
        }
        workflow_info_line1_->copy_label("No active job");
        workflow_info_line2_->copy_label("Queue a job, then click Run Pending / Run Selected.");
        workflow_info_line3_->copy_label("Workflow details will appear here while calculations are running.");
        return;
    }

    std::ostringstream header;
    header << "Job " << job->id << " | " << truncate_text(job->config.job_name, 22)
           << " | " << to_string(job->config.input_format) << " | " << workflow_display_name(job->config.workflow_kind);
    if (job->config.workflow_kind == WorkflowKind::Nmr) {
        header << ":" << job->config.nucleus;
    }
    workflow_info_line1_->copy_label(truncate_text(header.str(), 116).c_str());

    std::ostringstream current;
    if (!job->progress_message.empty()) {
        current << "Current: " << truncate_text(job->progress_message, 108);
    } else if (job->status == "done") {
        current << "Current: Prediction complete";
    } else if (job->status == "failed") {
        current << "Current: Failed - " << truncate_text(job->message, 88);
    } else if (job->status == "running") {
        current << "Current: Running...";
    } else {
        current << "Current: Waiting to run";
    }
    workflow_info_line2_->copy_label(current.str().c_str());

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

    std::ostringstream detail;
    if (active_idx >= 0 && active_idx < static_cast<int>(steps.size())) {
        detail << "Step " << (active_idx + 1) << "/" << steps.size()
               << ": " << steps[static_cast<std::size_t>(active_idx)].label
               << " | " << steps[static_cast<std::size_t>(active_idx)].method;
    } else if (job->status == "done") {
        if (job->config.workflow_kind == WorkflowKind::Compare && !job->reaction_summary_json.empty()) {
            // Extract overall indicator from JSON for a one-line summary.
            const std::string &rsj = job->reaction_summary_json;
            auto extract = [&](const std::string &key) -> std::string {
                const std::string search = "\"" + key + "\": \"";
                auto pos = rsj.find(search);
                if (pos == std::string::npos) {
                    return {};
                }
                pos += search.size();
                auto end = rsj.find('"', pos);
                return end == std::string::npos ? rsj.substr(pos) : rsj.substr(pos, end - pos);
            };
            const std::string indicator = extract("overall");
            const std::string description = extract("description");
            if (!indicator.empty()) {
                detail << "Result: " << indicator;
                if (!description.empty()) {
                    detail << " | " << description;
                }
            } else {
                detail << "Comparison complete. See output directory for reaction_summary.json.";
            }
        } else {
            detail << "Completed all workflow steps successfully.";
        }
    } else if (job->status == "failed") {
        detail << "Workflow ended with an error.";
    } else {
        detail << "Waiting for workflow start.";
    }
    if (job->status == "running" && job->progress_fraction > 0.0) {
        const int pct = static_cast<int>(job->progress_fraction * 100.0 + 0.5);
        detail << " (" << pct << "%)";
    }
    workflow_info_line3_->copy_label(truncate_text(detail.str(), 120).c_str());
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
                    if (jobs_[next_index].config.workflow_kind == WorkflowKind::Cd) {
                        jobs_[next_index].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> CD bands";
                    } else if (jobs_[next_index].config.workflow_kind == WorkflowKind::All) {
                        jobs_[next_index].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> NMR + CD products";
                    } else if (jobs_[next_index].config.workflow_kind == WorkflowKind::Compare) {
                        jobs_[next_index].message = "Reactant + product conformers -> NMR prediction -> spectral comparison";
                    } else {
                        const std::string nucleus_mode = (normalize_nucleus_label(jobs_[next_index].config.nucleus) == "1H"
                                                          && jobs_[next_index].config.nucleus != "auto")
                                                             ? "1H model"
                                                             : "multi-nucleus (1H/13C/19F/31P as present)";
                        jobs_[next_index].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> " + nucleus_mode;
                    }
                    jobs_[next_index].progress_stage = "launch";
                    jobs_[next_index].progress_message = "Launching backend process";
                    jobs_[next_index].progress_fraction = 0.0;
                }
            } else {
                for (std::size_t i = 0; i < jobs_.size(); ++i) {
                    if (jobs_[i].status == "pending") {
                        next_index = i;
                        jobs_[i].status = "running";
                        if (jobs_[i].config.workflow_kind == WorkflowKind::Cd) {
                            jobs_[i].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> CD bands";
                        } else if (jobs_[i].config.workflow_kind == WorkflowKind::All) {
                            jobs_[i].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> NMR + CD products";
                        } else if (jobs_[i].config.workflow_kind == WorkflowKind::Compare) {
                            jobs_[i].message = "Reactant + product conformers -> NMR prediction -> spectral comparison";
                        } else {
                            const std::string nucleus_mode = (normalize_nucleus_label(jobs_[i].config.nucleus) == "1H"
                                                              && jobs_[i].config.nucleus != "auto")
                                                                 ? "1H model"
                                                                 : "multi-nucleus (1H/13C/19F/31P as present)";
                            jobs_[i].message = "ETKDG conformers -> xTB/MMFF -> Boltzmann -> " + nucleus_mode;
                        }
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
                job.reaction_summary_json = result.reaction_summary_json;
                job.structure_product_svg = result.structure_product_svg;
                job.structure_atoms_product_csv = result.structure_atoms_product_csv;
                job.structure_bonds_product_csv = result.structure_bonds_product_csv;
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
    const int selected_index = (selected > 0) ? (selected - 1) : -1;
    int primary_index = selected_job_index_;
    if (primary_index < 0) {
        primary_index = selected_index;
    }
    queue_browser_->clear();

    for (std::size_t i = 0; i < jobs_.size(); ++i) {
        const int row_index = static_cast<int>(i);
        const bool is_selected = (row_index == selected_index);
        const bool is_primary_role = (row_index == primary_index);
        const bool is_compare_role = (row_index == comparison_job_index_);
        queue_browser_->add(format_label_for(jobs_[i], is_selected, is_primary_role, is_compare_role).c_str());
        const int browser_row = queue_browser_->size();
        if (is_primary_role) {
            queue_browser_->icon(browser_row, &kRoleBlueDot);
        } else if (is_compare_role) {
            queue_browser_->icon(browser_row, &kRoleGreenDot);
        } else {
            queue_browser_->icon(browser_row, nullptr);
        }
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
    if (status_box_ != nullptr) {
        const int max_chars = std::max(18, (status_box_->w() - 12) / 7);
        status_box_->copy_label(truncate_text(status.str(), static_cast<std::size_t>(max_chars)).c_str());
    }

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
