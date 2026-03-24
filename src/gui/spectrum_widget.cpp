#include "gui/spectrum_widget.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace easynmr {
namespace {

Fl_Color rgb(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

constexpr int kPadLeft = 20;
constexpr int kPadRight = 12;
constexpr int kPadTop = 30;
constexpr int kPadBottom = 28;

std::string to_lower_copy(const std::string &text) {
    std::string out = text;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

bool is_cd_label(const std::string &label) {
    const std::string lower = to_lower_copy(label);
    return lower.find("cd") != std::string::npos;
}

struct SpectrumPalette {
    Fl_Color panel_bg = rgb(250, 252, 255);
    Fl_Color panel_border = rgb(221, 229, 240);
    Fl_Color plot_bg = rgb(246, 250, 255);
    Fl_Color title = rgb(82, 94, 115);
    Fl_Color title_rule = rgb(196, 212, 232);
    Fl_Color grid_major = rgb(228, 236, 248);
    Fl_Color grid_minor = rgb(238, 244, 252);
    Fl_Color axis = rgb(187, 202, 222);
    Fl_Color trace_soft = rgb(166, 184, 216);
    Fl_Color trace_main = rgb(112, 129, 160);
    Fl_Color marker = rgb(206, 218, 235);
    Fl_Color marker_selected = rgb(143, 194, 186);
    Fl_Color marker_band = rgb(221, 239, 236);
    Fl_Color reference = rgb(179, 166, 208);
    Fl_Color reference_label_bg = rgb(241, 236, 250);
    Fl_Color experimental_soft = rgb(236, 192, 163);
    Fl_Color experimental_main = rgb(206, 149, 112);
    Fl_Color zoom_text = rgb(124, 121, 154);
    Fl_Color tick = rgb(117, 127, 146);
    Fl_Color axis_label = rgb(97, 107, 126);
    Fl_Color selection_fill = rgb(186, 210, 238);
    Fl_Color selection_border = rgb(140, 168, 204);
};

struct MultipletComponent {
    double offset_units = 0.0;
    double amplitude = 1.0;
};

SpectrumPalette palette_for_nucleus(const std::string &nucleus_label) {
    const std::string lower = to_lower_copy(nucleus_label);
    SpectrumPalette p;
    if (lower.find("13c") != std::string::npos) {
        p.trace_soft = rgb(224, 183, 156);
        p.trace_main = rgb(194, 146, 116);
        p.marker_selected = rgb(218, 176, 139);
        p.marker_band = rgb(249, 236, 223);
        p.reference = rgb(190, 157, 139);
        p.reference_label_bg = rgb(250, 242, 233);
        p.experimental_soft = rgb(183, 202, 232);
        p.experimental_main = rgb(130, 161, 204);
        p.selection_fill = rgb(235, 208, 188);
        p.selection_border = rgb(198, 163, 137);
        return p;
    }
    if (lower.find("19f") != std::string::npos) {
        p.trace_soft = rgb(157, 206, 188);
        p.trace_main = rgb(101, 165, 143);
        p.marker_selected = rgb(120, 192, 171);
        p.marker_band = rgb(224, 246, 239);
        p.reference = rgb(123, 176, 161);
        p.reference_label_bg = rgb(233, 248, 243);
        p.experimental_soft = rgb(219, 187, 220);
        p.experimental_main = rgb(179, 133, 182);
        p.selection_fill = rgb(188, 231, 218);
        p.selection_border = rgb(127, 188, 170);
        return p;
    }
    if (lower.find("31p") != std::string::npos) {
        p.trace_soft = rgb(210, 194, 229);
        p.trace_main = rgb(160, 137, 197);
        p.marker_selected = rgb(174, 153, 213);
        p.marker_band = rgb(241, 235, 250);
        p.reference = rgb(140, 122, 182);
        p.reference_label_bg = rgb(242, 236, 252);
        p.experimental_soft = rgb(191, 221, 184);
        p.experimental_main = rgb(140, 185, 130);
        p.selection_fill = rgb(224, 214, 244);
        p.selection_border = rgb(170, 150, 213);
        return p;
    }
    return p;
}

std::vector<MultipletComponent> multiplet_components(const std::string &multiplicity_raw) {
    const std::string mult = to_lower_copy(multiplicity_raw);
    if (mult == "doublet" || mult == "d") {
        return {{-0.5, 1.0}, {0.5, 1.0}};
    }
    if (mult == "triplet" || mult == "t") {
        return {{-1.0, 1.0}, {0.0, 2.0}, {1.0, 1.0}};
    }
    if (mult == "quartet" || mult == "q") {
        return {{-1.5, 1.0}, {-0.5, 3.0}, {0.5, 3.0}, {1.5, 1.0}};
    }
    if (mult == "quintet" || mult == "quint") {
        return {{-2.0, 1.0}, {-1.0, 4.0}, {0.0, 6.0}, {1.0, 4.0}, {2.0, 1.0}};
    }
    if (mult == "multiplet" || mult == "m") {
        return {{-1.2, 0.8}, {-0.6, 1.0}, {0.0, 1.2}, {0.6, 1.0}, {1.2, 0.8}};
    }
    return {{0.0, 1.0}};
}

double line_shape_value(const std::string &line_shape, double dx, double gamma, double sigma) {
    if (line_shape == "gaussian") {
        return std::exp(-((dx * dx) / (2.0 * sigma * sigma)));
    }
    if (line_shape == "voigt") {
        const double g = std::exp(-((dx * dx) / (2.0 * sigma * sigma)));
        const double l = (gamma * gamma) / (dx * dx + gamma * gamma);
        return 0.5 * (g + l);
    }
    return (gamma * gamma) / (dx * dx + gamma * gamma);
}

} // namespace

SpectrumWidget::SpectrumWidget(int x, int y, int w, int h, const char *label)
    : Fl_Widget(x, y, w, h, label) {}

void SpectrumWidget::set_points(std::vector<SpectrumPoint> points) {
    points_ = std::move(points);
    reset_zoom();
    selecting_zoom_ = false;
    redraw();
}

bool SpectrumWidget::load_from_csv(const std::string &csv_path) {
    auto loaded = load_spectrum_csv(csv_path);
    if (loaded.empty()) {
        return false;
    }
    set_points(std::move(loaded));
    return true;
}

void SpectrumWidget::set_peak_markers(std::vector<PeakMarker> markers) {
    peak_markers_ = std::move(markers);
    redraw();
}

void SpectrumWidget::set_selected_group(int group_id) {
    selected_group_ids_.clear();
    if (group_id > 0) {
        selected_group_ids_.insert(group_id);
    }
    redraw();
}

void SpectrumWidget::set_selected_groups(const std::vector<int> &group_ids) {
    selected_group_ids_.clear();
    for (int group_id : group_ids) {
        if (group_id > 0) {
            selected_group_ids_.insert(group_id);
        }
    }
    redraw();
}

void SpectrumWidget::set_reference_peaks(std::vector<ReferencePeak> peaks) {
    reference_peaks_ = std::move(peaks);
    if (highlighted_reference_index_ >= static_cast<int>(reference_peaks_.size())) {
        highlighted_reference_index_ = -1;
    }
    redraw();
}

void SpectrumWidget::set_highlighted_reference(int ref_index) {
    if (ref_index < 0 || ref_index >= static_cast<int>(reference_peaks_.size())) {
        highlighted_reference_index_ = -1;
    } else {
        highlighted_reference_index_ = ref_index;
    }
    redraw();
}

void SpectrumWidget::set_nucleus_label(const std::string &label) {
    nucleus_label_ = label.empty() ? "NMR Spectrum" : label;
    redraw();
}

void SpectrumWidget::set_render_settings(const std::string &line_shape, double fwhm_hz, double frequency_mhz) {
    std::string normalized = to_lower_copy(line_shape);
    if (normalized != "gaussian" && normalized != "voigt") {
        normalized = "lorentzian";
    }
    line_shape_ = normalized;
    fwhm_hz_ = (std::isfinite(fwhm_hz) && fwhm_hz > 0.0) ? fwhm_hz : 1.0;
    frequency_mhz_ = (std::isfinite(frequency_mhz) && frequency_mhz > 0.0) ? frequency_mhz : 400.0;
    redraw();
}

void SpectrumWidget::set_experimental_points(std::vector<SpectrumPoint> points) {
    experimental_points_ = std::move(points);
    redraw();
}

void SpectrumWidget::clear_experimental_points() {
    experimental_points_.clear();
    redraw();
}

void SpectrumWidget::set_on_peak_selected(std::function<void(int)> callback) {
    on_peak_selected_ = std::move(callback);
}

void SpectrumWidget::reset_zoom() {
    zoom_active_ = false;
    view_min_ppm_ = 0.0;
    view_max_ppm_ = 0.0;
    selecting_zoom_ = false;
    redraw();
}

std::pair<double, double> SpectrumWidget::data_ppm_bounds() const {
    if (points_.size() < 2 && experimental_points_.size() < 2) {
        return {0.0, 1.0};
    }
    bool has_value = false;
    double min_ppm = 0.0;
    double max_ppm = 1.0;
    auto consume = [&](const std::vector<SpectrumPoint> &pts) {
        for (const auto &p : pts) {
            if (!std::isfinite(p.ppm)) {
                continue;
            }
            if (!has_value) {
                min_ppm = p.ppm;
                max_ppm = p.ppm;
                has_value = true;
                continue;
            }
            min_ppm = std::min(min_ppm, p.ppm);
            max_ppm = std::max(max_ppm, p.ppm);
        }
    };
    consume(points_);
    consume(experimental_points_);
    if (!has_value) {
        return {0.0, 1.0};
    }
    return {min_ppm, max_ppm};
}

std::pair<double, double> SpectrumWidget::active_ppm_bounds() const {
    const auto [data_min, data_max] = data_ppm_bounds();
    if (!zoom_active_) {
        return {data_min, data_max};
    }

    const double lo = std::max(data_min, std::min(view_min_ppm_, view_max_ppm_));
    const double hi = std::min(data_max, std::max(view_min_ppm_, view_max_ppm_));
    if (hi - lo < 1e-6) {
        return {data_min, data_max};
    }
    return {lo, hi};
}

double SpectrumWidget::pixel_to_ppm(int pixel_x) const {
    if (points_.size() < 2 && experimental_points_.size() < 2) {
        return 0.0;
    }

    const auto ppm_bounds = active_ppm_bounds();
    const double min_ppm = ppm_bounds.first;
    const double max_ppm = ppm_bounds.second;
    const double ppm_range = (max_ppm - min_ppm == 0.0) ? 1.0 : (max_ppm - min_ppm);

    const int inner_left = x() + kPadLeft;
    const int inner_right = x() + w() - kPadRight;
    const int inner_w = std::max(1, inner_right - inner_left);
    const double x_norm = std::max(0.0, std::min(1.0, static_cast<double>(pixel_x - inner_left) / inner_w));
    if (is_cd_label(nucleus_label_)) {
        return min_ppm + x_norm * ppm_range;
    }
    return max_ppm - x_norm * ppm_range;
}

bool SpectrumWidget::point_in_plot(int px, int py) const {
    const int plot_x0 = x() + kPadLeft;
    const int plot_y0 = y() + kPadTop;
    const int plot_w = w() - kPadLeft - kPadRight;
    const int plot_h = h() - kPadTop - kPadBottom;
    return px >= plot_x0 && px <= plot_x0 + plot_w && py >= plot_y0 && py <= plot_y0 + plot_h;
}

void SpectrumWidget::draw() {
    const SpectrumPalette pal = palette_for_nucleus(nucleus_label_);

    fl_push_clip(x(), y(), w(), h());
    fl_color(pal.panel_bg);
    fl_rectf(x(), y(), w(), h());
    fl_color(pal.panel_border);
    fl_rect(x(), y(), w(), h());

    const int plot_x0 = x() + kPadLeft;
    const int plot_y0 = y() + kPadTop;
    const int plot_w = w() - kPadLeft - kPadRight;
    const int plot_h = h() - kPadTop - kPadBottom;
    const int plot_y_bottom = plot_y0 + plot_h;

    fl_color(pal.plot_bg);
    fl_rectf(plot_x0, plot_y0, plot_w, plot_h);

    fl_color(pal.title);
    fl_font(FL_HELVETICA_BOLD, 12);
    fl_draw(nucleus_label_.c_str(), plot_x0, y() + 18);
    fl_color(pal.title_rule);
    fl_rectf(plot_x0, y() + 21, 138, 1);

    fl_color(pal.grid_major);
    for (int i = 0; i <= 10; ++i) {
        const int gx = plot_x0 + i * plot_w / 10;
        fl_line(gx, plot_y0, gx, plot_y0 + plot_h);
    }

    fl_color(pal.grid_minor);
    for (int i = 0; i <= 4; ++i) {
        const int gy = plot_y0 + i * plot_h / 4;
        fl_line(plot_x0, gy, plot_x0 + plot_w, gy);
    }

    if (points_.size() < 2 && experimental_points_.size() < 2) {
        fl_color(pal.tick);
        fl_draw("No spectrum loaded", plot_x0 + 6, plot_y0 + 22);
        fl_pop_clip();
        return;
    }

    const auto ppm_bounds = active_ppm_bounds();
    const double min_ppm = ppm_bounds.first;
    const double max_ppm = ppm_bounds.second;
    const bool reverse_axis = !is_cd_label(nucleus_label_);
    double min_intensity = 0.0;
    double max_intensity = 1.0;
    if (!points_.empty()) {
        const auto [min_i_it, max_i_it] = std::minmax_element(points_.begin(), points_.end(), [](const auto &a, const auto &b) {
            return a.intensity < b.intensity;
        });
        min_intensity = min_i_it->intensity;
        max_intensity = max_i_it->intensity;
    }

    double exp_max_abs = 1.0;
    if (!experimental_points_.empty()) {
        exp_max_abs = 0.0;
        for (const auto &p : experimental_points_) {
            exp_max_abs = std::max(exp_max_abs, std::abs(p.intensity));
        }
        if (exp_max_abs <= 1e-9) {
            exp_max_abs = 1.0;
        }
    }

    const double ppm_range = (max_ppm - min_ppm == 0.0) ? 1.0 : (max_ppm - min_ppm);
    auto ppm_to_x_norm = [&](double x_value) {
        if (reverse_axis) {
            return (max_ppm - x_value) / ppm_range;
        }
        return (x_value - min_ppm) / ppm_range;
    };
    const double intensity_range = (max_intensity - min_intensity == 0.0) ? 1.0 : (max_intensity - min_intensity);
    const double render_min = experimental_points_.empty() ? 0.0 : -1.02;
    const double render_max = 1.02;
    const double render_range = std::max(1e-9, render_max - render_min);
    auto signal_to_py = [&](double signal_value) {
        const double y_norm = (signal_value - render_min) / render_range;
        return plot_y_bottom - static_cast<int>(y_norm * (plot_h - 6));
    };
    const int zero_line_y = signal_to_py(0.0);

    fl_color(pal.axis);
    fl_line(plot_x0, zero_line_y, plot_x0 + plot_w, zero_line_y);

    fl_push_clip(plot_x0, plot_y0, plot_w + 1, plot_h + 1);
    auto draw_trace = [&](Fl_Color color, int width) {
        fl_color(color);
        fl_line_style(FL_SOLID, width);
        int prev_x = 0;
        int prev_y = 0;
        bool has_prev = false;
        for (const auto &p : points_) {
            const double x_norm = ppm_to_x_norm(p.ppm);
            const double signal = (p.intensity - min_intensity) / intensity_range;

            const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
            const int py = signal_to_py(signal);
            if (has_prev) {
                fl_line(prev_x, prev_y, px, py);
            }
            prev_x = px;
            prev_y = py;
            has_prev = true;
        }
    };
    if (!points_.empty()) {
        draw_trace(pal.trace_soft, 3);
        draw_trace(pal.trace_main, 2);
    }

    if (!experimental_points_.empty()) {
        auto draw_experimental_trace = [&](Fl_Color color, int width) {
            fl_color(color);
            fl_line_style(FL_SOLID, width);
            int prev_x = 0;
            int prev_y = 0;
            bool has_prev = false;
            for (const auto &p : experimental_points_) {
                const double x_norm = ppm_to_x_norm(p.ppm);
                const double signal = -std::min(1.0, std::abs(p.intensity) / exp_max_abs);
                const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
                const int py = signal_to_py(signal);
                if (has_prev) {
                    fl_line(prev_x, prev_y, px, py);
                }
                prev_x = px;
                prev_y = py;
                has_prev = true;
            }
        };
        draw_experimental_trace(pal.experimental_soft, 3);
        draw_experimental_trace(pal.experimental_main, 2);
    }
    fl_line_style(FL_SOLID, 0);

    // Draw selected groups as a representative multiplet envelope so the
    // highlight follows the splitting pattern rather than a center line.
    if (!points_.empty() && !peak_markers_.empty() && !selected_group_ids_.empty()) {
        const double safe_freq_mhz = std::max(1e-6, std::abs(frequency_mhz_));
        const double fwhm_ppm = std::max(1e-6, std::abs(fwhm_hz_) / safe_freq_mhz);
        const double gamma = std::max(1e-6, fwhm_ppm / 2.0);
        const double sigma = std::max(1e-6, fwhm_ppm / 2.355);
        for (const auto &marker : peak_markers_) {
            if (selected_group_ids_.count(marker.group_id) == 0) {
                continue;
            }
            const auto components = multiplet_components(marker.multiplicity);
            double j_ppm = std::abs(marker.j_hz) / safe_freq_mhz;
            if (j_ppm <= 1e-8 && components.size() > 1) {
                j_ppm = std::max(fwhm_ppm * 0.9, ppm_range * 0.0025);
            }

            double min_component_offset = 0.0;
            double max_component_offset = 0.0;
            for (const auto &component : components) {
                const double off = component.offset_units * j_ppm;
                min_component_offset = std::min(min_component_offset, off);
                max_component_offset = std::max(max_component_offset, off);
            }

            const double core_half_span = std::max(std::abs(min_component_offset), std::abs(max_component_offset));
            const double half_window = std::max(core_half_span + 4.0 * fwhm_ppm, std::max(0.02, ppm_range * 0.01));
            const double lo = std::max(min_ppm, marker.center_ppm - half_window);
            const double hi = std::min(max_ppm, marker.center_ppm + half_window);
            if (hi - lo <= 1e-8) {
                continue;
            }

            double local_peak = min_intensity;
            for (const auto &p : points_) {
                if (p.ppm >= lo && p.ppm <= hi) {
                    local_peak = std::max(local_peak, p.intensity);
                }
            }
            if (local_peak <= min_intensity) {
                local_peak = min_intensity + 0.75 * intensity_range;
            }

            const int sample_count = std::max(64, static_cast<int>((hi - lo) / ppm_range * plot_w * 1.6));
            std::vector<std::pair<double, double>> model;
            model.reserve(static_cast<std::size_t>(sample_count + 1));
            double model_max = 0.0;
            for (int i = 0; i <= sample_count; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(sample_count);
                const double sample_ppm = hi - (hi - lo) * t;
                double sample_intensity = 0.0;
                for (const auto &component : components) {
                    const double center = marker.center_ppm + component.offset_units * j_ppm;
                    const double dx = sample_ppm - center;
                    sample_intensity += component.amplitude * line_shape_value(line_shape_, dx, gamma, sigma);
                }
                model.push_back({sample_ppm, sample_intensity});
                model_max = std::max(model_max, sample_intensity);
            }
            if (model_max <= 1e-12) {
                continue;
            }

            auto draw_model = [&](Fl_Color color, int width, double scale) {
                fl_color(color);
                fl_line_style(FL_SOLID, width);
                bool has_prev = false;
                int prev_x = 0;
                int prev_y = 0;
                for (const auto &sample : model) {
                    const double x_norm = ppm_to_x_norm(sample.first);
                    const double normalized = sample.second / model_max;
                    const double shaped_intensity = min_intensity + normalized * (local_peak - min_intensity) * scale;
                    const double signal = std::max(0.0, std::min(1.0, (shaped_intensity - min_intensity) / intensity_range));
                    const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
                    const int py = signal_to_py(signal);
                    if (has_prev) {
                        fl_line(prev_x, prev_y, px, py);
                    }
                    prev_x = px;
                    prev_y = py;
                    has_prev = true;
                }
            };

            draw_model(pal.marker_band, 7, 1.06);
            draw_model(pal.marker_selected, 3, 1.02);
            fl_line_style(FL_SOLID, 0);
        }
    }

    if (!peak_markers_.empty()) {
        for (const auto &marker : peak_markers_) {
            const double x_norm = ppm_to_x_norm(marker.center_ppm);
            const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
            if (px < plot_x0 || px > plot_x0 + plot_w) {
                continue;
            }
            if (selected_group_ids_.count(marker.group_id) > 0) {
                fl_color(pal.marker_selected);
                fl_line_style(FL_SOLID, 2);
                fl_pie(px - 4, zero_line_y - 4, 8, 8, 0, 360);
            } else {
                fl_color(pal.marker);
                fl_line_style(FL_DASH, 1);
                fl_line(px, plot_y0, px, zero_line_y);
            }
        }
        fl_line_style(FL_SOLID, 0);
    }

    if (!reference_peaks_.empty() && highlighted_reference_index_ >= 0
        && highlighted_reference_index_ < static_cast<int>(reference_peaks_.size())) {
        const auto &ref = reference_peaks_[static_cast<std::size_t>(highlighted_reference_index_)];
        fl_line_style(FL_DASH, 1);
        fl_color(pal.reference);
        fl_font(FL_HELVETICA, 10);
        const double x_norm = ppm_to_x_norm(ref.ppm);
        if (x_norm >= 0.0 && x_norm <= 1.0) {
            const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
            fl_line(px, plot_y0, px, zero_line_y);
            int lw = 0;
            int lh = 0;
            fl_measure(ref.label.c_str(), lw, lh, 0);
            const int label_w = lw + 8;
            const int label_h = 14;
            const int label_x = std::max(plot_x0 + 2, std::min(plot_x0 + plot_w - label_w - 2, px - label_w / 2));
            const int label_y = plot_y0 + 4;
            fl_color(pal.reference_label_bg);
            fl_rectf(label_x, label_y, label_w, label_h);
            fl_color(pal.reference);
            fl_rect(label_x, label_y, label_w, label_h);
            fl_draw(ref.label.c_str(), label_x + 4, label_y + 11);
        }
        fl_line_style(FL_SOLID, 0);
    }
    fl_pop_clip();

    if (selecting_zoom_) {
        const int x1 = std::max(plot_x0, std::min(plot_x0 + plot_w, drag_start_x_));
        const int x2 = std::max(plot_x0, std::min(plot_x0 + plot_w, drag_current_x_));
        const int left = std::min(x1, x2);
        const int right = std::max(x1, x2);
        if (right - left > 1) {
            fl_color(pal.selection_fill);
            fl_rectf(left, plot_y0, right - left, plot_h);
            fl_color(pal.selection_border);
            fl_line_style(FL_DASH, 2);
            fl_rect(left, plot_y0, right - left, plot_h);
            fl_line_style(FL_SOLID, 0);
        }
    }

    if (zoom_active_) {
        std::ostringstream zoom_note;
        if (reverse_axis) {
            zoom_note << "Zoom: " << std::fixed << std::setprecision(2) << max_ppm << " to " << min_ppm
                      << " ppm (double-click to reset)";
        } else {
            zoom_note << "Zoom: " << std::fixed << std::setprecision(2) << min_ppm << " to " << max_ppm
                      << " nm (double-click to reset)";
        }
        fl_color(pal.zoom_text);
        fl_font(FL_HELVETICA, 11);
        fl_draw(zoom_note.str().c_str(), plot_x0 + 130, y() + 18);
    }

    fl_color(pal.tick);
    fl_font(FL_HELVETICA, 10);
    for (int i = 0; i <= 10; ++i) {
        const int tick_x = plot_x0 + i * plot_w / 10;
        double ppm_tick = 0.0;
        if (reverse_axis) {
            ppm_tick = max_ppm - (static_cast<double>(i) / 10.0) * ppm_range;
        } else {
            ppm_tick = min_ppm + (static_cast<double>(i) / 10.0) * ppm_range;
        }
        std::ostringstream tick_text;
        const int precision = (ppm_range < 3.0) ? 2 : 1;
        tick_text << std::fixed << std::setprecision(precision) << ppm_tick;
        fl_draw(tick_text.str().c_str(), tick_x - 10, y() + h() - 12);
    }

    fl_color(pal.axis_label);
    fl_font(FL_HELVETICA, 11);
    if (reverse_axis) {
        fl_draw("delta [ppm]", plot_x0 + plot_w - 72, y() + h() - 8);
    } else {
        fl_draw("wavelength [nm]", plot_x0 + plot_w - 100, y() + h() - 8);
    }

    fl_pop_clip();
}

int SpectrumWidget::handle(int event) {
    if (event == FL_PUSH) {
        if (points_.size() < 2 && experimental_points_.size() < 2) {
            return Fl_Widget::handle(event);
        }

        const int mx = Fl::event_x();
        const int my = Fl::event_y();
        if (!point_in_plot(mx, my)) {
            return Fl_Widget::handle(event);
        }

        if (Fl::event_button() == FL_RIGHT_MOUSE) {
            reset_zoom();
            return 1;
        }

        if (Fl::event_button() == FL_LEFT_MOUSE) {
            if (Fl::event_clicks() > 0) {
                reset_zoom();
                return 1;
            }
            selecting_zoom_ = true;
            drag_start_x_ = mx;
            drag_current_x_ = mx;
            redraw();
            return 1;
        }
    }

    if (event == FL_DRAG && selecting_zoom_) {
        drag_current_x_ = Fl::event_x();
        redraw();
        return 1;
    }

    if (event == FL_RELEASE && selecting_zoom_) {
        const int release_x = Fl::event_x();
        drag_current_x_ = release_x;
        const int drag_dx = std::abs(drag_current_x_ - drag_start_x_);
        selecting_zoom_ = false;

        if (drag_dx >= 8) {
            const double ppm_a = pixel_to_ppm(drag_start_x_);
            const double ppm_b = pixel_to_ppm(drag_current_x_);
            const double lo = std::min(ppm_a, ppm_b);
            const double hi = std::max(ppm_a, ppm_b);
            if (hi - lo > 1e-6) {
                view_min_ppm_ = lo;
                view_max_ppm_ = hi;
                zoom_active_ = true;
            }
            redraw();
            return 1;
        }

        if (peak_markers_.empty()) {
            redraw();
            return 1;
        }

        const double clicked_ppm = pixel_to_ppm(release_x);
        int best_group = -1;
        double best_dist = 1e9;
        for (const auto &marker : peak_markers_) {
            const double d = std::abs(marker.center_ppm - clicked_ppm);
            if (d < best_dist) {
                best_dist = d;
                best_group = marker.group_id;
            }
        }

        const auto [min_ppm, max_ppm] = active_ppm_bounds();
        const double peak_pick_window = std::max(0.01, std::min(0.35, (max_ppm - min_ppm) * 0.03));
        if (best_group > 0 && best_dist <= peak_pick_window) {
            selected_group_ids_.clear();
            selected_group_ids_.insert(best_group);
            redraw();
            if (on_peak_selected_) {
                on_peak_selected_(best_group);
            }
            return 1;
        }

        selected_group_ids_.clear();
        if (on_peak_selected_) {
            on_peak_selected_(0);
        }
        redraw();
        return 1;
    }

    if (event == FL_RELEASE) {
        selecting_zoom_ = false;
    }

    return Fl_Widget::handle(event);
}

} // namespace easynmr
