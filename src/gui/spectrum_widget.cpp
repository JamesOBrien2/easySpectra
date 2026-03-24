#include "gui/spectrum_widget.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
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
    if (points_.size() < 2) {
        return {0.0, 1.0};
    }
    auto [min_ppm_it, max_ppm_it] = std::minmax_element(points_.begin(), points_.end(), [](const auto &a, const auto &b) {
        return a.ppm < b.ppm;
    });
    return {min_ppm_it->ppm, max_ppm_it->ppm};
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
    if (points_.size() < 2) {
        return 0.0;
    }

    const auto [min_ppm, max_ppm] = active_ppm_bounds();
    const double ppm_range = (max_ppm - min_ppm == 0.0) ? 1.0 : (max_ppm - min_ppm);

    const int inner_left = x() + kPadLeft;
    const int inner_right = x() + w() - kPadRight;
    const int inner_w = std::max(1, inner_right - inner_left);
    const double x_norm = std::max(0.0, std::min(1.0, static_cast<double>(pixel_x - inner_left) / inner_w));
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
    fl_push_clip(x(), y(), w(), h());
    fl_color(rgb(250, 252, 255));
    fl_rectf(x(), y(), w(), h());
    fl_color(rgb(221, 229, 240));
    fl_rect(x(), y(), w(), h());

    const int plot_x0 = x() + kPadLeft;
    const int plot_y0 = y() + kPadTop;
    const int plot_w = w() - kPadLeft - kPadRight;
    const int plot_h = h() - kPadTop - kPadBottom;

    fl_color(rgb(85, 95, 114));
    fl_font(FL_HELVETICA_BOLD, 12);
    fl_draw(nucleus_label_.c_str(), plot_x0, y() + 18);
    fl_color(rgb(186, 204, 225));
    fl_rectf(plot_x0, y() + 21, 138, 1);

    fl_color(rgb(236, 241, 248));
    for (int i = 0; i <= 10; ++i) {
        const int gx = plot_x0 + i * plot_w / 10;
        fl_line(gx, plot_y0, gx, plot_y0 + plot_h);
    }

    fl_color(rgb(240, 244, 250));
    for (int i = 0; i <= 4; ++i) {
        const int gy = plot_y0 + i * plot_h / 4;
        fl_line(plot_x0, gy, plot_x0 + plot_w, gy);
    }

    fl_color(rgb(208, 219, 233));
    const int baseline = plot_y0 + plot_h;
    fl_line(plot_x0, baseline, plot_x0 + plot_w, baseline);

    if (points_.size() < 2) {
        fl_color(rgb(127, 136, 152));
        fl_draw("No spectrum loaded", plot_x0 + 6, plot_y0 + 22);
        fl_pop_clip();
        return;
    }

    const auto [min_ppm, max_ppm] = active_ppm_bounds();
    const auto [min_i_it, max_i_it] = std::minmax_element(points_.begin(), points_.end(), [](const auto &a, const auto &b) {
        return a.intensity < b.intensity;
    });

    const double min_intensity = min_i_it->intensity;
    const double max_intensity = max_i_it->intensity;

    const double ppm_range = (max_ppm - min_ppm == 0.0) ? 1.0 : (max_ppm - min_ppm);
    const double intensity_range = (max_intensity - min_intensity == 0.0) ? 1.0 : (max_intensity - min_intensity);

    fl_push_clip(plot_x0, plot_y0, plot_w + 1, plot_h + 1);
    fl_color(rgb(104, 116, 138));
    fl_line_style(FL_SOLID, 1);
    int prev_x = 0;
    int prev_y = 0;
    bool has_prev = false;

    for (const auto &p : points_) {
        const double x_norm = (max_ppm - p.ppm) / ppm_range;
        const double y_norm = (p.intensity - min_intensity) / intensity_range;

        const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
        const int py = baseline - static_cast<int>(y_norm * (plot_h - 6));

        if (has_prev) {
            fl_line(prev_x, prev_y, px, py);
        }

        prev_x = px;
        prev_y = py;
        has_prev = true;
    }
    fl_line_style(FL_SOLID, 0);

    if (!peak_markers_.empty()) {
        for (const auto &marker : peak_markers_) {
            const double x_norm = (max_ppm - marker.center_ppm) / ppm_range;
            const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
            if (px < plot_x0 || px > plot_x0 + plot_w) {
                continue;
            }
            if (selected_group_ids_.count(marker.group_id) > 0) {
                fl_color(rgb(123, 183, 166));
                fl_line_style(FL_SOLID, 2);
            } else {
                fl_color(rgb(203, 214, 228));
                fl_line_style(FL_SOLID, 1);
            }
            fl_line(px, plot_y0, px, baseline);
        }
        fl_line_style(FL_SOLID, 0);
    }

    if (!reference_peaks_.empty() && highlighted_reference_index_ >= 0
        && highlighted_reference_index_ < static_cast<int>(reference_peaks_.size())) {
        const auto &ref = reference_peaks_[static_cast<std::size_t>(highlighted_reference_index_)];
        fl_line_style(FL_DASH, 1);
        fl_color(rgb(174, 160, 200));
        fl_font(FL_HELVETICA, 10);
        const double x_norm = (max_ppm - ref.ppm) / ppm_range;
        if (x_norm >= 0.0 && x_norm <= 1.0) {
            const int px = plot_x0 + static_cast<int>(x_norm * plot_w);
            fl_line(px, plot_y0, px, baseline);
            fl_draw(ref.label.c_str(), px - 34, plot_y0 + 12);
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
            fl_push_clip(left, plot_y0, right - left, plot_h);
            fl_color(rgb(178, 197, 221));
            for (int hx = left - plot_h; hx < right + plot_h; hx += 9) {
                fl_line(hx, plot_y0 + plot_h, hx + plot_h, plot_y0);
            }
            fl_pop_clip();
            fl_color(rgb(156, 178, 207));
            fl_line_style(FL_DASH, 1);
            fl_rect(left, plot_y0, right - left, plot_h);
            fl_line_style(FL_SOLID, 0);
        }
    }

    if (zoom_active_) {
        std::ostringstream zoom_note;
        zoom_note << "Zoom: " << std::fixed << std::setprecision(2) << max_ppm << " to " << min_ppm
                  << " ppm (double-click to reset)";
        fl_color(rgb(129, 125, 155));
        fl_font(FL_HELVETICA, 11);
        fl_draw(zoom_note.str().c_str(), plot_x0 + 130, y() + 18);
    }

    fl_color(rgb(123, 132, 147));
    fl_font(FL_HELVETICA, 10);
    for (int i = 0; i <= 10; ++i) {
        const int tick_x = plot_x0 + i * plot_w / 10;
        const double ppm_tick = max_ppm - (static_cast<double>(i) / 10.0) * ppm_range;
        std::ostringstream tick_text;
        tick_text << std::fixed << std::setprecision(1) << ppm_tick;
        fl_draw(tick_text.str().c_str(), tick_x - 10, y() + h() - 12);
    }

    fl_color(rgb(102, 112, 130));
    fl_font(FL_HELVETICA, 11);
    fl_draw("delta [ppm]", plot_x0 + plot_w - 72, y() + h() - 8);

    fl_pop_clip();
}

int SpectrumWidget::handle(int event) {
    if (event == FL_PUSH) {
        if (points_.size() < 2) {
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

        redraw();
        return 1;
    }

    if (event == FL_RELEASE) {
        selecting_zoom_ = false;
    }

    return Fl_Widget::handle(event);
}

} // namespace easynmr
