#include "gui/workflow_progress_widget.h"

#include <FL/fl_draw.H>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace easynmr {
namespace {

Fl_Color rgb(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

std::string truncate_copy(const std::string &text, std::size_t max_len) {
    if (text.size() <= max_len) {
        return text;
    }
    return text.substr(0, max_len - 1) + "...";
}

std::string truncate_to_width(const std::string &text, int max_px, int font_face, int font_size) {
    if (max_px <= 0) {
        return "";
    }
    fl_font(font_face, font_size);
    int tw = 0;
    int th = 0;
    fl_measure(text.c_str(), tw, th, 0);
    if (tw <= max_px) {
        return text;
    }
    std::string trimmed = text;
    while (!trimmed.empty()) {
        trimmed.pop_back();
        const std::string candidate = trimmed + "...";
        fl_measure(candidate.c_str(), tw, th, 0);
        if (tw <= max_px) {
            return candidate;
        }
    }
    return "...";
}

void draw_pill(int px, int py, int pw, int ph, Fl_Color color) {
    if (pw <= 0 || ph <= 0) {
        return;
    }
    const int diameter = std::min(ph, pw);
    const int radius = diameter / 2;
    fl_color(color);
    if (pw <= diameter) {
        fl_pie(px, py, pw, ph, 0, 360);
        return;
    }
    fl_rectf(px + radius, py, pw - 2 * radius, ph);
    fl_pie(px, py, diameter, diameter, 90, 270);
    fl_pie(px + pw - diameter, py, diameter, diameter, -90, 90);
}

void draw_pill_border(int px, int py, int pw, int ph, Fl_Color color) {
    if (pw <= 0 || ph <= 0) {
        return;
    }
    const int diameter = std::min(ph, pw);
    const int radius = diameter / 2;
    fl_color(color);
    if (pw <= diameter) {
        fl_arc(px, py, pw, ph, 0, 360);
        return;
    }
    fl_line(px + radius, py, px + pw - radius, py);
    fl_line(px + radius, py + ph, px + pw - radius, py + ph);
    fl_arc(px, py, diameter, diameter, 90, 270);
    fl_arc(px + pw - diameter, py, diameter, diameter, -90, 90);
}

} // namespace

WorkflowProgressWidget::WorkflowProgressWidget(int x, int y, int w, int h, const char *label)
    : Fl_Widget(x, y, w, h, label) {}

void WorkflowProgressWidget::set_steps(std::vector<std::string> steps) {
    steps_ = std::move(steps);
    redraw();
}

void WorkflowProgressWidget::set_progress_state(
    int active_step_index,
    int failed_step_index,
    bool running,
    bool done,
    bool failed,
    double fraction,
    std::string status_text) {
    active_step_index_ = active_step_index;
    failed_step_index_ = failed_step_index;
    running_ = running;
    done_ = done;
    failed_ = failed;
    fraction_ = std::max(0.0, std::min(1.0, fraction));
    status_text_ = std::move(status_text);
    redraw();
}

void WorkflowProgressWidget::draw() {
    fl_push_clip(x(), y(), w(), h());

    fl_color(rgb(249, 252, 255));
    fl_rectf(x(), y(), w(), h());
    fl_color(rgb(226, 234, 245));
    fl_rect(x(), y(), w(), h());

    const int n_steps = static_cast<int>(steps_.size());
    const int bar_x = x() + 14;
    const int bar_y = y() + 5;
    const int bar_w = std::max(20, w() - 28);
    const int bar_h = 10;

    if (n_steps <= 0) {
        fl_color(rgb(130, 140, 156));
        fl_font(FL_HELVETICA, 11);
        fl_draw("No workflow steps", x() + 10, y() + h() / 2 + 4);
        fl_pop_clip();
        return;
    }

    draw_pill(bar_x, bar_y, bar_w, bar_h, rgb(228, 236, 246));
    draw_pill_border(bar_x, bar_y, bar_w, bar_h, rgb(204, 217, 234));

    double shown_fraction = 0.0;
    if (done_) {
        shown_fraction = 1.0;
    } else if (running_ || failed_) {
        shown_fraction = fraction_;
    }
    shown_fraction = std::max(0.0, std::min(1.0, shown_fraction));

    const int progress_w = static_cast<int>(std::round(shown_fraction * bar_w));
    if (progress_w > 0) {
        const Fl_Color progress_color = failed_ ? rgb(234, 187, 198) : rgb(174, 214, 205);
        draw_pill(bar_x, bar_y, std::min(progress_w, bar_w), bar_h, progress_color);
    }

    const int node_y = bar_y + bar_h / 2;
    const int node_r = 4;
    for (int i = 0; i < n_steps; ++i) {
        const double t = (n_steps <= 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n_steps - 1);
        const int node_x = bar_x + static_cast<int>(std::round(t * bar_w));

        bool node_done = false;
        if (done_) {
            node_done = true;
        } else if (failed_ && failed_step_index_ >= 0 && i < failed_step_index_) {
            node_done = true;
        } else if (running_ && active_step_index_ >= 0 && i < active_step_index_) {
            node_done = true;
        } else if ((node_x - bar_x) <= progress_w - 2) {
            node_done = true;
        }

        const bool node_active = running_ && active_step_index_ >= 0 && i == active_step_index_;
        const bool node_failed = failed_ && failed_step_index_ >= 0 && i == failed_step_index_;

        Fl_Color fill = rgb(244, 248, 253);
        Fl_Color border = rgb(191, 204, 223);
        if (node_done) {
            fill = rgb(192, 225, 216);
            border = rgb(134, 182, 168);
        }
        if (node_active) {
            fill = rgb(201, 216, 239);
            border = rgb(137, 168, 204);
            fl_color(rgb(226, 236, 249));
            fl_pie(node_x - node_r - 2, node_y - node_r - 2, 2 * (node_r + 2), 2 * (node_r + 2), 0, 360);
        }
        if (node_failed) {
            fill = rgb(242, 205, 214);
            border = rgb(202, 150, 160);
        }

        fl_color(fill);
        fl_pie(node_x - node_r, node_y - node_r, 2 * node_r, 2 * node_r, 0, 360);
        fl_color(border);
        fl_arc(node_x - node_r, node_y - node_r, 2 * node_r, 2 * node_r, 0, 360);

        fl_color(rgb(206, 218, 235));
        fl_line(node_x, bar_y + bar_h + 1, node_x, bar_y + bar_h + 4);
    }

    if (n_steps > 1) {
        const int step_w = std::max(1, bar_w / n_steps);
        for (int i = 1; i < n_steps; ++i) {
            const int sx = bar_x + i * step_w;
            fl_color(rgb(236, 243, 251));
            fl_line(sx, bar_y + 1, sx, bar_y + bar_h - 1);
        }
    }

    std::string text = status_text_;
    if (text.empty()) {
        if (done_) {
            text = "Workflow complete";
        } else if (failed_) {
            text = "Workflow failed";
        } else if (running_) {
            text = "Workflow running";
        } else {
            text = "Workflow idle";
        }
    }
    if ((running_ || done_ || failed_) && text.find('%') == std::string::npos) {
        const int pct = static_cast<int>(std::round(shown_fraction * 100.0));
        text += " (" + std::to_string(pct) + "%)";
    }
    std::string caption;
    int caption_w = 0;
    if (running_ && active_step_index_ >= 0 && active_step_index_ < n_steps) {
        std::ostringstream step_text;
        step_text << "Step " << (active_step_index_ + 1) << "/" << n_steps << ": " << steps_[active_step_index_];
        caption = truncate_copy(step_text.str(), 64);
        fl_font(FL_HELVETICA, 10);
        int th = 0;
        fl_measure(caption.c_str(), caption_w, th, 0);
    }

    int left_text_max_px = w() - 22;
    if (!caption.empty() && caption_w + 48 < w()) {
        left_text_max_px = std::max(80, w() - caption_w - 34);
    }
    const std::string left_text = truncate_to_width(text, left_text_max_px, FL_HELVETICA, 11);

    fl_color(rgb(95, 108, 127));
    fl_font(FL_HELVETICA, 11);
    fl_draw(left_text.c_str(), x() + 10, y() + h() - 8);

    if (!caption.empty() && caption_w + 48 < w()) {
        fl_color(rgb(126, 137, 154));
        fl_font(FL_HELVETICA, 10);
        fl_draw(caption.c_str(), x() + w() - caption_w - 12, y() + h() - 8);
    }

    fl_pop_clip();
}

} // namespace easynmr
