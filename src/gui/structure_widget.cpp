#include "gui/structure_widget.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace easynmr {
namespace {

Fl_Color rgb(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

Fl_Color element_color(const std::string &symbol) {
    if (symbol == "O") {
        return rgb(225, 114, 114);
    }
    if (symbol == "N") {
        return rgb(121, 141, 220);
    }
    if (symbol == "S") {
        return rgb(222, 177, 113);
    }
    if (symbol == "P") {
        return rgb(216, 154, 105);
    }
    if (symbol == "F" || symbol == "Cl" || symbol == "Br" || symbol == "I") {
        return rgb(112, 178, 143);
    }
    if (symbol == "H") {
        return rgb(210, 216, 226);
    }
    return rgb(98, 106, 122);
}

bool has_overlap(const std::vector<int> &attached_hydrogens, const std::unordered_set<int> &highlighted) {
    for (int h : attached_hydrogens) {
        if (highlighted.count(h) > 0) {
            return true;
        }
    }
    return false;
}

void draw_wedge_bond(int ax, int ay, int bx, int by, Fl_Color color) {
    const double dx = static_cast<double>(bx - ax);
    const double dy = static_cast<double>(by - ay);
    const double len = std::max(1.0, std::hypot(dx, dy));
    const double nx = -dy / len;
    const double ny = dx / len;
    const double half_width = 5.5;

    const int x1 = ax;
    const int y1 = ay;
    const int x2 = static_cast<int>(bx + nx * half_width);
    const int y2 = static_cast<int>(by + ny * half_width);
    const int x3 = static_cast<int>(bx - nx * half_width);
    const int y3 = static_cast<int>(by - ny * half_width);

    fl_color(color);
    fl_begin_polygon();
    fl_vertex(static_cast<double>(x1), static_cast<double>(y1));
    fl_vertex(static_cast<double>(x2), static_cast<double>(y2));
    fl_vertex(static_cast<double>(x3), static_cast<double>(y3));
    fl_end_polygon();
}

void draw_hashed_wedge_bond(int ax, int ay, int bx, int by, Fl_Color color) {
    const double dx = static_cast<double>(bx - ax);
    const double dy = static_cast<double>(by - ay);
    const double len = std::max(1.0, std::hypot(dx, dy));
    const double nx = -dy / len;
    const double ny = dx / len;
    constexpr int kStripes = 7;
    const double max_half_width = 5.5;

    fl_color(color);
    fl_line_style(FL_SOLID, 1);
    for (int i = 1; i <= kStripes; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(kStripes + 1);
        const double cx = ax + dx * t;
        const double cy = ay + dy * t;
        const double hw = max_half_width * t;
        const int x1 = static_cast<int>(cx - nx * hw);
        const int y1 = static_cast<int>(cy - ny * hw);
        const int x2 = static_cast<int>(cx + nx * hw);
        const int y2 = static_cast<int>(cy + ny * hw);
        fl_line(x1, y1, x2, y2);
    }
}

} // namespace

StructureWidget::StructureWidget(int x, int y, int w, int h, const char *label)
    : Fl_Widget(x, y, w, h, label) {
    selected_fill_color_ = rgb(173, 214, 214);
    selected_border_color_ = rgb(84, 129, 128);
    attached_fill_color_ = rgb(210, 229, 246);
}

void StructureWidget::set_structure(std::vector<StructureAtom> atoms, std::vector<StructureBond> bonds) {
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);
    selected_atom_index_ = -1;
    highlighted_hydrogens_.clear();
    redraw();
}

void StructureWidget::set_empty_message(std::string message) {
    empty_message_ = std::move(message);
    if (empty_message_.empty()) {
        empty_message_ = "No 2D structure";
    }
    if (atoms_.empty()) {
        redraw();
    }
}

void StructureWidget::set_selected_atom(int atom_index) {
    selected_atom_index_ = atom_index;
    redraw();
}

void StructureWidget::set_highlight_hydrogens(const std::vector<int> &hydrogens) {
    highlighted_hydrogens_.clear();
    for (int h : hydrogens) {
        highlighted_hydrogens_.insert(h);
    }
    redraw();
}

void StructureWidget::set_highlight_palette(Fl_Color selected_fill, Fl_Color selected_border, Fl_Color attached_fill) {
    selected_fill_color_ = selected_fill;
    selected_border_color_ = selected_border;
    attached_fill_color_ = attached_fill;
    redraw();
}

void StructureWidget::clear_highlight() {
    highlighted_hydrogens_.clear();
    selected_atom_index_ = -1;
    redraw();
}

void StructureWidget::set_atom_overlay(std::map<int, double> values, std::string mode) {
    atom_overlay_values_ = std::move(values);
    atom_overlay_mode_ = std::move(mode);
    redraw();
}

void StructureWidget::clear_atom_overlay() {
    atom_overlay_values_.clear();
    atom_overlay_mode_ = "none";
    redraw();
}

void StructureWidget::set_on_atom_selected(std::function<void(int, const std::vector<int> &)> callback) {
    on_atom_selected_ = std::move(callback);
}

std::pair<int, int> StructureWidget::to_screen(
    double px,
    double py,
    double min_x,
    double max_x,
    double min_y,
    double max_y) const {
    constexpr int pad = 22;

    const double span_x = std::max(1e-6, max_x - min_x);
    const double span_y = std::max(1e-6, max_y - min_y);
    const double usable_w = std::max(1, w() - 2 * pad);
    const double usable_h = std::max(1, h() - 2 * pad);
    const double scale = std::min(usable_w / span_x, usable_h / span_y);

    const double cx = (min_x + max_x) * 0.5;
    const double cy = (min_y + max_y) * 0.5;

    const int sx = static_cast<int>(x() + w() * 0.5 + (px - cx) * scale);
    const int sy = static_cast<int>(y() + h() * 0.5 - (py - cy) * scale);
    return {sx, sy};
}

void StructureWidget::draw() {
    fl_push_clip(x(), y(), w(), h());

    fl_color(rgb(248, 250, 253));
    fl_rectf(x(), y(), w(), h());
    fl_color(rgb(219, 227, 238));
    fl_rect(x(), y(), w(), h());

    if (atoms_.empty()) {
        fl_color(rgb(118, 128, 142));
        fl_font(FL_HELVETICA, 12);
        fl_draw(empty_message_.c_str(), x() + 14, y() + 22);
        fl_pop_clip();
        return;
    }

    double min_x = atoms_.front().x;
    double max_x = atoms_.front().x;
    double min_y = atoms_.front().y;
    double max_y = atoms_.front().y;
    for (const auto &atom : atoms_) {
        min_x = std::min(min_x, atom.x);
        max_x = std::max(max_x, atom.x);
        min_y = std::min(min_y, atom.y);
        max_y = std::max(max_y, atom.y);
    }

    screen_positions_.clear();
    for (const auto &atom : atoms_) {
        screen_positions_[atom.atom_index] = to_screen(atom.x, atom.y, min_x, max_x, min_y, max_y);
    }

    fl_color(rgb(100, 111, 132));
    for (const auto &bond : bonds_) {
        const auto a_it = screen_positions_.find(bond.atom_a);
        const auto b_it = screen_positions_.find(bond.atom_b);
        if (a_it == screen_positions_.end() || b_it == screen_positions_.end()) {
            continue;
        }
        const auto [ax, ay] = a_it->second;
        const auto [bx, by] = b_it->second;

        int from_atom = bond.atom_a;
        int to_atom = bond.atom_b;
        if (bond.stereo_from_atom == bond.atom_b) {
            from_atom = bond.atom_b;
            to_atom = bond.atom_a;
        } else if (bond.stereo_from_atom == bond.atom_a) {
            from_atom = bond.atom_a;
            to_atom = bond.atom_b;
        }

        if (bond.stereo_style == "wedge" || bond.stereo_style == "dash") {
            const auto from_it = screen_positions_.find(from_atom);
            const auto to_it = screen_positions_.find(to_atom);
            if (from_it != screen_positions_.end() && to_it != screen_positions_.end()) {
                const auto [fx, fy] = from_it->second;
                const auto [tx, ty] = to_it->second;
                if (bond.stereo_style == "wedge") {
                    draw_wedge_bond(fx, fy, tx, ty, rgb(95, 107, 130));
                } else {
                    draw_hashed_wedge_bond(fx, fy, tx, ty, rgb(95, 107, 130));
                }
                continue;
            }
        }

        fl_line_style(FL_SOLID, 2);
        fl_line(ax, ay, bx, by);
        if (bond.order >= 2) {
            const double dx = static_cast<double>(bx - ax);
            const double dy = static_cast<double>(by - ay);
            const double len = std::max(1.0, std::hypot(dx, dy));
            const int ox = static_cast<int>(-dy / len * 3.0);
            const int oy = static_cast<int>(dx / len * 3.0);
            fl_line(ax + ox, ay + oy, bx + ox, by + oy);
        }
    }

    for (const auto &atom : atoms_) {
        const auto pos_it = screen_positions_.find(atom.atom_index);
        if (pos_it == screen_positions_.end()) {
            continue;
        }
        const auto [sx, sy] = pos_it->second;

        const bool is_selected = (atom.atom_index == selected_atom_index_);
        const bool touches_highlight = has_overlap(atom.attached_hydrogens, highlighted_hydrogens_);

        if (is_selected) {
            fl_color(selected_fill_color_);
            fl_pie(sx - 16, sy - 16, 32, 32, 0, 360);
        }
        if (touches_highlight) {
            fl_color(attached_fill_color_);
            fl_pie(sx - 13, sy - 13, 26, 26, 0, 360);
        }

        // Determine atom fill colour (overlay takes precedence over element colour).
        auto overlay_it = atom_overlay_values_.find(atom.atom_index);
        Fl_Color atom_fill;
        bool has_overlay = (atom_overlay_mode_ != "none" && overlay_it != atom_overlay_values_.end());
        if (has_overlay) {
            const double t = std::max(0.0, std::min(1.0, overlay_it->second));
            const unsigned char base_r = 245, base_g = 245, base_b = 245;
            unsigned char end_r, end_g, end_b;
            if (atom_overlay_mode_ == "gradient_red") {
                end_r = 210; end_g = 50; end_b = 50;
            } else if (atom_overlay_mode_ == "gradient_blue") {
                end_r = 50; end_g = 80; end_b = 200;
            } else {  // gradient_green (pKa)
                end_r = 50; end_g = 160; end_b = 90;
            }
            atom_fill = rgb(
                static_cast<unsigned char>(base_r + t * (end_r - base_r)),
                static_cast<unsigned char>(base_g + t * (end_g - base_g)),
                static_cast<unsigned char>(base_b + t * (end_b - base_b))
            );
        } else {
            atom_fill = element_color(atom.element);
        }

        const int radius = is_selected ? 7 : (touches_highlight ? 6 : has_overlay ? 6 : 5);
        fl_color(atom_fill);
        fl_pie(sx - radius, sy - radius, radius * 2, radius * 2, 0, 360);

        fl_color(is_selected ? selected_border_color_ : rgb(78, 86, 102));
        fl_line_style(FL_SOLID, is_selected ? 2 : 1);
        fl_arc(sx - radius, sy - radius, radius * 2, radius * 2, 0, 360);
    }

    fl_line_style(FL_SOLID, 0);
    fl_pop_clip();
}

int StructureWidget::handle(int event) {
    if (event == FL_PUSH && Fl::event_button() == FL_LEFT_MOUSE) {
        if (atoms_.empty()) {
            return Fl_Widget::handle(event);
        }

        const int mx = Fl::event_x();
        const int my = Fl::event_y();

        int best_atom = -1;
        double best_dist = 1e9;
        for (const auto &atom : atoms_) {
            const auto it = screen_positions_.find(atom.atom_index);
            if (it == screen_positions_.end()) {
                continue;
            }
            const auto [sx, sy] = it->second;
            const double dist = std::hypot(static_cast<double>(mx - sx), static_cast<double>(my - sy));
            if (dist < best_dist) {
                best_dist = dist;
                best_atom = atom.atom_index;
            }
        }

        if (best_atom > 0 && best_dist <= 18.0) {
            selected_atom_index_ = best_atom;
            redraw();

            if (on_atom_selected_) {
                for (const auto &atom : atoms_) {
                    if (atom.atom_index == best_atom) {
                        on_atom_selected_(best_atom, atom.attached_hydrogens);
                        break;
                    }
                }
            }
            return 1;
        }
    }

    return Fl_Widget::handle(event);
}

} // namespace easynmr
