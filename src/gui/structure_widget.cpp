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
        return rgb(215, 51, 51);
    }
    if (symbol == "N") {
        return rgb(56, 100, 210);
    }
    if (symbol == "S") {
        return rgb(214, 152, 20);
    }
    if (symbol == "F" || symbol == "Cl" || symbol == "Br" || symbol == "I") {
        return rgb(25, 132, 74);
    }
    return rgb(70, 75, 84);
}

bool has_overlap(const std::vector<int> &attached_hydrogens, const std::unordered_set<int> &highlighted) {
    for (int h : attached_hydrogens) {
        if (highlighted.count(h) > 0) {
            return true;
        }
    }
    return false;
}

} // namespace

StructureWidget::StructureWidget(int x, int y, int w, int h, const char *label)
    : Fl_Widget(x, y, w, h, label) {}

void StructureWidget::set_structure(std::vector<StructureAtom> atoms, std::vector<StructureBond> bonds) {
    atoms_ = std::move(atoms);
    bonds_ = std::move(bonds);
    selected_atom_index_ = -1;
    highlighted_hydrogens_.clear();
    redraw();
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

void StructureWidget::clear_highlight() {
    highlighted_hydrogens_.clear();
    selected_atom_index_ = -1;
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
        fl_draw("No 2D structure", x() + 14, y() + 22);
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

    fl_color(rgb(95, 104, 120));
    for (const auto &bond : bonds_) {
        const auto a_it = screen_positions_.find(bond.atom_a);
        const auto b_it = screen_positions_.find(bond.atom_b);
        if (a_it == screen_positions_.end() || b_it == screen_positions_.end()) {
            continue;
        }
        const auto [ax, ay] = a_it->second;
        const auto [bx, by] = b_it->second;

        fl_line_style(FL_SOLID, 1);
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

        const bool is_carbon = (atom.element == "C");
        const bool is_selected = (atom.atom_index == selected_atom_index_);
        const bool touches_highlight = has_overlap(atom.attached_hydrogens, highlighted_hydrogens_);

        if (touches_highlight) {
            fl_color(rgb(222, 237, 246));
            fl_pie(sx - 15, sy - 15, 30, 30, 0, 360);
        }

        const bool show_full_atom =
            !is_carbon || is_selected || touches_highlight;

        if (show_full_atom) {
            fl_color(rgb(255, 255, 255));
            fl_pie(sx - 11, sy - 11, 22, 22, 0, 360);

            fl_color(is_selected ? rgb(136, 175, 170) : rgb(196, 207, 221));
            fl_line_style(FL_SOLID, is_selected ? 3 : 1);
            fl_arc(sx - 11, sy - 11, 22, 22, 0, 360);

            fl_color(element_color(atom.element));
            fl_font(FL_HELVETICA_BOLD, 12);
            fl_draw(atom.element.c_str(), sx - 4, sy + 5);
        } else {
            fl_color(rgb(88, 97, 113));
            fl_pie(sx - 3, sy - 3, 6, 6, 0, 360);
        }

        std::ostringstream idx;
        idx << atom.atom_index;
        fl_color(rgb(122, 131, 146));
        fl_font(FL_HELVETICA, 9);
        fl_draw(idx.str().c_str(), sx + 8, sy - 8);

        if (!atom.attached_hydrogens.empty() && (is_selected || touches_highlight)) {
            std::ostringstream badge;
            badge << "H" << atom.attached_hydrogens.size();
            fl_color(rgb(236, 242, 249));
            fl_rectf(sx - 10, sy + 10, 20, 11);
            fl_color(rgb(190, 201, 216));
            fl_rect(sx - 10, sy + 10, 20, 11);
            fl_color(rgb(99, 108, 123));
            fl_font(FL_HELVETICA, 8);
            fl_draw(badge.str().c_str(), sx - 8, sy + 19);
        }
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
