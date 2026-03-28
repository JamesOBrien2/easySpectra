#include "gui/properties_widget.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <cstdio>
#include <cmath>
#include <map>
#include <string>

namespace easynmr {

// ── Palette helpers ───────────────────────────────────────────────────────────
static Fl_Color pw(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

// Per-overlay accent colours (matched to StructureWidget gradient end colours).
struct OverlayAccent {
    Fl_Color btn_active;    // button background when active
    unsigned char er, eg, eb;  // gradient end colour (matches structure_widget.cpp)
};

static OverlayAccent overlay_accent(const std::string &mode) {
    if (mode == "gradient_red")   return { pw(185, 65, 65),  200, 60,  60 };
    if (mode == "gradient_blue")  return { pw(60,  100, 190), 55,  85,  200 };
    if (mode == "gradient_zero")  return { pw(110, 55,  180), 130, 55,  195 };
    if (mode == "gradient_green") return { pw(50,  135, 75),  55,  150, 80 };
    return { pw(85, 105, 130), 0, 0, 0 };  // "none"
}

// ── Constructor ───────────────────────────────────────────────────────────────
PropertiesWidget::PropertiesWidget(int x, int y, int w, int h, const char *label)
    : Fl_Group(x, y, w, h, label)
{
    begin();

    // ── Overlay toggle buttons ────────────────────────────────────────────
    const int btn_y = y + 6;
    const int btn_h = 24;
    const int btn_w = 58;
    const int btn_gap = 4;
    int bx = x + 6;

    overlay_fplus_btn_  = new Fl_Button(bx, btn_y, btn_w, btn_h, "f+");
    bx += btn_w + btn_gap;
    overlay_fminus_btn_ = new Fl_Button(bx, btn_y, btn_w, btn_h, "f\xe2\x81\xbb");
    bx += btn_w + btn_gap;
    overlay_fzero_btn_  = new Fl_Button(bx, btn_y, btn_w, btn_h, "f\xe2\x81\xb0");
    bx += btn_w + btn_gap;
    overlay_pka_btn_    = new Fl_Button(bx, btn_y, btn_w, btn_h, "pKa");
    bx += btn_w + btn_gap;
    overlay_none_btn_   = new Fl_Button(bx, btn_y, btn_w, btn_h, "none");

    for (Fl_Button *btn : {overlay_fplus_btn_, overlay_fminus_btn_, overlay_fzero_btn_,
                           overlay_pka_btn_, overlay_none_btn_}) {
        btn->type(FL_TOGGLE_BUTTON);
        btn->callback(overlay_cb, this);
        btn->box(FL_ROUNDED_BOX);
        btn->color(pw(222, 230, 244));
        btn->labelcolor(pw(60, 72, 90));
        btn->labelfont(FL_HELVETICA_BOLD);
        btn->labelsize(11);
    }
    // "none" starts active — use its accent colour
    overlay_none_btn_->value(1);
    overlay_none_btn_->color(overlay_accent("none").btn_active);
    overlay_none_btn_->labelcolor(FL_WHITE);

    // ── Reserve space for gradient legend (drawn in draw(), not a widget) ─
    // legend_y = btn_y + btn_h + 5  → ~y+35
    // legend_h = 22 (8px bar + 12px label + 2px gap)
    const int legend_reserved_h = 24;
    const int content_y = btn_y + btn_h + 5 + legend_reserved_h;
    const int content_h = h - (content_y - y) - 4;

    // ── Two-column area ───────────────────────────────────────────────────
    const int left_w  = (w - 13) * 3 / 5;   // 60% for descriptors
    const int right_w = w - 13 - left_w;     // 40% for atom table

    // Left: main descriptors browser
    main_browser_ = new Fl_Hold_Browser(x + 4, content_y, left_w, content_h);
    main_browser_->textsize(12);
    main_browser_->box(FL_FLAT_BOX);
    main_browser_->color(pw(249, 252, 255));
    main_browser_->selection_color(pw(200, 215, 238));
    main_browser_->textcolor(pw(50, 65, 85));
    main_browser_->when(FL_WHEN_NEVER);

    // 1px vertical divider
    const int div_x = x + 4 + left_w + 2;
    divider_ = new Fl_Box(div_x, content_y, 1, content_h);
    divider_->box(FL_FLAT_BOX);
    divider_->color(pw(210, 220, 235));

    // Right: fixed column header + atom browser
    const int abr_x = div_x + 3;
    const int col_hdr_h = 18;

    atom_col_header_ = new Fl_Box(abr_x, content_y, right_w, col_hdr_h,
                                  "Atom   f+      f-      f0      pKa");
    atom_col_header_->box(FL_FLAT_BOX);
    atom_col_header_->color(pw(218, 228, 244));
    atom_col_header_->labelsize(11);
    atom_col_header_->labelfont(FL_COURIER_BOLD);
    atom_col_header_->labelcolor(pw(55, 70, 90));
    atom_col_header_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    atom_browser_ = new Fl_Hold_Browser(abr_x, content_y + col_hdr_h,
                                        right_w, content_h - col_hdr_h);
    atom_browser_->textsize(11);
    atom_browser_->textfont(FL_COURIER);
    atom_browser_->box(FL_FLAT_BOX);
    atom_browser_->color(pw(249, 252, 255));
    atom_browser_->selection_color(pw(155, 195, 188));
    atom_browser_->textcolor(pw(50, 65, 85));
    atom_browser_->when(FL_WHEN_CHANGED);
    atom_browser_->callback(atom_browser_cb, this);

    end();
    clear();
}

// ── draw() — background fill + gradient legend ────────────────────────────────
void PropertiesWidget::draw() {
    // Background
    fl_color(pw(249, 252, 255));
    fl_rectf(x(), y(), w(), h());

    // Thin rule below buttons
    const int btn_bottom = y() + 6 + 24;  // btn_y + btn_h
    fl_color(pw(210, 222, 238));
    fl_rectf(x() + 4, btn_bottom + 4, w() - 8, 1);

    // Gradient legend strip (only when overlay is active)
    if (active_overlay_ != "none") {
        const OverlayAccent ac = overlay_accent(active_overlay_);
        const int strip_x = x() + 6;
        const int strip_w = w() - 12;
        const int strip_y = btn_bottom + 7;
        const int strip_h = 7;

        for (int i = 0; i < strip_w; ++i) {
            const double t = double(i) / strip_w;
            fl_color(fl_rgb_color(
                static_cast<unsigned char>(245 + t * (ac.er - 245)),
                static_cast<unsigned char>(245 + t * (ac.eg - 245)),
                static_cast<unsigned char>(245 + t * (ac.eb - 245))
            ));
            fl_rectf(strip_x + i, strip_y, 1, strip_h);
        }
        // Subtle border around the strip
        fl_color(pw(195, 210, 230));
        fl_rect(strip_x, strip_y, strip_w, strip_h);

        // "Low" / "High" labels
        fl_font(FL_HELVETICA, 9);
        fl_color(pw(95, 110, 130));
        fl_draw("Low", strip_x + 2, strip_y + strip_h + 10);
        const char *hi = active_overlay_ == "gradient_red"   ? "High f+" :
                         active_overlay_ == "gradient_blue"  ? "High f\xe2\x81\xbb" :
                         active_overlay_ == "gradient_zero"  ? "High f\xe2\x81\xb0" : "High pKa";
        fl_draw(hi, strip_x + strip_w - static_cast<int>(fl_width(hi)) - 2, strip_y + strip_h + 10);
    }

    Fl_Group::draw();
}

// ── Callbacks ─────────────────────────────────────────────────────────────────
void PropertiesWidget::set_on_overlay_changed(std::function<void(std::string)> callback) {
    on_overlay_changed_ = std::move(callback);
}

void PropertiesWidget::set_on_atom_selected(std::function<void(int)> callback) {
    on_atom_selected_ = std::move(callback);
}

void PropertiesWidget::clear() {
    if (main_browser_) main_browser_->clear();
    if (atom_browser_) atom_browser_->clear();
    atom_index_to_row_.clear();
    row_to_atom_index_.clear();
}

void PropertiesWidget::select_atom_row(int atom_index) {
    if (!atom_browser_) return;
    if (atom_index < 0) {
        atom_browser_->value(0);
        atom_browser_->redraw();
        return;
    }
    auto it = atom_index_to_row_.find(atom_index);
    if (it == atom_index_to_row_.end()) return;
    atom_browser_->value(it->second);
    atom_browser_->middleline(it->second);
    atom_browser_->redraw();
}

void PropertiesWidget::set_active_overlay_button(Fl_Button *active_btn, const std::string &mode) {
    struct { Fl_Button *btn; const char *mode_str; } pairs[] = {
        {overlay_fplus_btn_,  "gradient_red"},
        {overlay_fminus_btn_, "gradient_blue"},
        {overlay_fzero_btn_,  "gradient_zero"},
        {overlay_pka_btn_,    "gradient_green"},
        {overlay_none_btn_,   "none"},
    };
    for (auto &p : pairs) {
        if (p.btn == active_btn) {
            p.btn->value(1);
            p.btn->color(overlay_accent(p.mode_str).btn_active);
            p.btn->labelcolor(FL_WHITE);
        } else {
            p.btn->value(0);
            p.btn->color(pw(222, 230, 244));
            p.btn->labelcolor(pw(60, 72, 90));
        }
        p.btn->redraw();
    }
    active_overlay_ = mode;
    redraw();  // trigger gradient legend redraw
    if (on_overlay_changed_) {
        on_overlay_changed_(mode);
    }
}

void PropertiesWidget::overlay_cb(Fl_Widget *w, void *data) {
    auto *self = static_cast<PropertiesWidget *>(data);
    auto *btn = static_cast<Fl_Button *>(w);
    std::string mode = "none";
    if (btn == self->overlay_fplus_btn_)  mode = "gradient_red";
    if (btn == self->overlay_fminus_btn_) mode = "gradient_blue";
    if (btn == self->overlay_fzero_btn_)  mode = "gradient_zero";
    if (btn == self->overlay_pka_btn_)    mode = "gradient_green";
    self->set_active_overlay_button(btn, mode);
}

void PropertiesWidget::atom_browser_cb(Fl_Widget *w, void *data) {
    auto *self = static_cast<PropertiesWidget *>(data);
    const int row = self->atom_browser_->value();
    auto it = self->row_to_atom_index_.find(row);
    if (it != self->row_to_atom_index_.end() && self->on_atom_selected_) {
        self->on_atom_selected_(it->second);
    }
}

// ── Content rendering ─────────────────────────────────────────────────────────
static std::string fmt_d(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    return buf;
}

void PropertiesWidget::set_properties(const MolecularProperties &props) {
    clear();
    if (!props.valid) {
        main_browser_->add("No properties available.");
        return;
    }

    // ── Main browser: molecular descriptors ──────────────────────────────
    main_browser_->add("@b@BMolecular Descriptors");
    if (!props.formula.empty()) {
        main_browser_->add(("  Formula:   " + props.formula).c_str());
    }
    main_browser_->add(("  MW:        " + fmt_d(props.mw, 2) + " g/mol").c_str());
    main_browser_->add(("  Exact MW:  " + fmt_d(props.exact_mw, 4) + " g/mol").c_str());
    main_browser_->add(("  LogP:      " + fmt_d(props.logp, 2)).c_str());
    main_browser_->add(("  TPSA:      " + fmt_d(props.tpsa, 1) + " \xc3\x85\xc2\xb2").c_str());

    char buf[256];
    std::snprintf(buf, sizeof(buf), "  HBD / HBA: %d / %d", props.hbd, props.hba);
    main_browser_->add(buf);
    std::snprintf(buf, sizeof(buf), "  RotBonds:  %d    ArRings: %d", props.rotbonds, props.ar_rings);
    main_browser_->add(buf);

    if (props.has_electronic) {
        main_browser_->add(" ");
        main_browser_->add("@b@BElectronic (xTB GFN2)");
        main_browser_->add(("  HOMO:      " + fmt_d(props.homo_ev, 2) + " eV").c_str());
        main_browser_->add(("  LUMO:      " + fmt_d(props.lumo_ev, 2) + " eV").c_str());
        main_browser_->add(("  Gap:       " + fmt_d(props.gap_ev, 2) + " eV").c_str());
        main_browser_->add(("  Dipole:    " + fmt_d(props.dipole_debye, 2) + " D").c_str());
    }

    if (props.has_redox) {
        main_browser_->add(" ");
        main_browser_->add("@b@BRedox  (gas-phase, +-0.5 V)");
        main_browser_->add(("  E_ox:   " + fmt_d(props.e_ox_nhe, 2) + " V vs NHE  ("
                            + fmt_d(props.e_ox_fc, 2) + " V vs Fc/Fc+)").c_str());
        main_browser_->add(("  E_red:  " + fmt_d(props.e_red_nhe, 2) + " V vs NHE  ("
                            + fmt_d(props.e_red_fc, 2) + " V vs Fc/Fc+)").c_str());
    }

    if (!props.pka_groups.empty()) {
        main_browser_->add(" ");
        main_browser_->add("@b@BpKa Estimates (heuristic)");
        for (const auto &pg : props.pka_groups) {
            std::snprintf(buf, sizeof(buf), "  Atom %d  %-18s ~%.1f (%s)",
                          pg.atom_index, pg.group_name.c_str(),
                          pg.pka_est, pg.site_type.c_str());
            main_browser_->add(buf);
        }
    }

    if (!props.warnings.empty()) {
        main_browser_->add(" ");
        main_browser_->add("@b@BWarnings");
        for (const auto &w : props.warnings) {
            main_browser_->add(("  \xe2\x9a\xa0 " + w).c_str());
        }
    }

    // ── Atom browser: per-atom Fukui + pKa ───────────────────────────────
    std::map<int, const PkaGroup *> pka_by_atom;
    for (const auto &pg : props.pka_groups) {
        pka_by_atom[pg.atom_index] = &pg;
    }

    if (props.has_fukui) {
        for (const auto &fa : props.fukui) {
            if (fa.element == "H") continue;
            auto pka_it = pka_by_atom.find(fa.atom_index);
            std::string pka_str = "  --  ";
            if (pka_it != pka_by_atom.end()) {
                pka_str = "~" + fmt_d(pka_it->second->pka_est, 1)
                        + " " + pka_it->second->site_type.substr(0, 4);
            }
            std::snprintf(buf, sizeof(buf), "%s(%d)  %-6s %-6s %-6s  %s",
                          fa.element.c_str(), fa.atom_index,
                          fmt_d(fa.f_plus,  3).c_str(),
                          fmt_d(fa.f_minus, 3).c_str(),
                          fmt_d(fa.f_zero,  3).c_str(),
                          pka_str.c_str());
            atom_browser_->add(buf);
            const int row = atom_browser_->size();
            atom_index_to_row_[fa.atom_index] = row;
            row_to_atom_index_[row] = fa.atom_index;
        }
    } else if (!props.pka_groups.empty()) {
        for (const auto &pg : props.pka_groups) {
            std::snprintf(buf, sizeof(buf), "Atom %d  %-18s ~%.1f %s",
                          pg.atom_index, pg.group_name.c_str(),
                          pg.pka_est, pg.site_type.c_str());
            atom_browser_->add(buf);
            const int row = atom_browser_->size();
            atom_index_to_row_[pg.atom_index] = row;
            row_to_atom_index_[row] = pg.atom_index;
        }
    } else {
        atom_browser_->add("  (no per-atom data available)");
    }
}

} // namespace easynmr
