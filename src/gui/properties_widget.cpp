#include "gui/properties_widget.h"

#include <FL/Fl.H>
#include <FL/fl_draw.H>

#include <cstdio>
#include <map>
#include <string>

namespace easynmr {

PropertiesWidget::PropertiesWidget(int x, int y, int w, int h, const char *label)
    : Fl_Group(x, y, w, h, label)
{
    begin();

    // ── Overlay toggle buttons (top strip, 28px high) ─────────────────────
    const int btn_y = y + 4;
    const int btn_h = 24;
    const int btn_w = 58;
    const int btn_gap = 4;
    int bx = x + 4;

    overlay_fplus_btn_  = new Fl_Button(bx,               btn_y, btn_w, btn_h, "f+");
    bx += btn_w + btn_gap;
    overlay_fminus_btn_ = new Fl_Button(bx,               btn_y, btn_w, btn_h, "f-");
    bx += btn_w + btn_gap;
    overlay_fzero_btn_  = new Fl_Button(bx,               btn_y, btn_w, btn_h, "f0");
    bx += btn_w + btn_gap;
    overlay_pka_btn_    = new Fl_Button(bx,               btn_y, btn_w, btn_h, "pKa");
    bx += btn_w + btn_gap;
    overlay_none_btn_   = new Fl_Button(bx,               btn_y, btn_w, btn_h, "none");

    for (Fl_Button *btn : {overlay_fplus_btn_, overlay_fminus_btn_, overlay_fzero_btn_,
                           overlay_pka_btn_, overlay_none_btn_}) {
        btn->type(FL_TOGGLE_BUTTON);
        btn->callback(overlay_cb, this);
        btn->box(FL_FLAT_BOX);
        btn->color(fl_rgb_color(210, 215, 220));
        btn->selection_color(fl_rgb_color(80, 120, 180));
        btn->labelcolor(FL_FOREGROUND_COLOR);
        btn->labelsize(11);
    }
    // "none" starts as active
    overlay_none_btn_->value(1);
    overlay_none_btn_->color(fl_rgb_color(80, 120, 180));
    overlay_none_btn_->labelcolor(FL_WHITE);

    // ── Two-column browsers ───────────────────────────────────────────────
    const int content_y = btn_y + btn_h + 4;
    const int content_h = h - (content_y - y) - 4;
    const int left_w  = (w - 12) / 2;
    const int right_w = w - 12 - left_w;

    main_browser_ = new Fl_Hold_Browser(x + 4, content_y, left_w, content_h);
    atom_browser_ = new Fl_Hold_Browser(x + 8 + left_w, content_y, right_w, content_h);

    main_browser_->textsize(12);
    atom_browser_->textsize(12);

    end();

    clear();
}

void PropertiesWidget::set_on_overlay_changed(std::function<void(std::string)> callback) {
    on_overlay_changed_ = std::move(callback);
}

void PropertiesWidget::clear() {
    if (main_browser_) main_browser_->clear();
    if (atom_browser_) {
        atom_browser_->clear();
        atom_browser_->add("@BAtom  f+     f-     f0     pKa");
    }
}

void PropertiesWidget::set_active_overlay_button(Fl_Button *active_btn, const std::string &mode) {
    for (Fl_Button *btn : {overlay_fplus_btn_, overlay_fminus_btn_, overlay_fzero_btn_,
                           overlay_pka_btn_, overlay_none_btn_}) {
        if (btn == active_btn) {
            btn->value(1);
            btn->color(fl_rgb_color(80, 120, 180));
            btn->labelcolor(FL_WHITE);
        } else {
            btn->value(0);
            btn->color(fl_rgb_color(210, 215, 220));
            btn->labelcolor(FL_FOREGROUND_COLOR);
        }
        btn->redraw();
    }
    active_overlay_ = mode;
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
        main_browser_->add(("  Formula:    " + props.formula).c_str());
    }
    main_browser_->add(("  MW:         " + fmt_d(props.mw, 2) + " g/mol").c_str());
    main_browser_->add(("  Exact MW:   " + fmt_d(props.exact_mw, 4) + " g/mol").c_str());
    main_browser_->add(("  LogP:       " + fmt_d(props.logp, 2)).c_str());
    main_browser_->add(("  TPSA:       " + fmt_d(props.tpsa, 1) + " \xc3\x85\xc2\xb2").c_str());

    char buf[128];
    std::snprintf(buf, sizeof(buf), "  HBD / HBA:  %d / %d", props.hbd, props.hba);
    main_browser_->add(buf);
    std::snprintf(buf, sizeof(buf), "  RotBonds:   %d  ArRings: %d", props.rotbonds, props.ar_rings);
    main_browser_->add(buf);

    if (props.has_electronic) {
        main_browser_->add(" ");
        main_browser_->add("@b@BElectronic (xTB GFN2)");
        main_browser_->add(("  HOMO:       " + fmt_d(props.homo_ev, 2) + " eV").c_str());
        main_browser_->add(("  LUMO:       " + fmt_d(props.lumo_ev, 2) + " eV").c_str());
        main_browser_->add(("  Gap:        " + fmt_d(props.gap_ev, 2) + " eV").c_str());
        main_browser_->add(("  Dipole:     " + fmt_d(props.dipole_debye, 2) + " D").c_str());
    }

    if (props.has_redox) {
        main_browser_->add(" ");
        main_browser_->add("@b@BRedox Estimates (gas-phase, +-0.5 V)");
        main_browser_->add(("  E_ox:  " + fmt_d(props.e_ox_nhe, 2) + " V vs NHE  ("
                            + fmt_d(props.e_ox_fc, 2) + " V vs Fc/Fc+)").c_str());
        main_browser_->add(("  E_red: " + fmt_d(props.e_red_nhe, 2) + " V vs NHE  ("
                            + fmt_d(props.e_red_fc, 2) + " V vs Fc/Fc+)").c_str());
    }

    if (!props.pka_groups.empty()) {
        main_browser_->add(" ");
        main_browser_->add("@b@BpKa Estimates (heuristic)");
        for (const auto &pg : props.pka_groups) {
            std::snprintf(buf, sizeof(buf), "  Atom %d  %-20s ~%.1f (%s)",
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
    // Build pKa lookup by 1-based atom index
    std::map<int, const PkaGroup *> pka_by_atom;
    for (const auto &pg : props.pka_groups) {
        pka_by_atom[pg.atom_index] = &pg;
    }

    if (props.has_fukui) {
        // Only show heavy atoms (filter H for cleaner display)
        for (const auto &fa : props.fukui) {
            if (fa.element == "H") continue;
            auto pka_it = pka_by_atom.find(fa.atom_index);
            std::string pka_str = "  —  ";
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
        }
    } else if (!props.pka_groups.empty()) {
        // No Fukui, but we still have pKa
        for (const auto &pg : props.pka_groups) {
            std::snprintf(buf, sizeof(buf), "Atom %d  %-20s ~%.1f %s",
                          pg.atom_index, pg.group_name.c_str(),
                          pg.pka_est, pg.site_type.c_str());
            atom_browser_->add(buf);
        }
    } else {
        atom_browser_->add("  (no xTB data available)");
    }
}

} // namespace easynmr
