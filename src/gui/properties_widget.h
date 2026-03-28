#pragma once

#include "core/properties.h"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Hold_Browser.H>

#include <functional>
#include <map>
#include <string>

namespace easynmr {

class PropertiesWidget : public Fl_Group {
  public:
    PropertiesWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_properties(const MolecularProperties &props);
    void clear();
    void draw() override;

    // Fires with overlay mode string when the user clicks an overlay button.
    // mode is one of: "gradient_red", "gradient_blue", "gradient_zero",
    //                 "gradient_green", "none"
    void set_on_overlay_changed(std::function<void(std::string)> callback);
    std::string active_overlay() const { return active_overlay_; }

    // Cross-highlight: select/deselect an atom row from the structure widget.
    // atom_index is 1-based; pass -1 to deselect.
    void select_atom_row(int atom_index);

    // Fires with 1-based atom_index when the user clicks a row in the atom browser.
    void set_on_atom_selected(std::function<void(int)> callback);

  private:
    static void overlay_cb(Fl_Widget *w, void *data);
    static void atom_browser_cb(Fl_Widget *w, void *data);
    void set_active_overlay_button(Fl_Button *btn, const std::string &mode);

    Fl_Hold_Browser *main_browser_     = nullptr;
    Fl_Hold_Browser *atom_browser_     = nullptr;
    Fl_Box          *atom_col_header_  = nullptr;  // fixed column header above atom_browser_
    Fl_Box          *divider_          = nullptr;   // 1px vertical divider
    Fl_Button *overlay_fplus_btn_      = nullptr;
    Fl_Button *overlay_fminus_btn_     = nullptr;
    Fl_Button *overlay_fzero_btn_      = nullptr;
    Fl_Button *overlay_pka_btn_        = nullptr;
    Fl_Button *overlay_none_btn_       = nullptr;

    std::string active_overlay_ = "none";
    std::function<void(std::string)> on_overlay_changed_;
    std::function<void(int)> on_atom_selected_;

    // Atom index ↔ browser row mappings (rows are 1-based in Fl_Hold_Browser).
    std::map<int, int> atom_index_to_row_;
    std::map<int, int> row_to_atom_index_;
};

} // namespace easynmr
