#pragma once

#include "core/properties.h"

#include <FL/Fl_Button.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Hold_Browser.H>

#include <functional>
#include <string>

namespace easynmr {

class PropertiesWidget : public Fl_Group {
  public:
    PropertiesWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_properties(const MolecularProperties &props);
    void clear();

    // Fires with overlay mode string when the user clicks an overlay button.
    // mode is one of: "gradient_red", "gradient_blue", "gradient_zero",
    //                 "gradient_green", "none"
    void set_on_overlay_changed(std::function<void(std::string)> callback);
    std::string active_overlay() const { return active_overlay_; }

  private:
    static void overlay_cb(Fl_Widget *w, void *data);
    void set_active_overlay_button(Fl_Button *btn, const std::string &mode);

    Fl_Hold_Browser *main_browser_  = nullptr;
    Fl_Hold_Browser *atom_browser_  = nullptr;
    Fl_Button *overlay_fplus_btn_   = nullptr;
    Fl_Button *overlay_fminus_btn_  = nullptr;
    Fl_Button *overlay_fzero_btn_   = nullptr;
    Fl_Button *overlay_pka_btn_     = nullptr;
    Fl_Button *overlay_none_btn_    = nullptr;

    std::string active_overlay_ = "none";
    std::function<void(std::string)> on_overlay_changed_;
};

} // namespace easynmr
