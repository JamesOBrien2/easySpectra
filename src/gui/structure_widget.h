#pragma once

#include <FL/Fl_Widget.H>

#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace easynmr {

struct StructureAtom {
    int atom_index = 0; // 1-based atom index in full molecule
    std::string element;
    double x = 0.0;
    double y = 0.0;
    std::vector<int> attached_hydrogens; // 1-based indices
};

struct StructureBond {
    int atom_a = 0; // 1-based atom index
    int atom_b = 0; // 1-based atom index
    int order = 1;
    std::string stereo_style = "none"; // none | wedge | dash
    int stereo_from_atom = 0;          // 1-based atom index at wedge apex (if any)
};

class StructureWidget : public Fl_Widget {
  public:
    StructureWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_structure(std::vector<StructureAtom> atoms, std::vector<StructureBond> bonds);
    void set_empty_message(std::string message);
    void set_selected_atom(int atom_index);
    void set_highlight_hydrogens(const std::vector<int> &hydrogens);
    void set_highlight_palette(Fl_Color selected_fill, Fl_Color selected_border, Fl_Color attached_fill);
    void clear_highlight();
    void set_on_atom_selected(std::function<void(int, const std::vector<int> &)> callback);

    // Atom overlay for property visualisation.
    // values: map from 1-based atom_index → normalised [0,1] value.
    // mode: "gradient_red" (f+), "gradient_blue" (f-), "gradient_zero" (f0),
    //       "gradient_green" (pKa), "none".
    void set_atom_overlay(std::map<int, double> values, std::string mode);
    void clear_atom_overlay();

    // Properties panel cross-highlight — amber ring drawn on top of the gradient overlay,
    // independent of the NMR selected_atom_index_ / filled-disc highlight.
    void set_props_selected_atom(int atom_index);  // -1 to clear
    void clear_props_selected_atom();

    void draw() override;
    int handle(int event) override;

  private:
    std::pair<int, int> to_screen(double px, double py, double min_x, double max_x, double min_y, double max_y) const;

    std::vector<StructureAtom> atoms_;
    std::vector<StructureBond> bonds_;
    int selected_atom_index_ = -1;
    std::unordered_set<int> highlighted_hydrogens_;
    std::function<void(int, const std::vector<int> &)> on_atom_selected_;
    std::unordered_map<int, std::pair<int, int>> screen_positions_;
    std::string empty_message_ = "No 2D structure";
    Fl_Color selected_fill_color_ = FL_WHITE;
    Fl_Color selected_border_color_ = FL_DARK_GREEN;
    Fl_Color attached_fill_color_ = FL_LIGHT2;
    // Atom overlay state
    std::map<int, double> atom_overlay_values_;  // 1-based atom index → normalised value [0,1]
    std::string atom_overlay_mode_ = "none";
    int props_selected_atom_ = -1;               // properties panel selection (amber ring)
};

} // namespace easynmr
