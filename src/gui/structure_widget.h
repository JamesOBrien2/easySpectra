#pragma once

#include <FL/Fl_Widget.H>

#include <functional>
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
};

class StructureWidget : public Fl_Widget {
  public:
    StructureWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_structure(std::vector<StructureAtom> atoms, std::vector<StructureBond> bonds);
    void set_selected_atom(int atom_index);
    void set_highlight_hydrogens(const std::vector<int> &hydrogens);
    void clear_highlight();
    void set_on_atom_selected(std::function<void(int, const std::vector<int> &)> callback);

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
};

} // namespace easynmr
