#pragma once

#include "core/spectrum.h"

#include <FL/Fl_Widget.H>

#include <functional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace easynmr {

struct PeakMarker {
    int group_id = 0;
    double center_ppm = 0.0;
    std::string label;
};

struct ReferencePeak {
    double ppm = 0.0;
    std::string label;
};

class SpectrumWidget : public Fl_Widget {
  public:
    SpectrumWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_points(std::vector<SpectrumPoint> points);
    void set_peak_markers(std::vector<PeakMarker> markers);
    void set_selected_group(int group_id);
    void set_selected_groups(const std::vector<int> &group_ids);
    void set_reference_peaks(std::vector<ReferencePeak> peaks);
    void set_highlighted_reference(int ref_index);
    void set_nucleus_label(const std::string &label);
    void set_on_peak_selected(std::function<void(int)> callback);
    bool load_from_csv(const std::string &csv_path);
    void reset_zoom();

    void draw() override;
    int handle(int event) override;

  private:
    std::pair<double, double> data_ppm_bounds() const;
    std::pair<double, double> active_ppm_bounds() const;
    double pixel_to_ppm(int pixel_x) const;
    bool point_in_plot(int px, int py) const;

    std::vector<SpectrumPoint> points_;
    std::vector<PeakMarker> peak_markers_;
    std::vector<ReferencePeak> reference_peaks_;
    int highlighted_reference_index_ = -1;
    std::string nucleus_label_ = "1H NMR Spectrum";
    std::unordered_set<int> selected_group_ids_;
    std::function<void(int)> on_peak_selected_;

    bool zoom_active_ = false;
    double view_min_ppm_ = 0.0;
    double view_max_ppm_ = 0.0;

    bool selecting_zoom_ = false;
    int drag_start_x_ = 0;
    int drag_current_x_ = 0;
};

} // namespace easynmr
