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
    std::string multiplicity = "singlet";
    double j_hz = 0.0;
};

struct ComparisonPeakMarker {
    int group_id = 0;
    double center_ppm = 0.0;
    std::string label;
};

struct ManualShiftPair {
    int primary_group_id = 0;
    int compare_group_id = 0;
    double primary_ppm = 0.0;
    double compare_ppm = 0.0;
    double delta_ppm = 0.0;
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
    void set_selected_comparison_groups(const std::vector<int> &group_ids);
    void set_reference_peaks(std::vector<ReferencePeak> peaks);
    void set_highlighted_reference(int ref_index);
    void set_nucleus_label(const std::string &label);
    void set_render_settings(const std::string &line_shape, double fwhm_hz, double frequency_mhz);
    void set_experimental_points(std::vector<SpectrumPoint> points);
    void set_comparison_points(std::vector<SpectrumPoint> points);
    void clear_comparison_points();
    void set_comparison_peak_markers(std::vector<ComparisonPeakMarker> markers);
    void set_manual_shift_pairs(std::vector<ManualShiftPair> pairs);
    void clear_manual_shift_pairs();
    void clear_experimental_points();
    void set_on_peak_selected(std::function<void(int)> callback);
    void set_on_manual_shift_pair(std::function<void(int, int)> callback);
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
    std::vector<SpectrumPoint> experimental_points_;
    std::vector<SpectrumPoint> comparison_points_;
    std::vector<PeakMarker> peak_markers_;
    std::vector<ComparisonPeakMarker> comparison_peak_markers_;
    std::vector<ManualShiftPair> manual_shift_pairs_;
    std::vector<ReferencePeak> reference_peaks_;
    int highlighted_reference_index_ = -1;
    std::string nucleus_label_ = "1H NMR Spectrum";
    std::unordered_set<int> selected_group_ids_;
    std::unordered_set<int> selected_comparison_group_ids_;
    std::function<void(int)> on_peak_selected_;
    std::function<void(int, int)> on_manual_shift_pair_;
    std::string line_shape_ = "lorentzian";
    double fwhm_hz_ = 1.0;
    double frequency_mhz_ = 400.0;
    int armed_primary_group_id_ = 0;

    bool zoom_active_ = false;
    double view_min_ppm_ = 0.0;
    double view_max_ppm_ = 0.0;

    bool selecting_zoom_ = false;
    int drag_start_x_ = 0;
    int drag_current_x_ = 0;
};

} // namespace easynmr
