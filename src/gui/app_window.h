#pragma once

#include "core/job.h"
#include "core/pipeline.h"
#include "core/spectral_product.h"
#include "gui/colored_input_editor.h"
#include "gui/spectrum_widget.h"
#include "gui/structure_widget.h"
#include "gui/workflow_progress_widget.h"

#include <FL/Fl_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Input.H>

#include <atomic>
#include <future>
#include <map>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace easynmr {

struct QueuedJob {
    JobConfig config;
    std::string id;
    std::string status;
    std::string output_dir;
    std::string spectrum_csv;
    std::string peaks_csv;
    std::string assignments_csv;
    std::string structure_svg;
    std::string structure_atoms_csv;
    std::string structure_bonds_csv;
    std::string structure_xyz;
    std::string spectra_manifest_csv;
    std::string reaction_summary_json;
    std::string structure_product_svg;
    std::string message;
    std::string progress_stage;
    std::string progress_message;
    double progress_fraction = 0.0;
};

class AppWindow : public Fl_Double_Window {
  public:
    AppWindow(int w, int h, const char *title);
    ~AppWindow() override;
    int handle(int event) override;

  private:
    static void on_queue_job_cb(Fl_Widget *, void *userdata);
    static void on_start_queue_cb(Fl_Widget *, void *userdata);
    static void on_cancel_cb(Fl_Widget *, void *userdata);
    static void on_run_selected_cb(Fl_Widget *, void *userdata);
    static void on_preview_cb(Fl_Widget *, void *userdata);
    static void on_edit_structure_cb(Fl_Widget *, void *userdata);
    static void on_select_job_cb(Fl_Widget *, void *userdata);
    static void on_select_peak_cb(Fl_Widget *, void *userdata);
    static void on_select_atom_cb(Fl_Widget *, void *userdata);
    static void on_select_reference_cb(Fl_Widget *, void *userdata);
    static void on_select_spectrum_nucleus_cb(Fl_Widget *, void *userdata);
    static void on_select_experimental_cb(Fl_Widget *, void *userdata);
    static void on_load_experimental_cb(Fl_Widget *, void *userdata);
    static void on_clear_experimental_cb(Fl_Widget *, void *userdata);
    static void on_export_spectrum_cb(Fl_Widget *, void *userdata);
    static void on_select_example_cb(Fl_Widget *, void *userdata);
    static void on_load_example_calc_cb(Fl_Widget *, void *userdata);
    static void on_load_example_bundle_cb(Fl_Widget *, void *userdata);
    static void on_debounced_preview_cb(void *userdata);
    static void on_worker_awake(void *userdata);
    static void on_ui_tick_cb(void *userdata);

    void queue_current_input();
    void start_queue();
    void cancel_queue();
    void run_selected();
    void preview_current_input(bool show_status);
    void schedule_preview(double delay_seconds, bool show_status);
    void edit_current_structure();
    void request_preview_for_job_index(std::size_t index, bool load_visuals_if_selected);
    void launch_preview_async(const JobConfig &cfg, bool show_status);
    void poll_preview_async();
    void on_select_job();
    void on_select_peak();
    void on_select_atom();
    void on_select_reference();
    void on_select_spectrum_nucleus();
    void on_select_experimental();
    void on_load_experimental();
    void on_clear_experimental();
    void on_export_spectrum();
    void on_select_example();
    void on_load_example(bool with_experimental);
    void on_peak_picked(int group_id);
    void on_manual_shift_pair(int primary_group_id, int compare_group_id);
    void on_structure_atom_picked(int atom_index, const std::vector<int> &attached_hydrogens);
    void on_compare_structure_atom_picked(int atom_index, const std::vector<int> &attached_hydrogens);
    bool maybe_show_queue_context_menu();
    void clear_comparison_state(bool clear_selected_job);
    void apply_comparison_visuals(const QueuedJob &active_job);
    void update_structure_compare_layout(bool compare_mode);
    void highlight_comparison_for_primary_group(int primary_group_id);

    void run_worker_loop();
    void refresh_queue_browser();
    void refresh_workflow_browser(const QueuedJob *job);
    void refresh_input_syntax_mode();
    void apply_nucleus_visuals(const QueuedJob &job, const std::string &nucleus);
    void load_selected_job_visuals();
    void highlight_hydrogen_rows(const std::vector<int> &highlighted_hydrogens);
    void apply_reference_peaks(const std::string &solvent, const std::string &nucleus);
    void refresh_experimental_choice();
    void refresh_example_choice();
    void apply_active_experimental_overlay();
    void maybe_apply_example_overlay_for_active_selection();

    ColoredInputEditor *input_box_ = nullptr;
    Fl_Box *compare_input_label_ = nullptr;
    ColoredInputEditor *compare_input_box_ = nullptr;
    Fl_Input *job_name_input_ = nullptr;
    Fl_Choice *workflow_choice_ = nullptr;
    Fl_Choice *solvent_choice_ = nullptr;
    Fl_Button *preview_button_ = nullptr;
    Fl_Button *edit_button_ = nullptr;
    Fl_Button *queue_button_ = nullptr;
    Fl_Button *start_button_ = nullptr;
    Fl_Button *cancel_button_ = nullptr;
    Fl_Button *run_selected_button_ = nullptr;
    Fl_Hold_Browser *queue_browser_ = nullptr;
    WorkflowProgressWidget *workflow_progress_widget_ = nullptr;
    Fl_Box *workflow_info_bg_ = nullptr;
    Fl_Box *workflow_info_line1_ = nullptr;
    Fl_Box *workflow_info_line2_ = nullptr;
    Fl_Box *workflow_info_line3_ = nullptr;
    Fl_Hold_Browser *peak_browser_ = nullptr;
    Fl_Hold_Browser *atom_browser_ = nullptr;
    Fl_Choice *spectrum_nucleus_choice_ = nullptr;
    Fl_Choice *experimental_choice_ = nullptr;
    Fl_Choice *example_choice_ = nullptr;
    StructureWidget *structure_widget_ = nullptr;
    StructureWidget *compare_structure_widget_ = nullptr;
    Fl_Box *structure_current_label_ = nullptr;
    Fl_Box *structure_compare_label_ = nullptr;
    Fl_Button *load_example_calc_button_ = nullptr;
    Fl_Button *load_example_bundle_button_ = nullptr;
    Fl_Choice *reference_choice_ = nullptr;
    Fl_Button *load_experimental_button_ = nullptr;
    Fl_Button *clear_experimental_button_ = nullptr;
    Fl_Button *export_spectrum_button_ = nullptr;
    Fl_Box *status_box_ = nullptr;
    SpectrumWidget *spectrum_widget_ = nullptr;

    Pipeline pipeline_;
    std::vector<QueuedJob> jobs_;
    std::mutex jobs_mutex_;
    std::thread worker_thread_;
    std::atomic<bool> worker_running_ = false;
    std::atomic<bool> cancel_requested_ = false;

    enum class RunScope {
        PendingQueue,
        SelectedOnly
    };

    RunScope run_scope_ = RunScope::PendingQueue;
    std::size_t selected_run_index_ = static_cast<std::size_t>(-1);
    int active_job_index_ = -1;
    int run_total_ = 0;
    int run_done_ = 0;
    bool preview_debounce_pending_ = false;
    bool preview_debounce_show_status_ = false;
    bool preview_future_active_ = false;
    bool preview_pending_ = false;
    bool preview_pending_show_status_ = false;
    std::future<JobOutputs> preview_future_;
    JobConfig preview_inflight_cfg_;
    bool preview_inflight_show_status_ = false;
    JobConfig preview_pending_cfg_;
    std::string active_nucleus_ = "1H";
    int selected_job_index_ = -1;
    int comparison_job_index_ = -1;
    std::map<std::string, SpectralProductFiles> active_spectral_products_;
    std::map<std::string, std::vector<SpectrumPoint>> experimental_overlays_;
    std::map<std::string, std::string> experimental_overlay_paths_;
    std::map<std::string, std::string> experimental_overlay_formats_;
    std::vector<std::string> experimental_choice_keys_;
    std::string active_experimental_overlay_key_;
    int structure_area_x_ = 0;
    int structure_area_y_ = 0;
    int structure_area_w_ = 0;
    int structure_area_h_ = 0;

    std::unordered_map<int, std::vector<int>> group_to_atoms_;
    std::unordered_map<int, std::vector<int>> comparison_group_to_atoms_;
    std::vector<ComparisonPeakMarker> comparison_peak_rows_;
    std::vector<ManualShiftPair> manual_shift_pairs_;
    std::map<int, double> primary_group_to_ppm_;
    std::map<int, double> comparison_group_to_ppm_;
    std::map<int, int> group_to_peak_row_;
    std::map<int, int> peak_row_to_group_;
    std::map<int, int> atom_to_row_;
    std::map<int, int> atom_row_to_atom_;
    std::map<int, std::string> atom_row_base_label_;
    std::vector<std::string> example_case_rows_;
    std::vector<std::string> example_bundle_names_;
    std::vector<std::vector<int>> example_bundle_row_indices_;
    std::vector<int> example_choice_case_indices_;
    std::unordered_map<int, std::map<std::string, std::string>> example_job_overlay_keys_;
};

} // namespace easynmr
