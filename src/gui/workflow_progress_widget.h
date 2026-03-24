#pragma once

#include <FL/Fl_Widget.H>

#include <string>
#include <vector>

namespace easynmr {

class WorkflowProgressWidget : public Fl_Widget {
  public:
    WorkflowProgressWidget(int x, int y, int w, int h, const char *label = nullptr);

    void set_steps(std::vector<std::string> steps);
    void set_progress_state(
        int active_step_index,
        int failed_step_index,
        bool running,
        bool done,
        bool failed,
        double fraction,
        std::string status_text);

    void draw() override;

  private:
    std::vector<std::string> steps_;
    int active_step_index_ = -1;
    int failed_step_index_ = -1;
    bool running_ = false;
    bool done_ = false;
    bool failed_ = false;
    double fraction_ = 0.0;
    std::string status_text_;
};

} // namespace easynmr

