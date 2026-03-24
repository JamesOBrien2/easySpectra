#pragma once

#include "core/job.h"

#include <functional>
#include <string>

namespace easynmr {

struct PipelineProgress {
    std::string stage;
    std::string message;
    double fraction = 0.0;
};

class Pipeline {
  public:
    JobOutputs run(const JobConfig &config) const;
    JobOutputs run(const JobConfig &config, const std::function<void(const PipelineProgress &)> &progress_cb) const;
    JobOutputs preview(const JobConfig &config) const;

  private:
    JobOutputs run_impl(
        const JobConfig &config,
        const std::string &mode,
        const std::function<void(const PipelineProgress &)> *progress_cb) const;
};

} // namespace easynmr
