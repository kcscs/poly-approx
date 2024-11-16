#pragma once
#include <vector>

#include <nlohmann/json.hpp>

#include "experiment_config.hpp"
#include "logging.hpp"

struct RunConfig {
  LogSettings logsettings;
  std::vector<ExperimentConfig> experiments;
};

void to_json(json &j, const RunConfig &c);
void from_json(const json &j, RunConfig &p);
