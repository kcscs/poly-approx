#pragma once
#include <vector>

#include <nlohmann/json.hpp>

#include "experiment_config.hpp"

struct RunConfig {
  std::vector<ExperimentConfig> experiments;
};

void to_json(json &j, const RunConfig &c);
void from_json(const json &j, RunConfig &p);
