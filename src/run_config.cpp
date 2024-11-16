#include "run_config.hpp"
#include "experiment_config.hpp"
#include <vector>

using json = nlohmann::json;

void to_json(json &j, const RunConfig &c) {
  j = json{{"experiments", c.experiments[0]},
    {"logging",c.logsettings}};
}

void from_json(const json &j, RunConfig &p) {
  // p.experiments = j.at("experiments").template
  p.experiments = j.at("experiments").get<std::vector<ExperimentConfig>>();
  p.logsettings = j["logging"];
}
