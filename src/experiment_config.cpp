#include "experiment_config.hpp"
#include "glm_parsers.hpp"

void to_json(json &j, const ExperimentConfig &c) {
  j = json{{"precision", c.precision == ExperimentConfig::FLOAT_TYPE::FLOAT
                             ? "float"
                             : "double"},
           {"surface", c.surface},
           {"camera", c.camera},
           {"view", c.view},
           {"trace_settings", c.trace_settings}};
}

void from_json(const json &j, ExperimentConfig &p) {
  std::string prec_str = j["precision"];
  p.surface = j["surface"];
  p.precision = prec_str == "float" ? ExperimentConfig::FLOAT_TYPE::FLOAT
                                    : ExperimentConfig::FLOAT_TYPE::DOUBLE;
  p.camera = j["camera"];
  p.view = j["view"];
  p.trace_settings = j["trace_settings"];
}

void to_json(json &j, const CameraConfig &c) {
  j = json{{"resolution", c.resolution}, {"fovy", c.fovy}};
}

void from_json(const json &j, CameraConfig &c) {
  c.resolution = j["resolution"];
  c.fovy = j["fovy"];
}

void to_json(json &j, const ViewConfig &v) {
  j = json{{"eye", v.eye}, {"target", v.at}};
}

void from_json(const json &j, ViewConfig &v) {
  v.eye = j["eye"];
  v.at = j["target"];
}
