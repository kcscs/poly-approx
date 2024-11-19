#include "experiment_config.hpp"
#include "glm_parsers.hpp"
#include <glm/gtx/transform2.hpp>

void to_json(json &j, const ExperimentConfig &c) {
  j = json{{"precision", c.precision == ExperimentConfig::FLOAT_TYPE::FLOAT
                             ? "float"
                             : "double"},
           {"surface", c.surface},
           {"camera", c.camera},
           {"view", c.view},
           {"trace_settings", c.trace_settings},
           {"title", c.title},
           {"enabled", c.enabled},
           {"light_dir", c.light_dir}};
}

void from_json(const json &j, ExperimentConfig &p) {
  std::string prec_str = j["precision"];
  p.surface = j["surface"];
  p.precision = prec_str == "float" ? ExperimentConfig::FLOAT_TYPE::FLOAT
                                    : ExperimentConfig::FLOAT_TYPE::DOUBLE;
  p.camera = j["camera"];
  p.view = j["view"];
  p.trace_settings = j["trace_settings"];
  p.title = j["title"];
  p.enabled = j.value("enabled", true);

  glm::vec4 default_light_dir = glm::vec4(p.view.at - p.view.eye, 1);
  default_light_dir =
      glm::rotate(glm::radians(40.0f), glm::vec3(0.0f, 1.0f, 0.0f)) *
      default_light_dir;
  default_light_dir =
      glm::rotate(glm::radians(-30.0f), glm::vec3(1.0f, 0.0f, 0.0f)) *
      default_light_dir;

  p.light_dir = glm::normalize(j.value("light_dir", default_light_dir));
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
