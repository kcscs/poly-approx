#pragma once
#include <nlohmann/json.hpp>
#include <string>
#include <glm/glm.hpp>

using json = nlohmann::json;

struct CameraConfig {
  glm::uvec2 resolution;
  // Degrees
  float fovy;
};

struct ViewConfig {
  glm::vec3 eye;
  glm::vec3 at;
};

struct ExperimentConfig {
  enum class FLOAT_TYPE { FLOAT, DOUBLE };
  FLOAT_TYPE precision;
  std::string surface;
  CameraConfig camera;
  ViewConfig view;
};

void to_json(json &j, const ExperimentConfig &c);
void from_json(const json &j, ExperimentConfig &p);

void to_json(json &j, const CameraConfig &c);
void from_json(const json &j, CameraConfig &p);

void to_json(json &j, const ViewConfig &c);
void from_json(const json &j, ViewConfig &p);
