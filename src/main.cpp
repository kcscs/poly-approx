#include "chebseg.hpp"
#include "experiment_config.hpp"
#include "fun.hpp"
#include "monseg.hpp"
#include "renderer.hpp"
#include "run_config.hpp"
#include "scenelib.hpp"
#include "seg.hpp"
#include "tracer.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

#include <argparse/argparse.hpp>
#include <glm/glm.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include <omp.h>

void run(const RunConfig &config);

int main(int argc, char *argv[]) {

  argparse::ArgumentParser program("Polynomial tracer");
  program.add_argument("config")
      .help("Path to configuration file describing tests")
      .required();
  ;
  program.add_argument("-o", "--outdir")
      .help("Output directory")
      .default_value("experiments");

  try {
    program.parse_args(argc, argv);

    std::string config_path = program.get("config");
    std::string output_directory = program.get("-o");
    std::cout << "Config: " << config_path << "\n";
    std::cout << "Output directory: " << config_path << "\n";

    // Load configuration
    if (std::filesystem::exists(config_path)) {
      std::ifstream config_file(config_path);
      json config;
      config_file >> config;
      config_file.close();

      RunConfig rc = config.get<RunConfig>();
      std::cout << "Parsed config:\n" << json{rc}.dump(4) << "\n";

      //////////////////////
      run(rc);
      //////////////////////
    } else {
      std::cerr << config_path << "does not exist\n";
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }
  return 0;
}

// Run experiments
void run(const RunConfig &config) {
  for (const ExperimentConfig &exp : config.experiments) {
    MonSeg<float> segf({}, 3,5);
    MonSeg<double> segd({},3,5);

    ChebSeg<float> csf({},3,5);
    ChebSeg<double> csd({},3,5);

    MonSeg<float> conv = MonSeg<float>::FitAtChebPoints(csf);

    PolynomialTracer<float> tracer;
    Renderer<float> ren;
    ren.render(Scenes<float>::GetSceneByName("torus"), exp.camera,
                        exp.view, glm::normalize(glm::vec3(1.0f, -3.0f, 1.0f)),
                        "Test.png", &tracer);
  }
}
