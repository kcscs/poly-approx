#include "chebseg.hpp"
#include "experiment_config.hpp"
#include "logging.hpp"
#include "monseg.hpp"
#include "renderer.hpp"
#include "run_config.hpp"
#include "scenelib.hpp"
#include "tracer.hpp"
#include <cassert>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>

#include <argparse/argparse.hpp>
#include <glm/glm.hpp>
#include <memory>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include <omp.h>

void run(const RunConfig &config);

int main(int argc, char *argv[]) {

  argparse::ArgumentParser program("Polynomial tracer");
  program.add_argument("config")
      .help("Path to configuration file describing tests")
      .required();

  program.add_argument("-o", "--outdir")
      .help("Output directory")
      .default_value("experiments");

  // try {
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

    Logger &logger = Logger::CreateOrRecreate(rc.logsettings);
    logger << Logger::cat("init");

    logger << "Parsed config:\n" << json{rc}.dump(4) << "\n";

    //////////////////////
    run(rc);
    //////////////////////
  } else {
    std::cerr << config_path << " does not exist\n";
  }
  // } catch (const std::exception &err) {
  //   std::cerr << err.what() << std::endl;
  //   std::cerr << program;
  //   return 1;
  // }
  return 0;
}

template <typename FT>
json cheb_exp(const ExperimentConfig &exp, std::string exp_dir) {
  Logger& log = Logger::Get();

  std::unique_ptr<PolynomialTracer<FT>> tracer;
  std::string trace_method = exp.trace_settings.value("algorithm", "global_chebyshev");
  log << "init"_cat << "Trace method: "<<trace_method<<"\n";
  if(trace_method == "global_chebyshev")
    tracer = std::make_unique<PolynomialTracer<FT>>(exp.trace_settings);
  else if(trace_method == "first_root_chebyshev")
    tracer = std::make_unique<FirstRootPolynomialTracer<FT>>(exp.trace_settings);
  else {
    log << "unknown trace algorithm\n";
    return json("ERROR - unknown trace algorithm");
  }

  Renderer<FT> ren;
  typename Types<FT>::SurfaceFunction surf =
      Scenes<FT>::GetSceneByName(exp.surface);

  // old light_dir: 1, -3, 1
  std::string image_name = "render.png";
  json render_data = ren.render(surf, exp.camera, exp.view,
                                exp.light_dir,
                                exp_dir + "/" + image_name, tracer.get());
  render_data["image"] = image_name;
  return render_data;
}

// Run experiments
void run(const RunConfig &config) {
  Logger &log = Logger::Get();

  auto now = std::chrono::system_clock::now();
  std::string timestamp_str = std::format("{:%Y-%m-%d_%H:%M:%S}", now);
  std::string workdir = config.output_dir + "/" + timestamp_str;
  std::filesystem::create_directories(workdir);

  std::string logoutput = config.logsettings.path;

  for (const ExperimentConfig &exp : config.experiments) {
    if (!exp.enabled) {
      log << "init"_cat << "skipping " << exp.title << "\n";
      continue;
    }

    std::filesystem::create_directory(workdir + "/" + exp.title);
    LogSettings logset = config.logsettings;
    std::string exp_dir = workdir + "/" + exp.title;
    logset.path = workdir + "/" + exp.title + "/" + logoutput;
    Logger &log = Logger::CreateOrRecreate(logset);

    json metadata;
    if (exp.precision == ExperimentConfig::FLOAT_TYPE::FLOAT) {
      log << "init"_cat << "\nUsing FLOAT precision\n";
      metadata = cheb_exp<float>(exp, exp_dir);
    } else if (exp.precision == ExperimentConfig::FLOAT_TYPE::DOUBLE) {
      log << "init"_cat << "\nUsing DOUBLE precision\n";
      metadata = cheb_exp<double>(exp, exp_dir);
    } else
      assert(false);

    std::ofstream config_copy(exp_dir + "/experiment_config.json");
    config_copy << json(exp).dump(4);
    config_copy.close();
    std::ofstream metadata_file(exp_dir + "/metadata.json");
    metadata_file << metadata.dump(4);
    metadata_file.close();
  }
}
