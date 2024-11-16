#include "logging.hpp"
#include <memory>
#include <stdexcept>

Logger::cat::cat(std::string name) : name(name) {}

Logger::Logger(LogSettings config) : conf(config) {
  // conf.printed_categories = config["print"];
  // conf.saved_categories = config["save"];
  // conf.save_all = config.value("save_all", true);
  // conf.print_all = config.value("print_all", true);

  print_enabled = conf.print_all;
  save_enabled = conf.save_all;

  file.open(conf.path);
}

Logger::~Logger() { file.close(); }

std::unique_ptr<Logger> Logger::global_logger = nullptr;

Logger &Logger::CreateOnce(LogSettings config) {
  if (global_logger == nullptr) {
    global_logger = std::unique_ptr<Logger>(new Logger(config));
    return *global_logger;
  } else {
    throw std::logic_error(
        "There can be only one logger and it is already initialized");
  }
}

Logger &Logger::Get() { return *global_logger; }

template <> Logger &Logger::operator<<(const cat &c) {
  if (!conf.print_all)
    print_enabled =
        conf.printed_categories.find(c.name) != conf.printed_categories.end();
  if (!conf.save_all)
    save_enabled =
        conf.saved_categories.find(c.name) != conf.saved_categories.end();
  return *this;
}

Logger::cat operator""_cat(const char *cstr, std::size_t len) {
  std::string str = std::string(cstr, len);
  return Logger::cat(str);
}

void to_json(json &j, const LogSettings &l) {
  j = json{{"print", l.printed_categories},
           {"save", l.saved_categories},
           {"print_all", l.print_all},
           {"save_all", l.save_all},
           {"path", l.path}};
}

void from_json(const json &j, LogSettings &l) {
  l.saved_categories = j["save"];
  l.printed_categories = j["print"];
  l.print_all = j.value("print_all", false);
  l.save_all = j.value("save_all", true);
  l.path = j.value("path", "render.log");
}
