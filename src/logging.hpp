#pragma once

#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_set>

using json = nlohmann::json;

struct LogSettings {
  std::unordered_set<std::string> printed_categories;
  std::unordered_set<std::string> saved_categories;
  bool print_all = false;
  bool save_all = true;
  std::string path;
};

class Logger : public std::enable_shared_from_this<Logger> {
public:
  struct cat {
    std::string name;
    cat(std::string name);
  };

  static Logger &CreateOnce(LogSettings config);
  static Logger &Get();

  Logger(const Logger &other) = delete;
  Logger(const Logger &&other) = delete;
  ~Logger();

  template <typename T> Logger &operator<<(const T &val) {
    if (print_enabled)
      std::cout << val;
    if (save_enabled)
      file << val;
    return *this;
  }

  template <typename T> Logger &operator<<(const std::vector<T> &val) {
    *this << "vec[";
    for (const auto &v : val)
      *this << v << ",";
    *this << "]";
    return *this;
  }

  template <> Logger &operator<<(const cat &c);

private:
  LogSettings conf;
  bool print_enabled = true;
  bool save_enabled = true;
  std::ofstream file;

  static std::unique_ptr<Logger> global_logger;
  Logger(LogSettings config);
};

Logger::cat operator""_cat(const char *cstr, std::size_t len);

void to_json(json &j, const LogSettings &l);
void from_json(const json &j, LogSettings &l);
