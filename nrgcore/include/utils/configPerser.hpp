#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>

/**
 * @brief Small configuration-file parser for simple key/value settings.
 *
 * The parser reads a plain-text configuration file of the form:
 * key = value
 *
 * Blank lines and comment lines beginning with '#' are ignored.
 */
class configPerser {
  std::map<std::string, std::string> config;

public:
  /// @brief File name used when opening the configuration file.
  const inline static std::string configFileName;

  /**
   * @brief Read configuration entries from the configured file.
   *
   * The file is scanned line by line, trimming whitespace around keys and values.
   * Entries without a '=' are ignored.
   */
  configPerser() {
    // this function is only called once
    if (configFileName.empty()) {
      throw std::runtime_error("configParser::gfilename is empty");
    }
    std::ifstream file(configFileName);
    if (!file.is_open()) {
      std::cout << configFileName << ":File not found" << std::endl;
      return;
    }
    std::string line;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      // remove empty spaces
      // line = std::regex_replace(line, std::regex("\\s+"), "");
      std::string            key;
      std::string            value;
      std::string::size_type pos = line.find('=');
      if (pos == std::string::npos) {
        continue;
      }
      key = line.substr(0, pos);
      // remove spaces from the beginning and the end of the key
      key   = std::regex_replace(key, std::regex("\\s+"), "");
      value = line.substr(pos + 1);
      // remove spaces from the beginning and the end of the key
      value       = std::regex_replace(value, std::regex("^ +| +$"), "$1");
      config[key] = value;
    }
  }
  /**
   * @brief Retrieve a configuration value converted to the requested type.
   *
   * @tparam T Output type: double, int, bool, or std::string.
   * @param key Name of the configuration entry.
   * @return Parsed value of that key.
   * @throws std::runtime_error If the key is missing or has an invalid format.
   */
  template <typename T> T get(const std::string &key) {
    if (config.find(key) == config.end()) {
      throw std::runtime_error(key + ": key not found!" + "from the file " +
                               configFileName);
    }
    std::cout << "configPerser:  " << key << "  " << config[key] << std::endl;
    if constexpr (std::is_same_v<T, double>) {
      return std::stod(config[key]);
    }
    if constexpr (std::is_same_v<T, int>) {
      return std::stoi(config[key]);
    }
    if constexpr (std::is_same_v<T, bool>) {
      bool check = false;
      if (config[key] == "true") {
        check = true;
        return true;
      }
      if (config[key] == "false") {
        check = true;
        return false;
      }
      if (!check) {
        throw std::runtime_error(key + ": key not found!" + "from the file " +
                                 configFileName);
      }
    }
    if constexpr (std::is_same_v<T, std::string>) {
      return config[key];
    }
  }

  /**
   * @brief Print all parsed configuration entries to standard output.
   */
  void print() {
    for (auto [key, value] : config) {
      std::cout << key << ":\t" << value << std::endl;
    }
  }
};
