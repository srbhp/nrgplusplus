#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
class configPerser {
  // # start starts with comment
  // spaces are ignored
  // We need to specify the configFileName before
  // we call the constructor
  // configParser::configFileName = "config.txt" ;
  //
  // get the values of a key as
  //
  // std::cout << "key: " << configParser().get<double>("key");
  std::map<std::string, std::string> config;

public:
  const inline static std::string configFileName;
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
  template <typename T> T get(const std::string &key) {
    if (config.find(key) == config.end()) {
      throw std::runtime_error(key + ": key not found!" + "from the file " +
                               configFileName);
    }
    std::cout << "configPerser:  " << key << "  " << config[key] << std::endl;
    // double
    if constexpr (std::is_same_v<T, double>) {
      return std::stod(config[key]);
    }
    // int
    if constexpr (std::is_same_v<T, int>) {
      return std::stoi(config[key]);
    }
    // bool
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
    // string
    if constexpr (std::is_same_v<T, std::string>) {
      return config[key];
    }
    // tuple
    //
  }
  void print() {
    for (auto [key, value] : config) {
      std::cout << key << ":\t" << value << std::endl;
    }
  }
};
