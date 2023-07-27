#pragma once
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
// Sign function declarations
template <typename T> T sign(T x) {
  if (x > 0) {
    return 1;
  }
  if (x < 0) {
    return -1;
  }
  return 0;
}
template <typename T> size_t minIndex(const std::vector<T> &vec) {
  size_t idx   = 0;
  auto   value = std::fabs(vec[0]); // make it auto
  for (size_t i = 0; i < vec.size(); i++) {
    if (std::fabs(vec[i]) < value) {
      idx   = i;
      value = std::fabs(vec[i]);
    }
  }
  return idx;
}
//
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
//
