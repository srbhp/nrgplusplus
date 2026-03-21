#pragma once
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Mathematical utility functions for NRG calculations.
 *
 * This header provides template utility functions including sign calculation,
 * vector operations, and string formatting commonly used in NRG algorithms.
 */

/**
 * @brief Compute the sign of a number.
 *
 * Returns:
 * - +1 if x > 0
 * - -1 if x < 0
 * - 0 if x == 0
 *
 * @tparam T Numeric type supporting relational operators
 * @param x Input value
 * @return Sign of x: -1, 0, or +1
 *
 * @note Useful for determining orientation and symmetry in quantum operators
 */
template <typename T> T sign(T x) {
  if (x > 0) {
    return 1;
  }
  if (x < 0) {
    return -1;
  }
  return 0;
}

/**
 * @brief Find index of element with minimum absolute value in vector.
 *
 * Scans vector and returns the index of element with smallest absolute value.
 * Useful for identifying smallest energy gaps or weakest couplings.
 *
 * @tparam T Element type from vector (must support std::fabs)
 * @param vec Input vector to search
 * @return Index of element with minimum |value|
 * @throw No bounds checking - behavior undefined for empty vectors
 *
 * @note Linear O(n) search; if performance critical use std::min_element with custom comparator
 */
template <typename T> size_t minIndex(const std::vector<T> &vec) {
  size_t idx   = 0;
  auto   value = std::fabs(vec[0]);
  for (size_t i = 0; i < vec.size(); i++) {
    if (std::fabs(vec[i]) < value) {
      idx   = i;
      value = std::fabs(vec[i]);
    }
  }
  return idx;
}

/**
 * @brief Convert number to string with specified decimal precision.
 *
 * Formats a number as fixed-point string with exact number of decimal places.
 * Used for consistent output formatting in data files and logging.
 *
 * @tparam T Numeric type to format
 * @param a_value Value to convert
 * @param n Number of decimal places (default: 6)
 * @return String representation with n decimal places
 *
 * @example
 * @code
 * double x = 3.14159265;
 * std::string s = to_string_with_precision(x, 3); // Returns "3.142"
 * @endcode
 */
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
