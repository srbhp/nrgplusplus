#pragma once
#include <chrono>
#include <iostream>
#include <string>
#include <thread>

/**
 * @brief RAII-style timer for measuring execution time of code blocks.
 *
 * Automatically measures elapsed time from construction to destruction.
 * Prints start message on creation and elapsed time on destruction.
 * Useful for profiling and logging performance-critical sections.
 *
 * @example
 * @code
 * {
 *   timer t("matrix inversion");  // Prints "Starting matrix inversion"
 *   // ... compute something ...
 * }  // Automatically prints "Elapsed time for matrix inversion: X.XXXsec"
 * @endcode
 */
class timer { // NOLINT
public:
  /**
   * @brief Constructor - starts timer with optional name.
   *
   * @param name Descriptive name for this timer instance (default: "main")
   *
   * @note Prints "Starting <name>" to stdout
   */
  explicit timer(const std::string &name = std::string("main"))
      : fn_name(name) {
    std::cout << "Starting " << name << std::endl;
    reset();
  }

  /**
   * @brief Destructor - prints elapsed time to stdout.
   *
   * Automatically called when timer goes out of scope.
   * Prints final elapsed time in seconds.
   */
  ~timer() {
    std::cout << "Elapsed time for " << fn_name << " : " << getDuration()
              << "sec" << std::endl;
  }

  /**
   * @brief Reset the start time to current moment.
   *
   * Useful for measuring intermediate segments of longer code blocks
   * or repeating measurements without creating new timer objects.
   */
  void reset() { start = std::chrono::high_resolution_clock::now(); }

  /**
   * @brief Get elapsed duration in seconds since construction/last reset().
   *
   * @return Elapsed time in seconds as double (fractional seconds supported)
   *
   * @note Does not reset the timer; use reset() to restart
   * @note Uses high_resolution_clock for maximum precision
   */
  double getDuration() {
    end = std::chrono::high_resolution_clock::now();
    double time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
            .count();
    return time_span;
  }

private:
  /// @brief Descriptive name for this timer
  std::string                                    fn_name;
  /// @brief Start time recorded at construction or last reset()
  std::chrono::high_resolution_clock::time_point start;
  /// @brief End time recorded at last getDuration() call
  std::chrono::high_resolution_clock::time_point end;
};
