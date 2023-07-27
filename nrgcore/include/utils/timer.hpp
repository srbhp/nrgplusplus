#pragma once
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
class timer { // NOLINT
public:
  explicit timer(const std::string &name = std::string("main"))
      : fn_name(name) {
    std::cout << "Starting " << name << std::endl;
    reset();
  }
  void reset() { start = std::chrono::high_resolution_clock::now(); }
  ~timer() {
    std::cout << "Elapsed time for " << fn_name << " : " << getDuration()
              << "sec" << std::endl;
  }
  double getDuration() {
    // Get the duration in second
    end = std::chrono::high_resolution_clock::now();
    double time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
            .count();
    return time_span;
  }

private:
  std::string                                    fn_name;
  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point end;
};
