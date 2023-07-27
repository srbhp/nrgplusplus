#pragma once
#include <chrono>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
// create a spinner class
class spinner {
public:
  spinner() = default;
  // constructor
  spinner(std::vector<std::string> const &values) : values_(values) {}
  // spin the spinner
  std::string spin(const std::string &st = "") {
    // std::cout << "Starting spinner " << lstring << std::endl << std::flush;
    lstring = st;
    if (!firstrun) {
      isrunning = false;
      fut.get();
      isrunning = true;
    }
    firstrun = false;
    fut = std::async(std::launch::async, [this] { rotate(this->isrunning); });
    return lstring;
  }
  // start the spinner
  void start() { isrunning = true; }
  // stop the spinner
  void stop() {
    // if (fut.wait_for(std::chrono::milliseconds(0)) ==
    //     std::future_status::timeout)
    { isrunning = false; }
  }
  // set the spinner values
  void set_values(std::vector<std::string> const &values) { values_ = values; }
  // get the spinner values
  std::vector<std::string> get_values() const { return values_; }
  bool                     get_status() const { return isrunning; }
  void                     rotate(std::atomic_bool &rstatus) {
    // std::cout << "Starting spinner " << lstring << std::endl;
    while (rstatus) {
      auto const value = values_[counter % values_.size()];
      counter++;
      std::cout << "\x1b[2K"
                << "\r" << value << " " << this->lstring << std::flush; //
      std::this_thread::sleep_for(std::chrono::milliseconds(ifps));
    }
  }

private:
  size_t                   counter = {0};
  std::future<void>        fut;
  int                      ifps{50};
  std::atomic_bool         isrunning{true};
  bool                     firstrun{true};
  std::string              lstring;
  std::vector<std::string> values_{"⠋", "⠙", "⠹", "⠸", "⠼",
                                   "⠴", "⠦", "⠧", "⠇", "⠏"};
  // 	https://github.com/sindresorhus/cli-spinners/blob/main/spinners.json
};
