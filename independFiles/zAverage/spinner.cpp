#include "spinner.hpp"
int main() {
  spinner sp;
  for (int i = 0; i < 500; i++) {
    sp.spin("Dots" + std::to_string(i)); //
    // sleep for 1 second
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  sp.stop();
  //
  std::cout << "Something happened" << std::endl;
  return 0;
}
