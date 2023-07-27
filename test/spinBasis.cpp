#include "../include/nrgcore/qOperator.hpp"
#include "../include/utils/qmatrix.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
int main() {
  size_t spin2S  = 2;          // 2 * S
  size_t spinDim = spin2S + 1; // 2 * S + 1
                               // S_x, S_y, S_z
  qmatrix<double> Sz(spinDim, spinDim, 0);
  qmatrix<double> Splus(spinDim, spinDim, 0);
  qmatrix<double> Sminus(spinDim, spinDim, 0);
  for (size_t i = 0; i < spinDim; i++) {
    Sz(i, i) = i - spin2S / 2.0;
  }
  double ss = spin2S / 2.0;
  for (size_t i = 0; i < spinDim - 1; i++) {
    double mz        = i - spin2S / 2.0;
    Splus(i + 1, i)  = std::sqrt(ss * (ss + 1.) - mz * (mz + 1.0));
    Sminus(i, i + 1) = std::sqrt(ss * (ss + 1.) - mz * (mz + 1.0));
  }
  std::cout << "Sz = " << Sz << std::endl;
  std::cout << "Splus = " << Splus << std::endl;
  std::cout << "Sminus = " << Sminus << std::endl;
  return 0;
}
