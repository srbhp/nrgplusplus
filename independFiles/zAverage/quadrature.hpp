#pragma once
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
// int itrCount = 0;
template <typename ntype>
ntype trapiziodal(std::function<ntype(ntype)> func, ntype min, ntype max) {
  // func: the function to integrate
  // max: the limit of integration
  // min: the lower limit of integration
  // calculate the integration using the trapiziodal algorithm.
  size_t itrCount = 100;
  ntype  dx       = (max - min) / itrCount;
  ntype  sum      = 0;
  for (size_t i = 0; i < itrCount; i++) {
    sum += func(min + i * dx);
  }
  return dx * sum;
}
template <typename ntype>
ntype quadrature(std::function<ntype(ntype)> func, ntype min, ntype max,
                 ntype tol = 1e-10) {
  // func: function to integrate
  // max: upper limit of integration
  // min: lower limit of integration
  // tol: tolerance
  if (max < min) {
    std::cerr << "Error: max < min" << std::endl;
    return 0;
  }
  if (max == min) {
    return 0;
  }
  // itrCount++;
  // Adaptive quadrature integration
  ntype a  = min;
  ntype b  = max;
  ntype c  = (a + b) / 2.0;
  ntype fa = func(a);
  ntype fb = func(b);
  ntype fc = func(c);
  ntype s  = (b - a) / 6.0 * (fa + 4.0 * fc + fb);
  ntype s1 = (b - a) / 12.0 *
             (fa + 4.0 * func((a + c) / 2.0) + 2.0 * fc +
              4.0 * func((c + b) / 2.0) + fb);
  if (abs(s1 - s) < tol) {
    return s1;
  }
  return quadrature(func, a, c, tol / 2) + quadrature(func, c, b, tol / 2);
}
// int main() {
//   std::function<double(double)> func = [](double x) {
//     return std::exp(-x * x);
//   };
//   double max    = 1.0;
//   double min    = 0.0;
//   double tol    = 1e-10;
//   double result = quadrature(func, max, min, tol);
//   std::cout << "Result: " << result << std::endl;
//   std::cout << "Iterations: " << itrCount << std::endl;
// }
