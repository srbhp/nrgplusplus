// This file computes the hopping elements for z-averaging
// Based on  `Vivaldo L. Campo, Jr. and Luiz N. Oliveira
// Phys. Rev. B 72, 104432(2015)
#include "../../include/utils/h5stream.hpp"
#include "hoppingElement.hpp"
#include <iomanip>
#include <sstream>
// Options
int main() {
  int                Nz     = 1;    // No of Z values to be calculated
  double             lambda = 3.00; // Lambda
  std::ostringstream ss;
  ss.str("");
  ss << "hoppingElementsL" << std::setprecision(3) << lambda << ".h5";
  h5stream::h5stream h5file(ss.str());
  ss.str("");
  for (int i = 1; i <= Nz; i++) {
    std::cout << "Z " << i << "/" << Nz << std::endl;
    double              zz = i / (1.0 * Nz);
    std::vector<double> tn = calcTn(zz, lambda);
    ss.str("");
    ss << "tnZ" << std::fixed << std::setprecision(3) << zz << "lmd" << lambda;
    h5file.write(tn, ss.str());
  }
  return 0;
}
