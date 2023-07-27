#include "../../include/utils/h5stream.hpp"
#include "zAverage.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
mpfr::mpreal psr = 1.0;
mpfr::mpreal dosPseudoGap(mpfr::mpreal x) {
  if (mpfr::abs(x) < 1) {
    return pow(abs(x), psr);
  } else {
    return 0;
  }
}
//
//
//
//
//
int main() {
  const int digits = 1000;
  mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
  //
  // Check the normalization of the dos
  std::cout << " Norm "
            << trapiziodal<mpfr::mpreal>(dosPseudoGap, -1, 1).toDouble()
            << std::endl;
  //
  //
  //
  double             lambda = 3.00; // Lambda
  std::ostringstream ss;
  ss.str("");
  ss << "pseudoGap" << std::setprecision(3) << lambda << ".h5";
  h5stream::h5stream h5file(ss.str());
  //
  // Flat band
  {
    auto nz3 = zAverage(dosPseudoGap, 1.0);
    nz3.solve(80);
    auto s3 = nz3.getHoppingElement();
    h5file.write(std::get<0>(s3), "tnPseudo");
    h5file.write(std::get<1>(s3), "enPseudo");
  }
  //
  //
  //
  //
  {
    psr      = 0.50;
    auto nz3 = zAverage(dosPseudoGap, 1.0);
    nz3.solve(80);
    auto s3 = nz3.getHoppingElement();
    h5file.write(std::get<0>(s3), "tnPseudo2");
    h5file.write(std::get<1>(s3), "enPseudo2");
  }
  // //
  // auto nz2 = zAverage(dosOdd, 1.0);
  // nz2.solve(80);  // auto s2 = nz2.getHoppingElement();
  // h5file.write(std::get<0>(s2), "tnOdd");
  // h5file.write(std::get<1>(s2), "enOdd");
  //
  return 0;
}
