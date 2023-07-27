#include "../../include/utils/h5stream.hpp"
#include "zAverage.hpp"
#include <iomanip>
#include <vector>
// const
mpfr::mpreal dosflat(mpfr::mpreal x) {
  if (mpfr::abs(x) < 1) {
    return 0.5;
  } else {
    return 0;
  }
}
//
//
//
//
const mpfr::mpreal pi       = cos(mpfr::mpreal(-1.));
const mpfr::mpreal kfR      = 0.650; // * pi
const mpfr::mpreal normEven = 2. * (1. + sin(2. * pi * kfR) / (2. * pi * kfR));
const mpfr::mpreal normOdd  = 2. * (1. - sin(2. * pi * kfR) / (2. * pi * kfR));
//
mpfr::mpreal dosEven(mpfr::mpreal x) {
  if (mpfr::abs(x) < 1) {
    return (1.0 + cos(kfR * pi * (1. + x))) / normEven;
  } else {
    return 0;
  }
}
mpfr::mpreal dosOdd(mpfr::mpreal x) {
  if (mpfr::abs(x) < 1) {
    return (1.0 - cos(kfR * pi * (1. + x))) / normOdd;
  } else {
    return 0;
  }
}
//
//
//
int main() {
  const int digits = 1000;
  mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
  //
  // Check the normalization of the dos
  std::cout << " Norm " << trapiziodal<mpfr::mpreal>(dosEven, -1, 1).toDouble()
            << " Even " << trapiziodal<mpfr::mpreal>(dosOdd, -1, 1).toDouble()
            << " Flat " << trapiziodal<mpfr::mpreal>(dosflat, -1, 1).toDouble()
            << std::endl;
  //
  //
  //
  double             lambda = 3.00; // Lambda
  std::ostringstream ss;
  ss.str("");
  ss << "tnTwoImpurityL" << std::setprecision(3) << lambda << ".h5";
  h5stream::h5stream h5file(ss.str());
  //
  // Flat band
  // auto nz3 = zAverage(dosflat, 1.0);
  // nz3.solve(80);
  // auto s3 = nz3.getHoppingElement();
  // h5file.write(std::get<0>(s3), "tnFlat");
  // h5file.write(std::get<1>(s3), "enFlat");
  //
  //
  //
  //
  auto nz = zAverage(dosEven, 1.0);
  nz.solve(80);
  auto s1 = nz.getHoppingElement();
  h5file.write(std::get<0>(s1), "tnEven");
  h5file.write(std::get<1>(s1), "enEven");
  // //
  // auto nz2 = zAverage(dosOdd, 1.0);
  // nz2.solve(80);
  // auto s2 = nz2.getHoppingElement();
  // h5file.write(std::get<0>(s2), "tnOdd");
  // h5file.write(std::get<1>(s2), "enOdd");
  return 0;
}
