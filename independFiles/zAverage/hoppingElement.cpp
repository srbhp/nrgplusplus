// This header file computes the hopping elements for z-averaging
// Based on  `Vivaldo L. Campo, Jr. and Luiz N. Oliveira
// Phys. Rev. B 72, 104432(2015)
// This file is part of the Hopping Element Library (HElib).
#include "hoppingElement.hpp"
#include <iomanip>
#include <sstream>
#include <vector>
mpfr::mpreal onsite(mpfr::mpreal lambda, const mpfr::mpreal &ZZ, int i) {
  int index = std::abs(i);
  if (index == 0)
    return (1. - pow(lambda, -ZZ)) / (ZZ * log(lambda));
  else
    return (1. - 1. / lambda) * pow(lambda, 1 - index - ZZ) / log(lambda);
}
std::vector<double> calcTn(double doubleZZ, double dlambda) {
  const int digits = 1000;
  mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
  int                       nmax   = 80; // Number of t_n requered
  int                       hcmax  = 10;
  mpfr::mpreal              lambda = dlambda;
  mpfr::mpreal              ZZ     = doubleZZ;
  std::vector<mpfr::mpreal> Hc;
  Hc.reserve(hcmax * 2 + 1);
  for (int i = 0; i <= hcmax; i++) {
    Hc.push_back(-onsite(lambda, ZZ, i));
  }
  for (int i = hcmax; i >= 0; i--) {
    Hc.push_back(onsite(lambda, ZZ, i));
  }
  std::vector<mpfr::mpreal> vec0;
  vec0.push_back(sqrt(0.5 * (1 - pow(lambda, -ZZ))));
  for (int i = 1; i <= hcmax; i++)
    vec0.push_back(sqrt(0.5 * (1. - 1. / lambda) * pow(lambda, 1. - i - ZZ)));
  for (int i = hcmax; i >= 1; i--)
    vec0.push_back(sqrt(0.5 * (1. - 1. / lambda) * pow(lambda, 1. - i - ZZ)));
  vec0.push_back(sqrt(0.5 * (1 - pow(lambda, -ZZ))));
  std::vector<mpfr::mpreal> vec1;
  std::vector<mpfr::mpreal> tnVec;
  int                       hdim = vec0.size();
  // std::cout << vec0.size() << " " << Hc.size() << std::endl;
  mpfr::mpreal aa = 0.0;
  for (int i = 0; i < hdim; i++) {
    aa += pow(vec0.at(i), 2);
    std::cout << "HA " << Hc[i] << " \t" << vec0[i] << std::endl;
  }
  std::cout << "Hdim" << vec0.size() << " " << Hc.size() << std::endl;
  /*
for(int i=0;i<hdim;i++)
vec0[i] = vec0[i] /sqrt(aa) ;
*/
  // std::cout << "REnorm " << sqrt(aa) << std::endl;
  aa = 0;
  for (int i = 0; i < hdim; i++)
    aa += Hc.at(i) * pow(vec0.at(i), 2);
  for (int i = 0; i < hdim; i++)
    vec1.push_back(vec0.at(i) * (Hc.at(i) - aa));
  aa = 0.0;
  for (int i = 0; i < hdim; i++)
    aa += pow(vec1.at(i), 2);
  for (int i = 0; i < hdim; i++)
    vec1[i] = vec1[i] / sqrt(aa);
  tnVec.push_back(sqrt(aa));
  std::vector<mpfr::mpreal> vecNM1 = vec0;
  std::vector<mpfr::mpreal> vecN   = vec1;
  for (int iN = 0; iN < nmax; iN++) {
    std::vector<mpfr::mpreal> vecNP1;
    mpfr::mpreal              a1 = 0;
    for (int i = 0; i < hdim; i++)
      a1 += Hc.at(i) * pow(vecN.at(i), 2);
    mpfr::mpreal a2 = 0;
    for (int i = 0; i < hdim; i++)
      a2 += Hc.at(i) * vecN.at(i) * vecNM1.at(i);
    for (int i = 0; i < hdim; i++)
      vecNP1.push_back(vecN.at(i) * (Hc.at(i) - a1) - a2 * vecNM1.at(i));
    aa = 0.0;
    for (int i = 0; i < hdim; i++)
      aa += pow(vecNP1.at(i), 2);
    for (int i = 0; i < hdim; i++)
      vecNP1[i] = vecNP1[i] / sqrt(aa);
    tnVec.push_back(sqrt(aa));
    vecN.swap(vecNP1);
    vecNM1.swap(vecNP1);
  }
  std::vector<double> hoppingElemnts;
  mpfr::mpreal        afactor =
      0.5 * log(lambda) * (1. + 1. / lambda) / (1. - 1. / lambda);
  for (int i = 0; i < nmax + 1; i++) {
    // Rescale to lambda^[ (n-1)/2 ]
    std::cout << "Tn" << tnVec[i] << std::endl;
    mpfr::mpreal tvar = tnVec.at(i) * pow(lambda, 0.5 * i);
    hoppingElemnts.push_back(tvar.toDouble());
    /*
     * compare with wilson's value :: from RMP paper
     if (int(doubleZZ) == 1) {
     int site = i;
     mpfr::mpreal tn = pow(lambda, -0.5 * site) * (1. / afactor) * 0.5 * (1.0
     + 1.0 / lambda) * (1.0 - pow(lambda, -site - 1)) / sqrt((1.0 - pow(lambda,
     -2 * site - 1)) * (1.0 - pow(lambda, -2 * site - 3))); std::cout << i << "
     Ex " << tn << " " << tnVec.at(i) << " Error : " << (tn - tnVec.at(i)) / tn
     << std::endl;
     }
     */
  }
  return hoppingElemnts;
}
