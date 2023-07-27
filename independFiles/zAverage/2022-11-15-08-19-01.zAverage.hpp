#pragma once
#include "mpreal.h" // This file is copied from https://github.com/advanpix/mpreal
#include "quadrature.hpp"
#include "spinner.hpp"
#include <iostream>
#include <string>
#include <vector>
/**
 * @brief [TODO:summary]
 */
class zAverage {
public:
  explicit zAverage(const std::function<mpfr::mpreal(mpfr::mpreal)> ndos,
                    double                                          nzz)
      :                //
        bandwidth(1.), //
        lambda(3.0),   //
        zValue(nzz),   //
        dos(ndos),     //
        dosByX([=](mpfr::mpreal x) { return this->dos(x) / x; })
  //
  {
    // set the precision
    const int digits = 1000;
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
    //
    //
    setlogmess();
    setFirstVector();
  }
  std::tuple<std::vector<double>, std::vector<double>> getHoppingElement() {
    // get rescaled coefficients tn
    std::vector<double> nvec(hoppingElement.size());
    for (size_t i = 0; i < nvec.size(); ++i) {
      nvec[i] = (hoppingElement[i] * pow(lambda, 0.5 * i)).toDouble();
    }
    std::vector<double> ovec(wilsonSiteEnergy.size());
    for (size_t i = 0; i < ovec.size(); ++i) {
      nvec[i] = (wilsonSiteEnergy[i] * pow(lambda, 0.5 * i)).toDouble();
    }
    return {nvec, ovec};
  };
  void solve(size_t nitr = 40) {
    std::cout << "Starting solve..." << std::endl;
    // spinner sp;
    for (size_t i = 0; i < nitr; i++) {
      std::cout << "+";
      // sp.spin(std::to_string(i));
      setNextVector();
    }
    // sp.stop();
    std::cout << std::endl;
  }
  void setZvalue(double zvalue) { zValue = zvalue; }
  void setBandwidth(double bw) { bandwidth = bw; }
  // Private class members
private:
  void setNextVector() {
    // <0|H_C |0>
    mpfr::mpreal aa  = 0;
    mpfr::mpreal aa2 = 0;
    for (size_t i = 0; i < ndim; i++) {
      aa += Hamiltonian[i] * nextVector[i] * nextVector[i];
      aa2 += Hamiltonian[i] * nextVector[i] * preVector[i];
    }
    wilsonSiteEnergy.push_back(aa);
    std::vector<mpfr::mpreal> tmpVector(ndim);
    // set the first vector
    for (size_t i = 0; i < ndim; i++) {
      tmpVector[i] = Hamiltonian[i] * nextVector[i] - aa * nextVector[i] -
                     aa2 * preVector[i];
    }
    normalize(tmpVector);
    preVector  = nextVector;
    nextVector = tmpVector;
  }
  void normalize(std::vector<mpfr::mpreal> &nvec) {
    mpfr::mpreal norm = 0;
    for (size_t i = 0; i < ndim; i++) {
      norm += nvec[i] * nvec[i];
    }
    norm = sqrt(norm);
    hoppingElement.push_back(mpfr::mpreal(norm));
    for (size_t i = 0; i < ndim; i++) {
      nvec[i] /= norm;
    }
  }
  void setFirstVector() {
    // std::cout << "Done " << energymess.size() << std::endl;
    ndim = energymess.size() - 2;
    std::vector<mpfr::mpreal> vec0(ndim);
    nextVector.resize(ndim);
    Hamiltonian.resize(ndim);
    std::cout << "Ndim " << ndim << std::endl;
    // f_0 operator
    // see the Oliveira paper for details
    // spinner     sp;
#pragma omp parallel for
    for (int i = -energyPoints + 1; i < -1; i++) {
      // std::cout << "Point\t" << i << std::endl;
      // sp.spin("setFirstVector" + std::to_string(i));
      mpfr::mpreal max         = energymess.at(i + 1);
      mpfr::mpreal min         = energymess.at(i);
      auto         numerator   = trapiziodal(dos, min, max);
      auto         denominator = trapiziodal(dosByX, min, max);
      // std::cout << "max" << max << " min" << min.toDouble() << std::endl;
      // std::cout << "numerator" << numerator << " denominator" << denominator
      std::size_t itr = i - (-energyPoints + 1);
      //           << std::endl;
      Hamiltonian[itr] = numerator / denominator;
      // vec0[itr]        = sqrt(abs(max - min) / (2. * bandwidth));
      vec0[itr] = sqrt(abs(numerator));
      itr++;
    }
#pragma omp parallel for
    for (size_t i = 1; i < energyPoints - 1; i++) {
      std::size_t itr = energyPoints - 3 + i;
      // sp.spin("setFirstVector" + std::to_string(i));
      mpfr::mpreal max         = energymess.at(i);
      mpfr::mpreal min         = energymess.at(i + 1);
      auto         numerator   = trapiziodal(dos, min, max);
      auto         denominator = trapiziodal(dosByX, min, max);
      Hamiltonian[itr]         = numerator / denominator;
      vec0[itr]                = sqrt(abs(numerator));
      // vec0[itr]                = sqrt(abs(max - min) / (2. * bandwidth));
      itr++;
    }
    // check the normalization of vec0
    for (size_t i = 0; i < ndim; i++) {
      std::cout << "Ha " << Hamiltonian[i] << "\t " << vec0[i] << std::endl;
    }
    // <0|H_C |0>
    mpfr::mpreal aa = 0;
    for (size_t i = 0; i < ndim; i++) {
      aa += Hamiltonian[i] * vec0[i] * vec0[i];
    }
    wilsonSiteEnergy.push_back(aa);
    // set the first vector
    for (size_t i = 0; i < ndim; i++) {
      nextVector[i] = Hamiltonian[i] * vec0[i] - aa * vec0[i];
    }
    normalize(nextVector);
    preVector = vec0;
    // for (auto const &ta : nextVector) {
    //   std::cout << "itr " << ta << std::endl;
    // }
    //
  }
  void setEnergyMess(const std::map<int, mpfr::mpreal> &enmess, size_t npoint) {
    energyPoints = npoint;
    energymess   = enmess;
  }
  void setlogmess() {
    // set the log mess  of the bath band
    energymess[1]  = bandwidth;
    energymess[-1] = -bandwidth;
    for (size_t i = 2; i < energyPoints; i++) {
      energymess[i]  = bandwidth * pow(lambda, 2.0 - 1. * i - zValue);
      energymess[-i] = -bandwidth * pow(lambda, 2.0 - 1. * i - zValue);
    }
    // for (auto const &aa : energymess) {
    //   std::cout << aa.first << "\t" << aa.second << std::endl;
    // }
  }
  // number of the energy points in the mess
  // and dimension of the Hamiltonian
  size_t                                    energyPoints = 10;
  size_t                                    ndim         = energyPoints * 2 - 2;
  mpfr::mpreal                              bandwidth    = 1.;
  mpfr::mpreal                              lambda       = 3.0;
  mpfr::mpreal                              zValue       = 1.0;
  std::vector<mpfr::mpreal>                 Hamiltonian;
  std::vector<mpfr::mpreal>                 nextVector;
  std::vector<mpfr::mpreal>                 preVector;
  std::vector<mpfr::mpreal>                 hoppingElement;
  std::vector<mpfr::mpreal>                 wilsonSiteEnergy;
  std::map<int, mpfr::mpreal>               energymess;
  std::function<mpfr::mpreal(mpfr::mpreal)> dos;
  std::function<mpfr::mpreal(mpfr::mpreal)> dosByX;
};
