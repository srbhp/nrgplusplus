#pragma once
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include "utils/timer.hpp"
#include <algorithm> // std::min_element
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <tuple>
#include <vector>
template <typename nrgcore_type> // nrgcore_type is a type of
// We don't want to inherit the nrgcore_type.
// Just want to use few of its objects as a pointer
// We don't want to reference to any specific model here.
class fdmThermodynamics {
  nrgcore_type *nrg_object;
  int           nrgMaxIterations;
  size_t        bathDimension{0};
  std::vector<double>
      temperatureArray; // Temperature of nrg system i.e., in FDM formalism
public:
  explicit fdmThermodynamics(nrgcore_type *t_nrg_object // nrgcore_type
                             ,
                             const std::vector<double> &tArray) {
    nrg_object       = t_nrg_object;
    nrgMaxIterations = t_nrg_object->nrg_iterations_cnt;
    // Set bathDimension
    for (auto &it : t_nrg_object->bath_eigenvaluesQ) {
      bathDimension += it.size();
    }
    // Resize the results array
    temperatureArray = tArray;
  }
  void calcThermodynamics(double energyRescale) {
    std::cout << "energyRescale" << energyRescale << std::endl;
    setCurrentKeptIndex();
    std::vector<double> eigVal;
    // Sum Discarded states only
    for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
      for (size_t j = currentKeptIndex[i].size();
           j < nrg_object->eigenvaluesQ[i].size(); j++) {
        eigVal.push_back(nrg_object->eigenvaluesQ[i][j] * energyRescale);
      }
    }
    discardedEigValues.push_back(eigVal);
    eigVal.clear();
    for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
      for (size_t j = 0; j < nrg_object->eigenvaluesQ[i].size(); j++) {
        eigVal.push_back(nrg_object->eigenvaluesQ[i][j] * energyRescale);
      }
    }
    allEigValues.push_back(eigVal);
  }
  //
  void setCurrentKeptIndex() {
    if (currentKeptIndex.size() ==
        0) { // Condition only satisfy at last iteration
      for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
        currentKeptIndex.push_back({});
      }
    } else {
      currentKeptIndex = previoudKeptIndex;
    }
    previoudKeptIndex = nrg_object->eigenvaluesQ_kept_indices;
  }
  template <typename filetype> void saveFinalData(filetype *pfile) {
    // do the Calculation
    // get the Shifted the eigen energies
    std::vector<double> groundStateEnergyArr;
    for (int it = 0; it < allEigValues.size(); it++) {
      double result =
          *std::min_element(allEigValues[it].begin(), allEigValues[it].end());
      // shift the energy of the discarded states wrt the ground state
      std::transform(discardedEigValues[it].begin(),
                     discardedEigValues[it].end(),
                     discardedEigValues[it].begin(),
                     [result](double x) { return x - result; });
      // shift the energy of the discarded states wrt the ground state
      groundStateEnergyArr.push_back(result);
      // std::cout << it << "  " << allEigValues[it];
      std::cout << it << "  " << result << "  "
                << *std::max_element(allEigValues[it].begin(),
                                     allEigValues[it].end())
                << std::endl;
      // Collect true groumd state
    }
    // correct wrt last iteration
    for (size_t it = 0; it < groundStateEnergyArr.size() - 1; it++) {
      groundStateEnergyArr[it] += groundStateEnergyArr[it + 1];
    }
    //
    std::cout << "|Ground State Energy|: " << groundStateEnergyArr << std::endl;
    double gren =
        groundStateEnergyArr[0]; // Ground state energy of the last iteration
    std::transform(groundStateEnergyArr.begin(), groundStateEnergyArr.end(),
                   groundStateEnergyArr.begin(),
                   [gren](double aa) { return aa - gren; });
    std::cout << "|Ground State Energy|: " << groundStateEnergyArr << std::endl;
    // finish the final Calc
    std::vector<double> specificHeat;
    std::vector<double> entropy;
    for (int itm = 0; itm < temperatureArray.size(); itm++) {
      double kBT = temperatureArray[itm];
      // calc partition Function
      auto nrgCount = nrgMaxIterations;
      // calculate w_m first
      std::vector<double> wm;
      std::vector<double> ZmPrime;
      for (size_t id = 0; id < discardedEigValues.size(); id++) {
        double prefac = std::pow(bathDimension, nrgMaxIterations - nrgCount) *
                        regulateExp(groundStateEnergyArr[id] / kBT);
        double wi = 0;
        double zi = 0;
        for (auto &aa : discardedEigValues[id]) {
          wi += prefac * regulateExp(-aa / kBT);
          zi += regulateExp((groundStateEnergyArr[id] - aa) / kBT);
        }
        nrgCount--;
        wm.push_back(wi);
        ZmPrime.push_back(zi);
      }
      // calculate partition function
      double lpart = std::accumulate(wm.begin(), wm.end(), 0.0);
      for (auto &aa : wm) {
        aa = aa / lpart;
      }
      // Normalize w_m
      double eAv = 0;
      // Average Energy
      double eSqAv = 0;
      for (size_t id = 0; id < discardedEigValues.size(); id++) {
        for (auto &aa : discardedEigValues[id]) {
          double blm =
              regulateExp((groundStateEnergyArr[id] - aa) / kBT) / ZmPrime[id];
          eAv += (aa - groundStateEnergyArr[id]) * wm[id] * blm;
          eSqAv += (aa - groundStateEnergyArr[id]) *
                   (aa - groundStateEnergyArr[id]) * wm[id] * blm;
        }
      }
      // calc entropy
      if (itm == 0) {
        // std::cout << "DiscardedEigValues: " << discardedEigValues[0]
        //          << std::endl;
        std::cout << "W_m : " << wm << std::endl;
        std::cout << "kbt" << kBT << std::endl;
        std::cout << "eAv: " << eAv << std::endl;
        std::cout << "eSqAv: " << eSqAv << std::endl;
        std::cout << "Entropy :" << std::log(lpart) + (eAv / kBT) << std::endl;
        std::cout << "partitionFunctionArray[i]: " << lpart << std::endl;
      }
      entropy.push_back(std::log(lpart) + (eAv / kBT));
      specificHeat.push_back((eSqAv - eAv * eAv) / (kBT * kBT));
    }
    // Save the data
    std::string hstr = "";
    pfile->write(temperatureArray, hstr + "temperatureArray");
    pfile->write(entropy, hstr + "entropy");
    pfile->write(specificHeat, hstr + "specificHeat");
    // pfile->write(vecPartitions, "vecPartitions");
  }

private:
  double regulateExp(double x) {
    if (x > 100) {
      return std::exp(100);
    } else if (x < -100) {
      return std::exp(-100);
    } else {
      return std::exp(x);
    }
  }
  // Results
  //
  //
  std::vector<std::vector<size_t>> previoudKeptIndex;
  // rho.B and B.rho
  std::vector<std::vector<double>> discardedEigValues;
  std::vector<std::vector<double>> allEigValues;
  // value
public: // Give access for openchain class
  std::vector<std::vector<size_t>> currentKeptIndex;
};
