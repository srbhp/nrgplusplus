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
#include <utility>
#include <vector>
//
template <typename nrgcore_type> //
class fdmThermodynamics {
  nrgcore_type                                      *nrg_object;
  int                                                nrgMaxIterations;
  size_t                                             bathDimension{0};
  std::vector<qOperator>                            *qsysOPerator{nullptr};
  std::vector<std::vector<std::map<size_t, double>>> qsysOPeratorValue;
  std::vector<double>                                temperatureArray;
  // Temperature of nrg system i.e., in FDM formalism
public:
  explicit fdmThermodynamics(nrgcore_type *t_nrg_object // nrgcore_type
                             ,
                             std::vector<double> tArray)
      : nrg_object(t_nrg_object),
        nrgMaxIterations(t_nrg_object->nrg_iterations_cnt),
        temperatureArray(std::move(tArray)) {
    // Set bathDimension
    for (auto &it : t_nrg_object->bath_eigenvaluesQ) {
      bathDimension += it.size();
    }
    std::cout << "bathDimension: " << bathDimension << std::endl;
    // Resize the results array
  }
  void setSystemOperator(std::vector<qOperator> *qsysOPerator_) {
    qsysOPerator = qsysOPerator_;
  }
  void calcThermodynamics(double energyRescale) { // NOLINT
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
    std::cout << "Number of eigen Values: " << eigVal.size() << std::endl;
    discardedEigValues.push_back(eigVal);
    eigVal.clear();
    // add  all states
    for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
      for (size_t j = 0; j < nrg_object->eigenvaluesQ[i].size(); j++) {
        eigVal.push_back(nrg_object->eigenvaluesQ[i][j] * energyRescale);
      }
    }
    allEigValues.push_back(eigVal);
    // Save the expectation value of the system operator
    if (qsysOPerator != nullptr) {
      std::vector<std::map<size_t, double>> qvalue(qsysOPerator->size());
      for (size_t iq = 0; iq < qsysOPerator->size(); iq++) {
        size_t icounter{0};
        for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
          auto qopr = qsysOPerator->at(iq).get(i, i);
          if (qopr) {
            auto *qmat = qopr.value();
            for (size_t j = currentKeptIndex[i].size();
                 // Only discarded states are saved
                 j < nrg_object->eigenvaluesQ[i].size(); j++) {
              if (icounter + j - currentKeptIndex[i].size() >=
                  discardedEigValues.back().size()) {
                std::cout << "icounter: " << icounter << std::endl;
              }
              qvalue[iq][icounter + j - currentKeptIndex[i].size()] =
                  qmat->at(j, j);
            }
          }
          icounter +=
              (nrg_object->eigenvaluesQ[i].size() - currentKeptIndex[i].size());
        }
      }
      qsysOPeratorValue.push_back(qvalue);
    }
  }
  //
  void setCurrentKeptIndex() {
    if (currentKeptIndex.empty()) { // Condition only satisfy at last iteration
      for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
        currentKeptIndex.emplace_back();
      }
    } else {
      currentKeptIndex = previoudKeptIndex;
    }
    previoudKeptIndex = nrg_object->eigenvaluesQ_kept_indices;
  }
  template <typename filetype> void saveFinalData(filetype *pfile) { // NOLINT
    // do the Calculation
    // get the Shifted the eigen energies
    // for (size_t iq = 0; iq < qsysOPerator->size(); iq++) {
    //   for (size_t id = 0; id < discardedEigValues.size(); id++) {
    //     std::cout << "Itr: ";
    //     for (auto const &[ie, qval] : qsysOPeratorValue[id][iq]) {
    //       std::cout << qval << " ";
    //     }
    //     std::cout << std::endl;
    //   }
    // }
    //
    std::vector<double> groundStateEnergyArr(allEigValues.size(), 0);
    for (int it = 0; it < allEigValues.size(); it++) {
      double result =
          *std::min_element(allEigValues[it].begin(), allEigValues[it].end());
      // shift the energy of the discarded states wrt the ground state
      std::transform(discardedEigValues[it].begin(),
                     discardedEigValues[it].end(),
                     discardedEigValues[it].begin(),
                     [result](double x) { return x - result; });
      // shift the energy of the discarded states wrt the ground state
      groundStateEnergyArr[it] = result;
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
    double gren = groundStateEnergyArr[0]; // Ground state energy of the
                                           // last iteration
    std::transform(groundStateEnergyArr.begin(), groundStateEnergyArr.end(),
                   groundStateEnergyArr.begin(),
                   [gren](double aa) { return aa - gren; });
    std::cout << "|Ground State Energy|: " << groundStateEnergyArr << std::endl;
    // finish the final Calc
    std::vector<double>              specificHeat(temperatureArray.size(), 0.0);
    std::vector<double>              entropy(temperatureArray.size(), 0.0);
    std::vector<std::vector<double>> qsysOPeratorExpValue;
    if (qsysOPerator != nullptr) {
      qsysOPeratorExpValue.resize(
          qsysOPerator->size(),
          std::vector<double>(temperatureArray.size(), 0.0));
    }
    std::cout << "Done " << std::endl;
    for (int itm = 0; itm < temperatureArray.size(); itm++) {
      double kBT = temperatureArray[itm];
      // calc partition Function
      auto nrgCount = nrgMaxIterations;
      // calculate w_m first
      std::vector<double> wm(discardedEigValues.size(), 0);
      std::vector<double> ZmPrime(discardedEigValues.size(), 0);
      for (size_t id = 0; id < discardedEigValues.size(); id++) {
        double prefac = std::pow(bathDimension, nrgMaxIterations - nrgCount);
        double zi     = 0;
        for (auto &aa : discardedEigValues[id]) {
          zi += regulateExp((groundStateEnergyArr[id] - aa) / kBT);
        }
        nrgCount--;
        wm[id]      = zi * prefac;
        ZmPrime[id] = zi;
      }
      // calculate partition function
      double lpart = std::accumulate(wm.begin(), wm.end(), 0.0);
      for (auto &aa : wm) {
        aa = aa / lpart;
      }
      std::cout << "Done 2" << std::endl;
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
      // 			// calc qsysOPeratorExpValue
      if (qsysOPerator != nullptr) {
        for (size_t iq = 0; iq < qsysOPerator->size(); iq++) {
          double tvalue = 0;
          for (size_t id = 0; id < discardedEigValues.size(); id++) {
            std::cout << "Itr: " << iq << " "
                      << qsysOPeratorValue[id][iq].size() << std::endl;
            for (auto const &[ie, qval] : qsysOPeratorValue[id][iq]) {
              if (ie >= discardedEigValues[id].size()) {
                std::cout << "Error: " << ie << " "
                          << discardedEigValues[id].size() << std::endl;
              }
              double blm = regulateExp((groundStateEnergyArr[id] -
                                        discardedEigValues[id][ie]) /
                                       kBT) /
                           ZmPrime[id];
              tvalue += qval * wm[id] * blm;
            }
          }
          qsysOPeratorExpValue[iq][itm] = tvalue;
        }
      }
      std::cout << "Done 3" << std::endl;
      // calc entropy
      if (itm == 0) {
        // std::cout << "DiscardedEigValues: "
        // << discardedEigValues[0]           << std::endl;
        std::cout << "W_m : " << wm << std::endl;
        std::cout << "kbt" << kBT << std::endl;
        std::cout << "eAv: " << eAv << std::endl;
        std::cout << "eSqAv: " << eSqAv << std::endl;
        std::cout << "Entropy :" << std::log(lpart) + (eAv / kBT) << std::endl;
        std::cout << "partitionFunctionArray[i]: " << lpart << std::endl;
      }
      entropy[itm]      = std::log(lpart) + (eAv / kBT);
      specificHeat[itm] = (eSqAv - eAv * eAv) / (kBT * kBT);
    }
    // Save the data
    std::string hstr;
    pfile->write(temperatureArray, hstr + "temperatureArray");
    pfile->write(entropy, hstr + "entropy");
    pfile->write(specificHeat, hstr + "specificHeat");
    if (qsysOPerator != nullptr) {
      pfile->write(qsysOPeratorExpValue, hstr + "qsysOPeratorExpValue");
    }
    // pfile->write(vecPartitions, "vecPartitions");
  }

private:
  double maxExpNumber = 100;
  double regulateExp(double x) {
    if (x > 100) {
      return std::exp(100);
    }
    if (x < -100) {
      // return 0;
      return std::exp(-100);
    }
    return std::exp(x);
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
