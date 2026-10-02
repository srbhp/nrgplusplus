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
/**
 * @class fdmThermodynamics
 * @brief Evaluate thermodynamic quantities from discarded-shell FDM data.
 *
 * This helper accumulates the discarded eigenvalues and the full spectrum for
 * each NRG iteration, builds the corresponding thermodynamic weights, and stores
 * the resulting entropy, specific heat, and optionally observable expectation
 * values for a supplied temperature grid.
 *
 * @tparam nrgcore_type Type of the underlying NRG core object.
 */
template <typename nrgcore_type>
class fdmThermodynamics {
  nrgcore_type                                      *nrg_object;
  int                                                nrgMaxIterations;
  size_t                                             bathDimension{0};
  std::vector<qOperator>                            *qsysOPerator{nullptr};
  std::vector<std::vector<std::map<size_t, double>>> qsysOPeratorValue;
  std::vector<double>                                temperatureArray;

public:
  /**
   * @brief Construct a thermodynamic evaluator for a given NRG run.
   *
   * @param t_nrg_object Pointer to the NRG model object.
   * @param tArray Temperature grid used to evaluate entropy and specific heat.
   */
  explicit fdmThermodynamics(nrgcore_type *t_nrg_object,
                             std::vector<double> tArray)
      : nrg_object(t_nrg_object),
        nrgMaxIterations(t_nrg_object->nrg_iterations_cnt),
        temperatureArray(std::move(tArray)) {
    for (auto &it : t_nrg_object->bath_eigenvaluesQ) {
      bathDimension += it.size();
    }
    std::cout << "bathDimension: " << bathDimension << std::endl;
  }

  /**
   * @brief Register a system operator whose expectation value should be stored.
   *
   * @param qsysOPerator_ Pointer to the vector of operators to evaluate.
   */
  void setSystemOperator(std::vector<qOperator> *qsysOPerator_) {
    qsysOPerator = qsysOPerator_;
  }

  /**
   * @brief Accumulate one thermodynamic shell contribution for the current NRG
   * iteration.
   *
   * The routine records the discarded eigenvalues, the full eigenvalue set, and
   * the diagonal matrix element of the registered operator on the discarded
   * states for later thermal averaging.
   *
   * @param energyRescale Energy scaling factor applied to the stored eigenvalues.
   */
  void calcThermodynamics(double energyRescale) {
    std::cout << "energyRescale" << energyRescale << std::endl;
    setCurrentKeptIndex();
    std::vector<double> eigVal;
    for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
      for (size_t j = currentKeptIndex[i].size();
           j < nrg_object->eigenvaluesQ[i].size(); j++) {
        eigVal.push_back(nrg_object->eigenvaluesQ[i][j] * energyRescale);
      }
    }
    std::cout << "Number of eigen Values: " << eigVal.size() << std::endl;
    discardedEigValues.push_back(eigVal);
    eigVal.clear();
    for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
      for (size_t j = 0; j < nrg_object->eigenvaluesQ[i].size(); j++) {
        eigVal.push_back(nrg_object->eigenvaluesQ[i][j] * energyRescale);
      }
    }
    allEigValues.push_back(eigVal);
    if (qsysOPerator != nullptr) {
      std::vector<std::map<size_t, double>> qvalue(qsysOPerator->size());
      for (size_t iq = 0; iq < qsysOPerator->size(); iq++) {
        size_t icounter{0};
        for (int i = 0; i < nrg_object->eigenvaluesQ.size(); i++) {
          auto qopr = qsysOPerator->at(iq).get(i, i);
          if (qopr) {
            auto *qmat = qopr.value();
            for (size_t j = currentKeptIndex[i].size();
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

  /**
   * @brief Update the kept-state index set for the current shell.
   */
  void setCurrentKeptIndex() {
    if (currentKeptIndex.empty()) {
      for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
        currentKeptIndex.emplace_back();
      }
    } else {
      currentKeptIndex = previoudKeptIndex;
    }
    previoudKeptIndex = nrg_object->eigenvaluesQ_kept_indices;
  }

  /**
   * @brief Compute and write the final thermodynamic observables to a file.
   *
   * The function builds the shell weights, entropy, and specific heat as a
   * function of the supplied temperature grid and writes the results to the
   * provided file wrapper.
   *
   * @tparam filetype Type of the output file object exposing a `write` method.
   * @param pfile Pointer to the output file handle.
   */
  template <typename filetype> void saveFinalData(filetype *pfile) {
    std::vector<double> groundStateEnergyArr(allEigValues.size(), 0);
    for (int it = 0; it < allEigValues.size(); it++) {
      double result =
          *std::min_element(allEigValues[it].begin(), allEigValues[it].end());
      std::transform(discardedEigValues[it].begin(),
                     discardedEigValues[it].end(),
                     discardedEigValues[it].begin(),
                     [result](double x) { return x - result; });
      groundStateEnergyArr[it] = result;
      std::cout << it << "  " << result << "  "
                << *std::max_element(allEigValues[it].begin(),
                                     allEigValues[it].end())
                << std::endl;
    }
    for (size_t it = 0; it < groundStateEnergyArr.size() - 1; it++) {
      groundStateEnergyArr[it] += groundStateEnergyArr[it + 1];
    }
    std::cout << "|Ground State Energy|: " << groundStateEnergyArr << std::endl;
    double gren = groundStateEnergyArr[0];
    std::transform(groundStateEnergyArr.begin(), groundStateEnergyArr.end(),
                   groundStateEnergyArr.begin(),
                   [gren](double aa) { return aa - gren; });
    std::cout << "|Ground State Energy|: " << groundStateEnergyArr << std::endl;
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
      auto   nrgCount = nrgMaxIterations;
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
      double lpart = std::accumulate(wm.begin(), wm.end(), 0.0);
      for (auto &aa : wm) {
        aa = aa / lpart;
      }
      std::cout << "Done 2" << std::endl;
      double eAv = 0;
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
      if (itm == 0) {
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
    std::string hstr;
    pfile->write(temperatureArray, hstr + "temperatureArray");
    pfile->write(entropy, hstr + "entropy");
    pfile->write(specificHeat, hstr + "specificHeat");
    if (qsysOPerator != nullptr) {
      pfile->write(qsysOPeratorExpValue, hstr + "qsysOPeratorExpValue");
    }
  }

private:
  /**
   * @brief Clamp the exponent argument to avoid overflow in thermal weights.
   *
   * @param x Exponential argument to regularize.
   * @return Stabilized value of `exp(x)`.
   */
  double maxExpNumber = 100;
  double regulateExp(double x) {
    if (x > 100) {
      return std::exp(100);
    }
    if (x < -100) {
      return std::exp(-100);
    }
    return std::exp(x);
  }

  /**
   * @brief Previously selected kept-state indices for the preceding shell.
   */
  std::vector<std::vector<size_t>> previoudKeptIndex;

  /**
   * @brief Discarded eigenvalues accumulated for each NRG iteration.
   */
  std::vector<std::vector<double>> discardedEigValues;

  /**
   * @brief Full eigenvalue list accumulated for each NRG iteration.
   */
  std::vector<std::vector<double>> allEigValues;

public:
  /**
   * @brief Kept-state indices for the current shell.
   */
  std::vector<std::vector<size_t>> currentKeptIndex;
};
