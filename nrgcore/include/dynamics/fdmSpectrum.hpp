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
#include <stdexcept>
#include <tuple>
#include <vector>
/**
 * @class fdmSpectrum
 * @brief Computes FDM spectral weights from a full density matrix.
 *
 * This helper accumulates the operator-weighted spectral contributions for each
 * frequency bin and keeps track of the kept/discarded state structure needed to
 * propagate the reduced density matrix across NRG iterations.
 *
 * @tparam nrgcore_type Type of the underlying NRG core object.
 */
template <typename nrgcore_type> class fdmSpectrum {
  nrgcore_type *nrg_object;
  double        kBT{0}; // Temperature of nrg system i.e., in FDM formalism

public:
  /**
   * @brief Construct a spectral evaluator bound to an NRG core object.
   *
   * @param t_nrg_object Pointer to the NRG core whose eigenstates and operators
   * will be used to build the FDM spectral weights.
   */
  explicit fdmSpectrum(nrgcore_type *t_nrg_object) { setup(t_nrg_object); }

  /**
   * @brief Store the NRG core pointer and initialize the iteration state.
   *
   * @param t_nrg_object Pointer to the current NRG model object.
   */
  void setup(nrgcore_type *t_nrg_object) {
    lastiteration           = true;
    nrg_object              = t_nrg_object;
    globalGroundStateEnergy = 0; //= nrg_object->all_eigenvalue[0];
  }

  /**
   * @brief Compute one spectral step for the current Wilson shell.
   *
   * The routine updates the kept-state index set, initializes the density matrix,
   * evaluates the operator contraction, and reduces the result to the impurity
   * sector.
   *
   * @param energyScale Rescaling factor applied to the energy grid.
   */
  void calcSpectrum(double energyScale) {
    energyRescale = energyScale;
    setCurrentIndex();
    setRhoZero();
    rhoDotOperators();
    setReduceDensityMatrix();
    lastiteration = false;
  }

  /**
   * @brief Build the reduced density matrix from the current full density matrix.
   *
   * This traces out the bath degrees of freedom after rotating `rhoZero` into the
   * current eigenbasis, leaving the reduced matrix for the impurity or kept
   * subspace in `reducedRho`.
   */
  void setReduceDensityMatrix() {
    reducedRho.clear();
    for (size_t i = 0; i < nrg_object->pre_sysmQ.size(); i++) {
      reducedRho.push_back( //
          qmatrix<>(nrg_object->eigenvaluesQ_kept_indices[i].size(),
                    nrg_object->eigenvaluesQ_kept_indices[i].size(), 0));
    }
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      rhoZero[i] = nrg_object->current_hamiltonQ[i].dot(
          rhoZero[i].dot(nrg_object->current_hamiltonQ[i].cTranspose()));
    }
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      size_t kidx = 0;
      for (auto kindex : nrg_object->coupled_nQ_index[i]) {
        auto ii = kindex / nrg_object->nq_bath.size();
        auto bb = kindex % nrg_object->nq_bath.size();
        for (size_t it : nrg_object->eigenvaluesQ_kept_indices[ii]) {
          for (size_t it_p : nrg_object->eigenvaluesQ_kept_indices[ii]) {
            double aa{0};
            for (size_t il = 0; il < nrg_object->bath_eigenvaluesQ[bb].size();
                 il++) {
              aa += rhoZero[i].at(
                  kidx + it +
                      (nrg_object->eigenvaluesQ_kept_indices[ii].size() * il),
                  kidx + it_p +
                      (nrg_object->eigenvaluesQ_kept_indices[ii].size() * il));
            }
            reducedRho[ii].at(it, it_p) += aa;
          }
        }
        kidx += nrg_object->eigenvaluesQ_kept_indices[ii].size() *
                nrg_object->bath_eigenvaluesQ[bb].size();
      }
    }
  }

  /**
   * @brief Compute the local partition function and Boltzmann weights.
   *
   * The routine identifies the ground-state energy for the current shell and
   * assigns a unit weight to states within the degeneracy tolerance defined by
   * `energyErrorBar`.
   */
  void setLocalPartitionFunction() {
    localGroundStateEnergy = 0;
    localPartitionFunction = 0; // Ground state degenarecy
    for (auto &aa : nrg_object->eigenvaluesQ) {
      if (aa.size() != 0) {
        double result          = *std::min_element(aa.begin(), aa.end());
        localGroundStateEnergy = std::min(localGroundStateEnergy, result);
      }
    }
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
        double energy =
            std::fabs(nrg_object->eigenvaluesQ[i][ie] - localGroundStateEnergy);
        if (energy < energyErrorBar) {
          localPartitionFunction += 1.;
        }
      }
    }
    BoltzmannFactor = nrg_object->eigenvaluesQ;
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
        double energy =
            std::fabs(nrg_object->eigenvaluesQ[i][ie] - localGroundStateEnergy);
        if (energy < energyErrorBar) {
          BoltzmannFactor[i][ie] = 1. / localPartitionFunction;
        } else {
          BoltzmannFactor[i][ie] = 0;
        }
      }
    }
  }

  /**
   * @brief Boltzmann weights for the states in the current shell.
   */
  std::vector<std::vector<double>> BoltzmannFactor;

  /**
   * @brief Build the full density matrix for the current shell.
   *
   * At the final iteration the density matrix is initialized from the local
   * partition function; otherwise it is overwritten with the previously reduced
   * density matrix on the kept states.
   */
  void setRhoZero() {
    if (lastiteration) {
      setLocalPartitionFunction();
    }
    double rhoTrace{0};
    vecPartitions.push_back(localPartitionFunction);
    rhoZero.clear();
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      size_t kpdim = nrg_object->eigenvaluesQ[i].size();
      qmatrix<> tmat(kpdim, kpdim, 0);
      if (lastiteration) {
        for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
          tmat(ie, ie) = BoltzmannFactor[i][ie];
        }
      }
      if (!reducedRho.empty()) {
        for (auto ik : currentKeptIndex[i]) {
          for (auto ikp : currentKeptIndex[i]) {
            tmat(ik, ikp) = reducedRho[i](ik, ikp);
          }
        }
      }
      rhoTrace += tmat.trace();
      rhoZero.push_back(tmat);
    }
    std::cout << "rhoTrace: " << rhoTrace << std::endl;
  }

  /**
   * @brief Accumulate the spectral weight contributions for the current shell.
   *
   * The routine evaluates $\rho B$ and $B \rho$ contractions and accumulates the
   * resulting positive and negative frequency contributions in the internal
   * weight arrays.
   */
  void rhoDotOperators() {
    double specSum = 0.0;
    for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
      for (size_t j = 0; j < nrg_object->current_sysmQ.size(); j++) {
        size_t kpdim   = nrg_object->eigenvaluesQ[i].size();
        size_t kpdim_p = nrg_object->eigenvaluesQ[j].size();
        for (size_t ip = 0; ip < bOperator->size(); ip++) {
          auto sys_opr_opt = (*bOperator)[ip].get(i, j);
          if (sys_opr_opt) {
            auto     *sys_opr = sys_opr_opt.value();
            qmatrix<> aMatrix;
            if (aOperator == nullptr) {
              aMatrix = sys_opr->cTranspose();
            } else {
              auto a_opr_opt = (*aOperator)[ip].get(i, j);
              if (!a_opr_opt) {
                continue;
              }
              aMatrix = *a_opr_opt.value();
            }
            for (auto iv : currentKeptIndex[i]) {
              for (auto iv_p : currentKeptIndex[j]) {
                aMatrix.at(iv_p, iv) = 0;
              }
            }
            auto rhoA = sys_opr->dot(rhoZero[j]);
            auto ARho = rhoZero[i].dot(*sys_opr);
            for (size_t iv = 0; iv < kpdim; iv++) {
              for (size_t iv_p = 0; iv_p < kpdim_p; iv_p++) {
                double aa{0};
                double bbv{0};
                aa          = rhoA(iv, iv_p) * aMatrix.at(iv_p, iv);
                bbv         = ARho(iv, iv_p) * aMatrix.at(iv_p, iv);
                int tmindex = int(
                    std::log(1 +
                             (std::fabs(energyRescale *
                                        (nrg_object->eigenvaluesQ[i][iv] -
                                         nrg_object->eigenvaluesQ[j][iv_p])) /
                              minEnergy)) *
                    delE);
                if (tmindex < energyPts && tmindex >= 0) {
                  if (std::fabs(aa) > spWeightErrorBar) {
                    positiveWeight[ip][tmindex] += std::fabs(aa);
                    specSum += std::fabs(aa);
                  }
                  if (std::fabs(bbv) > spWeightErrorBar) {
                    negativeWeight[ip][tmindex] += std::fabs(bbv);
                    specSum += std::fabs(bbv);
                  }
                }
              }
            }
          }
        }
      }
    }
    std::cout << nrg_object->nrg_iterations_cnt << "specSum: " << specSum
              << " Scale: " << energyRescale << std::endl;
  }

  /**
   * @brief Set the temperature used by the FDM evaluation.
   *
   * @param at Temperature value in the same units as the NRG eigenvalues.
   */
  void setTemperature(double at) { kBT = at; }

  /**
   * @brief Update the current kept-state indices for the active shell.
   */
  void setCurrentIndex() {
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
   * @brief Register the creation and annihilation operators used in the
   * spectral weight calculation.
   *
   * @param bopr Pointer to the creation-like operator set.
   * @param aopr Optional pointer to the annihilation-like operator set. If null,
   * the same operator set is used for both directions.
   * @throws std::invalid_argument If the input operator pointer is null or the
   * operator counts mismatch.
   */
  void setOperator(std::vector<qOperator> *bopr,
                   std::vector<qOperator> *aopr = nullptr) {
    if (bopr == nullptr) {
      throw std::invalid_argument("The creation-like operator cannot be null");
    }
    if (aopr != nullptr && aopr->size() != bopr->size()) {
      throw std::invalid_argument(
          "Creation-like and annihilation-like operator counts differ");
    }
    aOperator = aopr;
    bOperator = bopr;
    for (size_t i = 0; i < bOperator->size(); i++) {
      positiveWeight.emplace_back(energyPts, 0);
      negativeWeight.emplace_back(energyPts, 0);
    }
  }

  /**
   * @brief Save the computed spectral weights to a file object.
   *
   * @tparam filetype Type of the output file wrapper exposing a `write` method.
   * @param pfile Pointer to the file object receiving the spectral data.
   */
  template <typename filetype> void saveFinalData(filetype *pfile) {
    std::vector<double> energyPoints(energyPts, 0);
    for (int i = 0; i < energyPts; i++) {
      energyPoints[i] = (minEnergy * std::exp(i * 1. / delE));
    }
    std::string hstr = "GreenFn";
    pfile->write(energyPoints, hstr + "EnergyPoints");
    pfile->write(positiveWeight, hstr + "PositiveWeight");
    pfile->write(negativeWeight, hstr + "NegativeWeight");
    for (size_t i = 0; i < bOperator->size(); i++) {
      positiveWeight[i].clear();
      negativeWeight[i].clear();
    }
  }

private:
  double localGroundStateEnergy{0};
  double localPartitionFunction{0}; // Ground state degenarecy
  std::vector<std::vector<size_t>> previoudKeptIndex;
  double spWeightErrorBar{1e-20}; // We dont care for the lower value
  double globalGroundStateEnergy{0};
  std::vector<qmatrix<>> reducedRho;
  std::vector<double>    vecPartitions;
  std::vector<std::vector<double>> positiveWeight;
  std::vector<std::vector<double>> negativeWeight;
  int    energyPts = 100000;
  double maxEnergy = 10; // One decade more
  double minEnergy = 1e-10;
  double delE      = (energyPts - 1.0) / (std::log(maxEnergy / minEnergy));

public: // Give access for openchain class
  std::vector<std::vector<size_t>> currentKeptIndex;
  std::vector<qmatrix<>>           rhoZero;
  std::vector<qOperator>          *aOperator{};
  std::vector<qOperator>          *bOperator{};
  bool                             lastiteration{true};
  double                           energyRescale{1};
  double                           energyErrorBar{1e-5};
};
