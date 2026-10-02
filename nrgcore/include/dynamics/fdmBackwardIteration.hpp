#pragma once
#include "nrgcore/nrgData.hpp"
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
/**
 * @class fdmBackwardIteration
 * @brief Propagates the reduced density matrix through backward NRG iterations.
 *
 * This helper class tracks the kept/discarded state indices for the current
 * Wilson shell and builds the density matrix used to evaluate FDM-based
 * spectral quantities such as local Green's functions and static responses.
 *
 * @tparam nrgcore_type Type of the underlying NRG core object.
 */
template <typename nrgcore_type> class fdmBackwardIteration {
public:
  /// @brief Pointer to the active NRG core object used by the backward iteration.
  nrgcore_type *nrgObject;
  /// @brief Temperature of the NRG system in the FDM formalism.
  double        kBT{0};

  /**
   * @brief Construct an FDM backward-iteration helper for a given NRG core.
   *
   * @param t_nrgObject Pointer to the NRG model state used for the current
   * calculation.
   */
  explicit fdmBackwardIteration(nrgcore_type *t_nrgObject) {
    setup(t_nrgObject);
    this->clearKeptIndex();
  }

  /**
   * @brief Bind the helper to a specific NRG core object.
   *
   * This resets the iteration state to the initial backward-iteration regime and
   * prepares the object to work with the supplied model data.
   *
   * @param t_nrgObject Pointer to the NRG core object.
   */
  void setup(nrgcore_type *t_nrgObject) {
    lastiteration = true;
    nrgObject     = t_nrgObject;
  }

  /**
   * @brief Compute the density-matrix ingredients for the current shell.
   *
   * The routine updates the kept-state index set, reconstructs the density
   * matrix for the current shell, and reduces it to the impurity sector in the
   * order required for the backward iteration.
   *
   * @param energyScale Energy scale associated with the current NRG iteration,
   * typically proportional to \f$\Lambda^{-(N-1)/2}\f$.
   */
  void calcSpectrum(double energyScale) {
    // Clear the operator
    setCurrentIndex();
    // Order of these functions are important
    setRhoZero(energyScale);
    // rhoDotOperators();
    setReduceDensityMatrix();
  }

  /**
   * @brief Build the reduced density matrix for the impurity sector.
   *
   * The full density matrix is rotated into the instantaneous eigenbasis and
   * then traced over the bath degrees of freedom to produce the reduced density
   * matrix stored in `reducedRho` for the current shell.
   */
  void setReduceDensityMatrix() {
    // Set reducedRho
    reducedRho.clear();
    for (size_t i = 0; i < nrgObject->pre_sysmQ.size(); i++) {
      reducedRho.push_back( //
          qmatrix<>(nrgObject->eigenvaluesQ_kept_indices[i].size(),
                    nrgObject->eigenvaluesQ_kept_indices[i].size(), 0));
    }
    // Rotate the eigen basis
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      // U. rhoZero . U.T : TODO: check
      rhoZero[i] = nrgObject->current_hamiltonQ[i].dot(
          rhoZero[i].dot(nrgObject->current_hamiltonQ[i].cTranspose()));
    }
    // Set reducedRho
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      size_t kidx = 0;
      for (auto kindex : nrgObject->coupled_nQ_index[i]) {
        auto ii = kindex / nrgObject->nq_bath.size(); // impurity nqi index
        auto bb = kindex % nrgObject->nq_bath.size(); // bath nqi index
        // create previous bath id matrix
        for (size_t it : nrgObject->eigenvaluesQ_kept_indices[ii]) {
          for (size_t it_p : nrgObject->eigenvaluesQ_kept_indices[ii]) {
            double aa{0};
            for (size_t il = 0; il < nrgObject->bath_eigenvaluesQ[bb].size();
                 il++) {
              aa += rhoZero[i].at(
                  kidx + it +
                      (nrgObject->eigenvaluesQ_kept_indices[ii].size() * il),
                  kidx + it_p +
                      (nrgObject->eigenvaluesQ_kept_indices[ii].size() * il));
            }
            reducedRho[ii].at(it, it_p) += aa;
          }
        }
        kidx += nrgObject->eigenvaluesQ_kept_indices[ii].size() *
                nrgObject->bath_eigenvaluesQ[bb].size();
      }
      // End of matrix generation.
    }
    // once the reduced density matrix is defined we
    // set the last iteration flag for the next shell.
    lastiteration = false;
  }

  /**
   * @brief Compute the local partition function and ground-state energy.
   *
   * The routine identifies the lowest-energy state in each shell and assigns
   * Boltzmann weights to the states that are considered degenerate within the
   * numerical tolerance defined by `energyErrorBar`.
   */
  void setLocalPartitionFunction() {
    localGroundStateEnergy = 0;
    localPartitionFunction = 0; // Ground state degenarecy
    for (const auto &aa : nrgObject->eigenvaluesQ) {
      if (aa.size() != 0) {
        double result          = *std::min_element(aa.begin(), aa.end());
        localGroundStateEnergy = std::min(localGroundStateEnergy, result);
      }
    }
    for (size_t i = 0; i < nrgObject->eigenvaluesQ.size(); i++) {
      for (size_t ie = 0; ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
        double energy =
            std::fabs(nrgObject->eigenvaluesQ[i][ie] - localGroundStateEnergy);
        if (energy < energyErrorBar) {
          localPartitionFunction += 1.;
        }
      }
    }
    std::cout << "localGroundStateEnergy" << localGroundStateEnergy
              << " localPartitionFunction: " << localPartitionFunction
              << std::endl;
    BoltzmannFactor = nrgObject->eigenvaluesQ;
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      for (size_t ie = 0; ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
        double energy =
            std::fabs(nrgObject->eigenvaluesQ[i][ie] - localGroundStateEnergy);
        if (energy < energyErrorBar) {
          BoltzmannFactor[i][ie] = 1. / localPartitionFunction;
        } else {
          BoltzmannFactor[i][ie] = 0;
        }
      }
    }
  }

  /**
   * @brief Contract the density matrix with a set of static operators.
   *
   * @param bOperator Pointer to the operator set used in the contraction.
   * @return Vector of scalar traces \f$\mathrm{Tr}[\rho B]\f$ for each operator.
   */
  auto rhoDotStaticOperators(std::vector<qOperator> *bOperator) {
    // timer               t1("rhoDotStaticOperators");
    std::vector<double> specSum(bOperator->size(), 0.0);
    for (size_t ip = 0; ip < bOperator->size(); ip++) {
      for (size_t i = 0; i < nrgObject->eigenvaluesQ.size(); i++) {
        size_t kpdim = nrgObject->eigenvaluesQ[i].size();
        auto sys_opr_opt = (*bOperator)[ip].get(i, i);
        if (sys_opr_opt) {
          auto *sys_opr = sys_opr_opt.value();
          for (auto iv : currentKeptIndex[i]) {
            for (auto iv_p : currentKeptIndex[i]) {
              sys_opr->at(iv, iv_p) = 0;
            }
          }
          for (size_t iv = 0; iv < kpdim; iv++) {
            for (size_t iv_p = 0; iv_p < kpdim; iv_p++) {
              specSum[ip] += sys_opr->at(iv, iv_p) * rhoZero[i].at(iv_p, iv);
            }
          }
        }
      }
    }
    return specSum;
  }

  /**
   * @brief Construct the density matrix from a provided set of Boltzmann factors.
   *
   * The discarded states are filled with the supplied weights while the kept
   * states are overwritten with the reduced density matrix from the previous
   * Wilson shell whenever the backward iteration is not in the final step.
   *
   * @param tBoltzmannFactor Boltzmann factors of the form
   * \f$\exp(-\beta E_n)\f$ for each shell state.
   */
  void setRhoZero(const std::vector<std::vector<double>> &tBoltzmannFactor) {
    rhoZero.clear();
    double rhoTrace = 0;
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      size_t kpdim = nrgObject->eigenvaluesQ[i].size();
      qmatrix<> tmat(kpdim, kpdim, 0);
      for (size_t ie = currentKeptIndex[i].size();
           ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
        tmat(ie, ie) = tBoltzmannFactor[i][ie];
      }
      if (!lastiteration) {
        for (auto ik : currentKeptIndex[i]) {
          for (auto ikp : currentKeptIndex[i]) {
            tmat(ik, ikp) = reducedRho[i](ik, ikp);
          }
        }
      }
      rhoTrace += tmat.trace();
      rhoZero.push_back(tmat);
    }
    std::cout << "NRG Itr: " << nrgObject->nrg_iterations_cnt
              << "rhoTrace: " << rhoTrace << std::endl;
  }

  /**
   * @brief Construct the density matrix using the internally stored Boltzmann
   * factors and reduced density matrix.
   *
   * This overload is used when the class has already computed the local
   * partition function in the current shell.
   */
  void setRhoZero() {
    if (lastiteration) {
      setLocalPartitionFunction();
    }
    double rhoTrace{0};
    vecPartitions.push_back(localPartitionFunction);
    rhoZero.clear();
    std::cout << "Size : " << reducedRho.size() << " "
              << currentKeptIndex.size() << std::endl;
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      size_t kpdim = nrgObject->eigenvaluesQ[i].size();
      qmatrix<> tmat(kpdim, kpdim, 0);
      if (lastiteration) {
        for (size_t ie = 0; ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
          tmat(ie, ie) = BoltzmannFactor[i][ie];
        }
      }
      if (!lastiteration) {
        for (auto ik : currentKeptIndex[i]) {
          for (auto ikp : currentKeptIndex[i]) {
            tmat(ik, ikp) = reducedRho[i](ik, ikp);
          }
        }
      }
      rhoTrace += tmat.trace();
      rhoZero.push_back(tmat);
    }
    std::cout << "NrgItr: " << nrgObject->nrg_iterations_cnt
              << "rhoTrace: " << rhoTrace << std::endl;
  }

  /**
   * @brief Set the temperature used for the calculation.
   *
   * @param mkBT Temperature in the same units as the NRG energy scale.
   */
  void setTemperature(double mkBT) { kBT = mkBT; }

  /**
   * @brief Update the kept states for the current shell.
   *
   * The implementation reuses the previous kept-state set during later
   * iterations and initializes an empty set when the current shell is the last
   * one in the backward pass.
   */
  void setCurrentIndex() {
    if (lastiteration) {
      for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
        currentKeptIndex.emplace_back();
      }
    } else {
      currentKeptIndex = previoudKeptIndex;
    }
    previoudKeptIndex = nrgObject->eigenvaluesQ_kept_indices;
  }

  /**
   * @brief Reset the stored shell information and clear the density-matrix cache.
   *
   * This is useful when restarting the backward iteration from the final shell.
   */
  void clearKeptIndex() {
    lastiteration = true;
    BoltzmannFactor.clear();
    vecPartitions.clear();
    rhoZero.clear();
    reducedRho.clear();
    currentKeptIndex.clear();
    previoudKeptIndex.clear();
  }

private:
  std::vector<std::vector<double>> BoltzmannFactor;
  double                           localGroundStateEnergy{0};
  double localPartitionFunction{0}; // Ground state degenarecy
  std::vector<std::vector<size_t>> previoudKeptIndex;
  double spWeightErrorBar{1e-20}; // We dont care for the lower value
  std::vector<qmatrix<>> reducedRho;
  std::vector<double>    vecPartitions;

public: // Give access for openchain class
  /// @brief Kept-state indices for the current Wilson shell.
  std::vector<std::vector<size_t>> currentKeptIndex;
  /// @brief Density matrices for the current shell in the active basis.
  std::vector<qmatrix<>>           rhoZero;
  /// @brief True when the current shell is the final backward-iteration step.
  bool                             lastiteration{true};
  /// @brief Numerical tolerance used when identifying degenerate shell states.
  double                           energyErrorBar{1e-5};
};
