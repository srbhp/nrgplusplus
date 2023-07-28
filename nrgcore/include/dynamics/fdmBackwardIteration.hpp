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
 * @brief  This class is responsible for the backward iteration
 * of the NRG algorithm. This class has three useful functions
 * to calculate the energy dependent quantities such as $A(\omega)$.
 * These functions are :
 *   - 1. setCurrentIndex : set the indices  for the kept and Discarded states
 *   - 2. setRhoZero : Creates the density matrix of the current Wilson Chain.
 *   - 3. setReduceDensityMatrix: Reduce the density matrix after summing over
 * the Enviorentment degree.
 *
 *
 * @tparam nrgcore_type : Type of the nrgcore object
 * @param t_nrgObject : norgcore object created from a Impurity and bath class.
 * @return [TODO:return]
 */
template <typename nrgcore_type> // nrgcore_type is a type of
class fdmBackwardIteration {
public:
  nrgcore_type *nrgObject;
  double        kBT{0}; // Temperature of nrg system i.e., in FDM formalism
  // fdmBackwardIteration() {}
  explicit fdmBackwardIteration(nrgcore_type *t_nrgObject // nrgcore_type
  ) {
    setup(t_nrgObject);
    this->clearKeptIndex();
  }
  /**
   * @brief Explicitly set the `nrgobject`
   *
   * @param t_nrgObject
   */
  void setup(nrgcore_type *t_nrgObject) {
    lastiteration = true;
    nrgObject     = t_nrgObject;
  }
  /**
   * @brief
   * These function calls the following functions :
   *   - 1. setCurrentIndex : set the indices  for the kept and Discarded states
   *   - 2. setRhoZero : Creates the density matrix of the current Wilson Chain.
   *   - 3. setReduceDensityMatrix: Reduce the density matrix after summing over
   * the Enviorentment degree.
   *
   * @param energyScale: Energy scale of the currect NRG iteration. \f$
   * \Lambda^{-(N-1)/2} \f$
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
   * @brief Creates the reduced density matrix.
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
                      nrgObject->eigenvaluesQ_kept_indices[ii].size() * il,
                  kidx + it_p +
                      nrgObject->eigenvaluesQ_kept_indices[ii].size() * il);
            }
            reducedRho[ii].at(it, it_p) += aa;
          }
        }
        kidx += nrgObject->eigenvaluesQ_kept_indices[ii].size() *
                nrgObject->bath_eigenvaluesQ[bb].size();
      }
      // End of matrix generation.
    }
    // double trace = std::accumulate(
    //                    reducedRho.begin(), reducedRho.end(), 0.0,
    // [](double a, const qmatrix<> &b) {
    //     return a + b.trace();
    // });
    // std::cout << "NrgItr: " << nrgObject->nrg_iterations_cnt
    //           << "rhoTrace: " << trace << std::endl;
    // once the reduce density matrix is defined We
    // set the lastiteration to be false for the next iteration
    lastiteration = false;
  }
  /**
   * @brief Sets the Local partition function for the current Wilson chain.
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
    } //
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
    } //
  }
  /**
   * @brief $\rho \cdot B$ : Dot product between $\rho$ and B operator.
   *
   * @param bOperator: Pointer to a std::vector of  B operators.
   * @return returns the dot product value as std::vector
   */
  auto rhoDotStaticOperators(std::vector<qOperator> *bOperator) {
    // timer               t1("rhoDotStaticOperators");
    std::vector<double> specSum(bOperator->size(), 0.0);
    //
    for (size_t ip = 0; ip < bOperator->size(); ip++) {
      for (size_t i = 0; i < nrgObject->eigenvaluesQ.size(); i++) {
        size_t kpdim = nrgObject->eigenvaluesQ[i].size();
        // std::cout << "--------------------------";
        // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
        // TODO(sp): This
        auto sys_opr_opt = (*bOperator)[ip].get(i, i);
        if (sys_opr_opt) {
          auto *sys_opr = sys_opr_opt.value();
          // set kept-kept part operator
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
        // all the matrices are set
      }
    }
    // set the matrix elements
    // end of lm loop
    // End of matrix generation.
    // Rotate the c operator in the eigen basis
    return specSum;
  }
  /**
   * @brief Set the density matrix based on the BoltzmannFactors
   * and reduced density matrix of the previous Wilson Chain if available.
   *
   * @param tBoltzmannFactor: Boltzmann Factors  of the form of $exp(- \beta
   * E_n)$
   */
  void setRhoZero(const std::vector<std::vector<double>> &tBoltzmannFactor) {
    // Clear the operator
    // Order of these functions are important
    rhoZero.clear();
    double rhoTrace = 0;
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      size_t kpdim = nrgObject->eigenvaluesQ[i].size();
      // std::cout << "--------------------------";
      // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
      qmatrix<> tmat(kpdim, kpdim, 0);
      // Discarded states
      for (size_t ie = currentKeptIndex[i].size();
           ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
        tmat(ie, ie) = tBoltzmannFactor[i][ie];
      }
      // Kept states
      if (!lastiteration) { // Condition for the last Wilson site
        // std::cout << "Not lastiteration" << std::endl;
        //  Just Override the kept states
        for (auto ik : currentKeptIndex[i]) {
          for (auto ikp : currentKeptIndex[i]) {
            tmat(ik, ikp) = reducedRho[i](ik, ikp);
          }
        }
      }
      rhoTrace += tmat.trace();
      // Save the matrix
      rhoZero.push_back(tmat);
    }
    std::cout << "NRG Itr: " << nrgObject->nrg_iterations_cnt
              << "rhoTrace: " << rhoTrace << std::endl;
    // move the operator
  } // End of  update_system_operatorQ
  void setRhoZero() {
    // Clear the operator
    // Order of these functions are important
    // This function is called
    // for setting up current rho
    //
    //
    //
    if (lastiteration) {
      setLocalPartitionFunction();
    }
    double rhoTrace{0};
    // Calc partition function of the shell ::
    vecPartitions.push_back(localPartitionFunction);
    rhoZero.clear();
    std::cout << "Size : " << reducedRho.size() << " "
              << currentKeptIndex.size() << std::endl;
    for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
      size_t kpdim = nrgObject->eigenvaluesQ[i].size();
      // std::cout << "--------------------------";
      // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
      qmatrix<> tmat(kpdim, kpdim, 0);
      // Discarded states
      if (lastiteration) { // Only for the laast iterationT= 0
        for (size_t ie = 0; ie < nrgObject->eigenvaluesQ[i].size(); ie++) {
          tmat(ie, ie) = BoltzmannFactor[i][ie];
        }
      }
      // Kept states
      if (!lastiteration) { // Condition for the last Wilson site
        // Just Override the kept states
        // std::cout << reducedRho[i].size() << " " <<
        // currentKeptIndex[i].size()
        //           << std::endl;
        for (auto ik : currentKeptIndex[i]) {
          for (auto ikp : currentKeptIndex[i]) {
            tmat(ik, ikp) = reducedRho[i](ik, ikp);
          }
        }
      }
      rhoTrace += tmat.trace();
      // Save the matrix
      rhoZero.push_back(tmat);
    }
    std::cout << "NrgItr: " << nrgObject->nrg_iterations_cnt
              << "rhoTrace: " << rhoTrace << std::endl;
    // move the operator
  } // End of  update_system_operatorQ
  //
  /**
   * @brief Set the temperature of the system.
   *
   * @param mkBT: Temperature
   */
  void setTemperature(double mkBT) { kBT = mkBT; }
  void setCurrentIndex() {
    // This only work for the backward iteration
    if (lastiteration) {
      // every state is Discarded
      for (size_t i = 0; i < nrgObject->current_sysmQ.size(); i++) {
        currentKeptIndex.emplace_back();
      }
    } else {
      currentKeptIndex = previoudKeptIndex;
    }
    previoudKeptIndex = nrgObject->eigenvaluesQ_kept_indices;
  }
  /**
   * @brief Clear everything. Useful for the last iteration.
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
  std::vector<std::vector<size_t>> currentKeptIndex;
  std::vector<qmatrix<>>           rhoZero;
  bool                             lastiteration{true};
  double                           energyErrorBar{1e-5};
};
