#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>

/**
 * @brief Free-fermion reference model used for static thermodynamic benchmarks.
 *
 * This model builds a non-interacting fermionic basis with charge and spin
 * quantum numbers and diagonalizes the local Hamiltonian in each symmetry sector.
 * It is useful as a reference when comparing with interacting models or computing
 * entropy and specific-heat contributions in a noninteracting limit.
 */
class freeModel : public fermionBasis {
public:
  /**
   * @brief Construct the free-fermion benchmark model.
   *
   * Builds the fermion basis, identifies the conserved quantum numbers, and
   * computes the diagonalized eigenvalues for each block.
   */
  freeModel() : fermionBasis(2, fermionBasis::chargeAndSpin) {
    createBasis();
  }

  /**
   * @brief Return the list of quantum numbers for each symmetry sector.
   *
   * @return Vector of charge/spin sector labels.
   */
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; }

  /**
   * @brief Return the eigenvalues of the Hamiltonian in each sector.
   *
   * @return Block-diagonal eigenvalues associated with each quantum-number sector.
   */
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }

  /**
   * @brief Return the fermionic parity factor for each sector.
   *
   * @return Vector of parity signs used in the fermionic sector structure.
   */
  [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }
  //
  std::vector<std::vector<double>> eigenvalues_Q;
  std::vector<double>              chi_Q;
  std::vector<std::vector<int>>    n_Q;
  //    ########################################
private:
  void createBasis() {
    //
    createFermionBasis(2);
    std::cout << "FermionBasis Size" << fermionOprMat.size() << "\n";
    auto Hamiltonian = fermionOprMat[0] * 0; // Set Hamiltonian to Zero
    //
    create_QuantumNspinCharge();
    create_Block_structure();
    // ####################################################################
    n_Q = get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    //
    // set foperator
    auto h_blocked = get_block_Hamiltonian(Hamiltonian);
    //    std::cout << "h_blocked: " << h_blocked << std::endl;
    //    std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
    }
    //    std::cout << "Eigenvalues: " << eigenvalues_Q << std::endl;
    // TODO(sp): rotate the f operator
    // ####################################################################
    f_dag_operator = get_block_operators({fermionOprMat[0], fermionOprMat[1]});
    std::cout << "f_dag_operators: " << f_dag_operator.size() << std::endl;
    std::vector<qOperator> topr(f_dag_operator.size(), qOperator());
    for (size_t ip = 0; ip < f_dag_operator.size(); ip++) {
      for (size_t i = 0; i < n_Q.size(); i++) {
        for (size_t j = 0; j < n_Q.size(); j++) {
          auto tfopr = f_dag_operator[ip].get(i, j);
          if (tfopr) {
            topr[ip].set((h_blocked.get(i, i))
                             .value()
                             ->cTranspose()
                             .dot(*tfopr.value())
                             .dot(*(h_blocked.get(j, j)).value()),
                         i, j);
          }
        }
      }
    }
    f_dag_operator = topr;
  }
  //    ######################################
};
