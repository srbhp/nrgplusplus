#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>
class freeModel : public fermionBasis {
  /** This is single channel free model class
   * to calculate the entropy/ specific heat
   * and other static quantitities.
   *
   *
   *
   */
public:
  /**
   * @brief Construct a new rabiSpinless object
   *
   * @param JKondo: J value for the Kondo interaction
   * @param spinS: Spin value i.e,. 1/2 or 3/2. Integer Spin may not work
   */
  freeModel() : fermionBasis(2, fermionBasis::chargeAndSpin) {
    createBasis(); // create the basis in nstates x nstates
  }
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; }
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }
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
