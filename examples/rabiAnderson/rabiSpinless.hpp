#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>
class rabiSpinless : public fermionBasis {
  // See PHYSICAL REVIEW B 101, 085110 (2020)
  /** This class is for a single orbital with spin up and down
   * f operator. SIAM can be made entirely from this class.
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
  rabiSpinless(double JKondo, double omega)
      : fermionBasis(5, fermionBasis::chargeOnly) {
    int dof = 1 + // for the trion State (up)
              2 + // for the electron state
              2;  // for the wilson site
    // Dimension of the system
    double nstates = std::pow(2, dof);
    // Dimension of the system
    std::cout << "nstates: " << nstates << "\n"
              << "dof: " << dof << "\n";
    createBasis(JKondo, omega); // create the basis in nstates x nstates
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
  void createBasis(double JKondo, double Omega) {
    //
    // createFermionBasis();
    std::cout << "FermionBasis Size" << fermionOprMat.size() << "\n";
    // Create the spinOperator electron[Impurity]
    auto iSpinx = (fermionOprMat[3].dot(fermionOprMat[2].cTranspose()) +
                   fermionOprMat[2].dot(fermionOprMat[3].cTranspose())) *
                  0.5;
    auto iSpiny = (fermionOprMat[3].dot(fermionOprMat[2].cTranspose()) -
                   fermionOprMat[2].dot(fermionOprMat[3].cTranspose())) *
                  0.5;
    auto iSpinz = (fermionOprMat[2].dot(fermionOprMat[2].cTranspose()) -
                   fermionOprMat[3].dot(fermionOprMat[3].cTranspose())) *
                  0.5;
    // Create the spinOperator electron[Wilson Site ]
    auto wSpinx = (fermionOprMat[1].dot(fermionOprMat[0].cTranspose()) +
                   fermionOprMat[0].dot(fermionOprMat[1].cTranspose())) *
                  0.5;
    auto wSpiny = (fermionOprMat[1].dot(fermionOprMat[0].cTranspose()) -
                   fermionOprMat[0].dot(fermionOprMat[1].cTranspose())) *
                  0.5;
    auto wSpinz = (fermionOprMat[0].dot(fermionOprMat[0].cTranspose()) -
                   fermionOprMat[1].dot(fermionOprMat[1].cTranspose())) *
                  0.5;
    // ########################################################
    //    std::cout << "iSpinx" << iSpinx << std::endl;
    //    std::cout << "iSpiny: " << iSpiny << std::endl;
    //    std::cout << "iSpinz" << iSpinz << std::endl;
    //    //
    //    std::cout << "wSpinx: " << wSpinx << std::endl;
    //    std::cout << "wSpiny: " << wSpiny << std::endl;
    //    std::cout << "wSpinz: " << wSpinz << std::endl;
    // ############################################
    auto Hamiltonian =
        ((iSpinz.dot(wSpinz) + iSpinx.dot(wSpinx) //
          - iSpiny.dot(wSpiny)) * // Imaginary part is taken care here
         JKondo)
        // Trion coupling
        + ((fermionOprMat[4].dot(fermionOprMat[3].cTranspose()) +
            fermionOprMat[3].dot(fermionOprMat[4].cTranspose())) *
           Omega);
    // End
    //
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
    // TODO: rotate the f operator
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
