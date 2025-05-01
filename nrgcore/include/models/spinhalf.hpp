#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>
/**
 * @class spinhalf
 * @brief Represents a single orbital with spin-up and spin-down operators.
 *
 * This model can be used to construct the Single Impurity Anderson Model (SIAM).
 * It includes methods to calculate eigenvalues, quantum numbers, and operators.
 */
class spinhalf {
  std::vector<qmatrix<>> f_dag_raw;
  double                 Uint;
  double                 epsilon_d;
  double                 magnetic_field;

public:
  /**
   * @brief Constructs the `spinhalf` model.
   *
   * @param teps Onsite energy of the SIAM.
   * @param tUint Coulomb interaction energy of the SIAM.
   * @param tmag Magnetic field of the SIAM (default is 0).
   */
  spinhalf(double teps, double tUint, double tmag = 0) // NOLINT
      : Uint(tUint), epsilon_d(teps), magnetic_field(tmag) {
    // set up the basis
    fermionBasis spinhalfBasis(2, fermionBasis::chargeAndSpin);
    // Number of fermion channels/spins
    f_dag_operator = spinhalfBasis.get_f_dag_operator();
    f_dag_raw      = spinhalfBasis.get_raw_f_dag_operator();
    // set f_operator
    // {
    //   auto f_raw = f_dag_raw;
    //   for (auto &aa : f_raw) {
    //     aa = aa.cTranspose();
    //   }
    //   f_operator = spinhalfBasis.get_block_operators(f_raw);
    // }
    // n_Q
    n_Q = spinhalfBasis.get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    // std::cout << "chi_Q: " << chi_Q << std::endl;
    // set up the Hamiltonian
    auto n_up   = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
    auto n_down = f_dag_raw[1].dot(f_dag_raw[1].cTranspose());
    // std::cout << "n_up" << n_up << "n_down" << n_down;
    auto H = (n_up + n_down) * epsilon_d        // onsite energy
             + (n_up - n_down) * magnetic_field // Magnetic Field
             + (n_up.dot(n_down)) * Uint;       // Columb Energy
    // get the hamiltonian  in the blocked structure
    auto h_blocked  = spinhalfBasis.get_block_Hamiltonian(H);
    doubleOccupancy = spinhalfBasis.get_block_operators({n_up.dot(n_down)});
    //  h_blocked.display();
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag(); // NOLINT
    }
    // TODO(sp): rotate the f operator
    // End of the constructor
  }
  // TODO(sp): Remove the copy of the objects
  /**
   * @brief Returns the quantum numbers of the basis states.
   *
   * @return A vector of quantum numbers.
   */
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; };

  /**
   * @brief Returns the eigenvalues of the Hamiltonian.
   *
   * @return A vector of eigenvalues for each quantum number block.
   */
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }

  /**
   * @brief Returns the fermion sign for each basis state.
   *
   * @return A vector of fermion signs.
   */
  [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }

  /**
   * @brief \f$ f^\dagger\f$ operator
   */
  std::vector<qOperator> f_dag_operator;

  /**
   * @brief Eigenvalues of the Hamiltonian
   */
  std::vector<std::vector<double>> eigenvalues_Q;

  /**
   * @brief Fermion sign of the each basis states. \f$ (-1)^ {n_{particle}}\f$
   */
  std::vector<double> chi_Q;

  /**
   * @brief Quantum numbers of the basis states.
   */
  std::vector<std::vector<int>> n_Q;

  std::vector<qOperator>        doubleOccupancy;

  /**
   * @brief Returns the double occupancy operator.
   *
   * This operator is useful for static or dynamic calculations.
   *
   * @return A vector of `qOperator` objects representing double occupancy.
   */
  [[nodiscard]] auto getDoubleOcc() const { return doubleOccupancy; }
};
