#pragma once
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

/**
 * @class twoChannel
 * @brief Represents a two-channel spin-1/2 impurity model.
 *
 * This class models a quantum impurity with two independent conduction channels,
 * each with spin-up and spin-down electrons. It can be used to construct
 * multi-channel Single Impurity Anderson Models (SIAM).
 *
 * The model conserves total particle number and total spin z-component across both channels.
 *
 * @see spinhalf - The single-channel version of this model.
 * @see twoChannelSiam - An example of a SIAM model built with this class.
 */
class twoChannel {
public:
  /**
   * @brief Constructs a twoChannel object.
   *
   * Initializes the basis and operators for a two-channel spin-1/2 impurity.
   * The Hamiltonian for the impurity itself is zero, as interactions are typically
   * defined when coupling to a bath in a larger model.
   */
  twoChannel();

  /**
   * @brief Get the basis quantum numbers for all eigenstates.
   *
   * Returns all valid (N_ch1, S_z_ch1, N_ch2, S_z_ch2) quantum number combinations.
   *
   * @return Vector of quantum number vectors.
   */
  std::vector<std::vector<int>> get_basis();

  /**
   * @brief Get the ground state energy for each quantum number sector.
   *
   * For the bare impurity, all eigenvalues are zero.
   *
   * @return Vector of eigenvalue vectors.
   */
  std::vector<std::vector<double>> get_eigenvaluesQ();

  /**
   * @brief Get the fermion parity factor for each quantum number sector.
   *
   * Returns (-1)^N for each sector, where N is the total particle number.
   *
   * @return Vector of fermion signs: ±1 for each quantum number block.
   */
  std::vector<double> get_chi_Q();

  /**
   * @brief Fermion creation operators f† in the quantum-number block basis.
   *
   * A vector of qOperator objects for each spin and channel combination.
   */
  std::vector<qOperator> f_dag_operator;

  /**
   * @brief Eigenvalues of the Hamiltonian in each quantum number sector.
   *
   * For the bare impurity, all eigenvalues are zero.
   */
  std::vector<std::vector<double>> eigenvalues_Q;

  /**
   * @brief Fermion parity sign for each basis quantum number sector.
   */
  std::vector<double> chi_Q;

  /**
   * @brief Quantum numbers labeling each Hamiltonian block.
   */
  std::vector<std::vector<int>> n_Q;
};
