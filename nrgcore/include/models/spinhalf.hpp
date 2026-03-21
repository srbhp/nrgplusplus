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
 * @brief Single orbital with spin-up and spin-down electrons (spin-1/2 impurity).
 *
 * Represents a single-site impurity with both spin-up and spin-down orbitals.
 * This is the fundamental quantum dot model used to construct NRG problems
 * like the Single Impurity Anderson Model (SIAM) and Kondo model.
 *
 * The Hamiltonian includes:
 * - Onsite energy (ε_d) for the orbital
 * - Coulomb interaction energy (U) between spin-up and spin-down
 * - Magnetic field coupling (h) for Zeeman splitting
 *
 * Quantum numbers:
 * - Conserved: particle number and spin z-component (N, S_z)
 * - Symmetry basis block-diagonalizes the space
 *
 * @example
 * @code
 * double eps = -1.0;        // Onsite energy in eV
 * double U = 2.0;           // Coulomb energy in eV
 * double B = 0.1;           // Magnetic field in eV (optional)
 * spinhalf impurity(eps, U, B);
 * @endcode
 *
 * @see fermionBasis - Base class using charge/spin conservation
 * @see SIAM - Single Impurity Anderson Model built with this impurity
 */
class spinhalf {
  /// @brief Creation operators in the single-particle basis before block-diagonalization
  std::vector<qmatrix<>> f_dag_raw;
  /// @brief Coulomb (Hubbard) interaction strength: c†↑c↑ c†↓c↓
  double                 Uint;
  /// @brief Onsite orbital energy: single-particle level
  double                 epsilon_d;
  /// @brief Magnetic field for Zeeman splitting: (n↑ - n↓) * h
  double                 magnetic_field;

public:
  /**
   * @brief Construct spin-1/2 impurity with specified parameters.
   *
   * Sets up the two-site fermion basis (spin-up, spin-down) and constructs
   * the single-impurity Hamiltonian with charge and spin conservation symmetries.
   *
   * The Hamiltonian is: H = εd(n↑ + n↓) + h(n↑ - n↓) + U·n↑·n↓
   * where:
   * - n↑, n↓ are occupation numbers for spin-up and spin-down
   * - U is Hubbard interaction (Coulomb repulsion)
   * - h is magnetic field strength
   *
   * @param teps Onsite orbital energy in eV
   * @param tUint Coulomb interaction parameter U in eV (charge repulsion)
   * @param tmag Magnetic field coupling h in eV (Zeeman splitting, default: 0)
   *
   * @note Automatically diagonalizes Hamiltonian in the (n_Q, s_z) basis
   * @note Conserves: particle number (0-2 electrons) and spin projection
   * @note Member variables populated: n_Q, chi_Q, eigenvalues_Q, f_dag_operator
   *
   * @see get_eigenvaluesQ() - Ground state energies for each quantum number sector
   * @see get_unique_Qnumbers() - Valid (N, S_z) quantum number combinations
   *
   * @example
   * @code
   * // SIAM parameters: single level at -1eV with U=2eV
   * spinhalf quantum_dot(-1.0, 2.0, 0.0);
   * @endcode
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
   * @brief Get the basis quantum numbers for all eigenstates.
   *
   * Returns all valid (N, S_z) quantum number pairs where:
   * - N is total particle number (0, 1, or 2)
   * - S_z is z-component of spin (-1, 0, +1 in half-unit spacing)
   *
   * @return Vector of quantum number vectors: each element is [N, S_z]
   * @note Ordered by construction; corresponds to Hamiltonian block structure
   * @see get_eigenvaluesQ() - energies for each quantum number sector
   */
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; };

  /**
   * @brief Get the ground state energy for each quantum number sector.
   *
   * Returns lowest eigenvalue of the spin-1/2 Hamiltonian in each (N, S_z) sector.
   * Multiple energies per sector indicate degeneracy.
   *
   * @return Vector of eigenvalue vectors: eigenvalues_Q[i] = all eigenvalues in sector i
   * @note Eigenvalues within each sector are sorted in ascending order
   * @see get_basis() - quantum numbers corresponding to each sector
   */
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }

  /**
   * @brief Get the fermion parity factor for each quantum number sector.
   *
   * Fermionic operators anticommute, resulting in sign factors ±1 when ordering.
   * Returns (-1)^N for each sector, where N is particle number.
   *
   * @return Vector of fermion signs: ±1 for each quantum number block
   * @note Used for correct anticommutation relations in fermionic algebra
   */
  [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }

  /**
   * @brief Get fermion creation operators f† in quantum-number block basis.
   *
   * Returns f† operators block-diagonalized by (N, S_z) quantum numbers.
   * Useful for computing spectral functions and Green's functions.
   *
   * @return Vector of qOperator: f†_spin for spin ∈ {up, down} channels
   */
  std::vector<qOperator> f_dag_operator;

  /**
   * @brief Eigenvalues of the Hamiltonian in each quantum number sector.
   *
   * eigenvalues_Q[i] contains all eigenvalues for the i-th quantum number block.
   * Used for FDM calculations and spectral density computations.
   */
  std::vector<std::vector<double>> eigenvalues_Q;

  /**
   * @brief Fermion parity sign for each basis quantum number sector.
   *
   * Returns (-1)^N where N is the particle number in that sector.
   * This sign is essential for fermionic anticommutation relations.
   */
  std::vector<double> chi_Q;

  /**
   * @brief Quantum numbers (N, S_z) labeling each Hamiltonian block.
   *
   * Used to organize eigenstates and operators by conserved quantum numbers.
   * Allows block-diagonal structure exploited for computational efficiency.
   */
  std::vector<std::vector<int>> n_Q;

  /// @brief Double occupancy n↑·n↓ operator in quantum number basis
  std::vector<qOperator>        doubleOccupancy;

  /**
   * @brief Get the double occupancy operator for computing correlations.
   *
   * The double occupancy operator D = n↑·n↓ projects onto states with
   * both spin-up and spin-down electrons present simultaneously.
   * Useful for:
   * - Computing correlation functions and susceptibilities
   * - Static properties like charge susceptibility
   * - Comparison with experimental charge compressibility
   *
   * @return Vector of qOperator objects representing D in quantum number basis
   * @note Only has weight in sectors with N=2 (doubly-occupied state)
   */
  [[nodiscard]] auto getDoubleOcc() const { return doubleOccupancy; }
};
