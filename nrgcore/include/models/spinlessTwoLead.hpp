#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"
/**
 * @class spinnlessTwoLead
 * @brief Represents a spinless two-lead model.
 *
 * This class models a system with one spinless orbital coupled to two leads.
 * It is a simplified model often used to study transport properties in quantum dots.
 * The model only conserves charge.
 *
 * The constructor builds the fermion basis, sets up operators, and
 * diagonalizes the local Hamiltonian for each charge sector.
 */
class spinnlessTwoLead {
public:
  /**
   * @brief Constructs a spinnlessTwoLead object.
   *
   * Initializes the basis and operators for a spinless two-lead model.
   *
   * @param teps The onsite energy of the orbital.
   */
  explicit spinnlessTwoLead(double teps = 0) {
    fermionBasis spinnlessTwoLeadBasis(
        2,                         //  two for the lead
        fermionBasis::chargeOnly); // Number of fermion channels/spins
    f_dag_operator = spinnlessTwoLeadBasis.get_f_dag_operator();
    auto f_dag_raw = spinnlessTwoLeadBasis.get_raw_f_dag_operator();
    // n_Q
    n_Q = spinnlessTwoLeadBasis.get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    // std::cout << "chi_Q: " << chi_Q << std::endl;
    // set up the Hamiltonian
    auto n_label = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
    // std::cout << "n_up" << n_up << "n_down" << n_down;
    auto H = n_label * teps // onsite energy
        ;                   // Columb Energy
    // get the hamiltonian  in the blocked structure
    auto h_blocked = spinnlessTwoLeadBasis.get_block_Hamiltonian(H);
    //  h_blocked.display();
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
    }
    // TODO(sp): rotate the f operator
    // End of the constructor
    // left and right operators are only needed
    // TODO(sp): rotate the f operator  and n operator
  }
  /**
   * @brief Get the basis quantum numbers for all eigenstates.
   *
   * Returns the sector basis `n_Q` as a list of charge vectors.
   * For spinnlessTwoLead, each inner vector corresponds to a unique
   * charge configuration used by the fermionBasis blocking.
   *
   * @return A vector of vectors, where each inner vector represents the quantum numbers of a state.
   */
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const {
    return n_Q;
  }

  /**
   * @brief Get the ground state energy for each quantum number sector.
   *
   * Returns eigenvalues of the local Hamiltonian for each `n_Q` sector.
   * The outer vector index matches the Q sector from `get_basis()`.
   *
   * @return A vector of vectors, where each inner vector contains the eigenvalues for a quantum number sector.
   */
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }
  /**
   * @brief Get the fermion parity factor for each quantum number sector.
   *
   * Returns (-1)^N for each sector, where N is the total particle number.
   *
   * @return A vector of doubles, where each element is the fermion parity sign for a quantum number sector.
   */
  [[nodiscard]] std::vector<double> get_chi_Q() const {
    return chi_Q;
  }
  /**
   * @brief Fermion creation operators f† in the quantum-number block basis.
   */
  std::vector<qOperator>           f_dag_operator;
  /**
   * @brief Eigenvalues of the Hamiltonian in each quantum number sector.
   */
  std::vector<std::vector<double>> eigenvalues_Q;
  /**
   * @brief Fermion parity sign for each basis quantum number sector.
   */
  std::vector<double>              chi_Q;
  /**
   * @brief Quantum numbers labeling each Hamiltonian block.
   */
  std::vector<std::vector<int>>    n_Q;
};
