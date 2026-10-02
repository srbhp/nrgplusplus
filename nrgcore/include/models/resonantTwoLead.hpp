#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"

/**
 * @brief Resonant two-lead impurity model with charge conservation.
 *
 * The model describes a single impurity site coupled to left and right leads,
 * conserving total charge while allowing tunneling and local interaction terms.
 */
class resonantTwoLead {
public:
  /**
   * @brief Construct the resonant two-lead impurity model.
   *
   * @param teps Onsite impurity energy.
   * @param tleftGamma Left-lead tunneling amplitude.
   * @param trightGamma Right-lead tunneling amplitude.
   * @param Uinterct Local interaction parameter between the impurity occupancy and lead occupancy.
   */
  resonantTwoLead(double teps, double tleftGamma, double trightGamma,
                  double Uinterct) {
    fermionBasis resonantTwoLeadBasis(
        3,                         // 1 for the label and two for the lead
        fermionBasis::chargeOnly); // Number of fermion channels/spins
    f_dag_operator = resonantTwoLeadBasis.get_f_dag_operator();
    auto f_dag_raw = resonantTwoLeadBasis.get_raw_f_dag_operator();
    // n_Q
    n_Q = resonantTwoLeadBasis.get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    // std::cout << "chi_Q: " << chi_Q << std::endl;
    // set up the Hamiltonian
    auto n_label = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
    auto n_left  = f_dag_raw[1].dot(f_dag_raw[1].cTranspose());
    auto n_right = f_dag_raw[2].dot(f_dag_raw[2].cTranspose());
    auto half_id = n_label.id();
    // std::cout << "n_up" << n_up << "n_down" << n_down;
    auto H = n_label * teps // onsite energy
             + (f_dag_raw[0].dot(f_dag_raw[1].cTranspose()) +
                f_dag_raw[1].dot(f_dag_raw[0].cTranspose())) *
                   tleftGamma // left lead coupling
             // Right lead coupling
             + (f_dag_raw[0].dot(f_dag_raw[2].cTranspose()) +
                f_dag_raw[2].dot(f_dag_raw[0].cTranspose())) *
                   trightGamma // left lead coupling
             + ((n_label - half_id).dot(n_left + n_right - half_id)) *
                   Uinterct; // Columb Energy
    // get the hamiltonian  in the blocked structure
    auto h_blocked = resonantTwoLeadBasis.get_block_Hamiltonian(H);
    //  h_blocked.display();
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
    }
    // TODO : rotate the f operator
    // End of the constructor
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
    f_dag_operator = {topr[1], topr[2]};
    impurityFdag   = topr[0];
    impurityNparticle.clear();
    for (size_t i = 0; i < n_Q.size(); i++) {
      for (size_t j = 0; j < n_Q.size(); j++) {
        auto tfopr = impurityFdag.get(i, j);
        if (tfopr) {
          impurityNparticle.set(tfopr.value()->dot(tfopr.value()->cTranspose()),
                                i, i);
        }
      }
    }
    // left and right operators are only needed
    // TODO: rotate the f operator  and n operator
  }
  /// @brief Creation operator for the impurity site in the rotated basis.
  qOperator                                   impurityFdag;
  /// @brief Number operator for the impurity site.
  qOperator                                   impurityNparticle;
  /// @brief Return the charge quantum-number sectors of the basis.
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const {
    /** returns the basis vector
     *
     */
    return n_Q;
  }
  /// @brief Return the Hamiltonian eigenvalues for each charge sector.
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }
  /// @brief Return the parity factor associated with each charge sector.
  [[nodiscard]] std::vector<double> get_chi_Q() const {
    return chi_Q;
  }
  /// @brief Creation operators for the impurity and lead modes.
  std::vector<qOperator>           f_dag_operator;
  /// @brief Eigenvalues in each block of the charge-conserving basis.
  std::vector<std::vector<double>> eigenvalues_Q;
  /// @brief Charge-parity sign for each sector.
  std::vector<double>              chi_Q;
  /// @brief Charge quantum numbers for each many-body state.
  std::vector<std::vector<int>>    n_Q;
};
