#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"
class spinnlessTwoLead {
  /** This class is for a single orbital with
   * coupled to Two lead. Only chrge is conserved
   *
   *
   */
public:
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
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const {
    /** returns the basis vector
     *
     */
    return n_Q;
  }
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    /** returns eigenvalues_Q
     *
     *
     */
    return eigenvalues_Q;
  }
  [[nodiscard]] std::vector<double> get_chi_Q() const {
    /** This functions returns
     *  `vector<vector>` of
     *  \f$ \chi_Q  = e^{n_Q} \f$
     *
     *
     */
    return chi_Q;
  }
  // protected:
  // parameter
  // functions
  //
  std::vector<qOperator>           f_dag_operator;
  std::vector<std::vector<double>> eigenvalues_Q;
  std::vector<double>              chi_Q;
  std::vector<std::vector<int>>    n_Q;
};
