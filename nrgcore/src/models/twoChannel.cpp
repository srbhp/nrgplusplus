#include "models/twoChannel.hpp"
#include "models/fermionBasis.hpp"
#include <cmath>
#include <cstddef>
twoChannel::twoChannel() {
  //
  // create_basis();
  // set_foperator();
  // set_chi_Q();
  // create_hamiltonian();
  fermionBasis tChlBasis(
      4, fermionBasis::chargeAndSpin); // Number of fermion channels/spins
  f_dag_operator = tChlBasis.get_f_dag_operator();
  n_Q            = tChlBasis.get_unique_Qnumbers();
  // set f_operator
  auto f_dag_raw = tChlBasis.get_raw_f_dag_operator();
  auto f_raw     = f_dag_raw;
  for (auto &aa : f_raw) {
    aa = aa.cTranspose();
  }
  // Set Chi_Q
  chi_Q.clear();
  for (auto ai : n_Q) {
    double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
    chi_Q.push_back(std::pow(-1., t_charge));
  }
  // set up the Hamiltonian
  auto n_up   = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
  auto n_down = f_dag_raw[1].dot(f_dag_raw[1].cTranspose());
  // std::cout << "n_up" << n_up << "n_down" << n_down;
  auto H         = (n_up + n_down) * 0.0;
  auto h_blocked = tChlBasis.get_block_Hamiltonian(H);
  //  h_blocked.display();
  // Diagonalize the hamilton
  eigenvalues_Q.clear();
  eigenvalues_Q.resize(n_Q.size(), {});
  for (size_t i = 0; i < n_Q.size(); i++) {
    eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
  }
  // TODO: rotate the f operator
}
std::vector<std::vector<int>> twoChannel::get_basis() {
  /** returns the basis vector
   *
   */
  return n_Q;
}
std::vector<std::vector<double>> twoChannel::get_eigenvaluesQ() {
  /** returns eigenvalues_Q
   *
   *
   */
  return eigenvalues_Q;
}
std::vector<double> twoChannel::get_chi_Q() {
  /** This functions returns
   *  `vector<vector>` of
   *  \f$ \chi_Q  = e^{n_Q} \f$
   *
   *
   */
  return chi_Q;
}
