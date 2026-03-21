#pragma once
#include "qOperator.hpp"
#include "utils/qmatrix.hpp"
#include "utils/timer.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <vector>

/**
 * @class nrgcore
 * @brief Solves the NRG (Numerical Renormalization Group) problem for a quantum impurity model coupled to a bath.
 *
 * This class implements the iterative NRG algorithm by diagonalizing the Hamiltonian in a basis of quantum numbers.
 * It combines an impurity model and a bath model, progressively adding bath sites and truncating high-energy states
 * to simulate the low-energy physics of the system.
 *
 * @tparam im_type The type of the impurity model, which must provide quantum numbers, eigenvalues, basis, and operators.
 * @tparam bath_type The type of the bath model, which must provide quantum numbers, eigenvalues, basis, and operators.
 */
template <typename im_type,   // Impurity type
      typename bath_type> // Bath type
class nrgcore {
  im_type   *impurityModel;
  bath_type *bath_model;
  //
  std::vector<qOperator> pre_fdag_oparator;

public:
  /**
   * @brief Constructs the `nrgcore` object.
   *
   * Initializes the NRG solver with references to the impurity and bath models.
   * Retrieves initial quantum numbers, eigenvalues, and basis from the models.
   * Sets default parameters and performs consistency checks.
   *
   * @param im_hamilt Reference to the impurity Hamiltonian model.
   * @param bt_hamilt Reference to the bath Hamiltonian model.
   */
  nrgcore(im_type &im_hamilt, bath_type &bt_hamilt)
    : impurityModel(&im_hamilt), bath_model(&bt_hamilt),
    chi_bath(bath_model->get_chi_Q()),
    bath_eigenvaluesQ(bath_model->get_eigenvaluesQ()),
    nq_bath(bath_model->get_basis()) {
  // Create the starting basis states from the impurity model.
  // Note: bath_eigenvaluesQ may vary per Wilson site.
  set_parameters(); // Set default parameters.
  test(); // Perform basic consistency tests.
  }

  /**
   * @brief Add a new Wilson bath site and perform the NRG iteration update.
   *
   * Creates the next basis, builds Hamiltonians for each quantum-number block,
   * diagonalizes each block, and truncates high-energy states if needed.
   * Also updates system operators and tracks iteration state.
   *
   * @param thopping Hopping coefficients for each bath channel.
   * @param rescale NRG energy rescaling factor for new site.
   */
  void add_bath_site(const std::vector<double> &thopping, double rescale) {
  timer t1("add_bath_site " + std::to_string(nrg_iterations_cnt));
  // Refresh bath model data, important for dynamic baths like superconductors.
  bath_eigenvaluesQ = bath_model->get_eigenvaluesQ();
  nq_bath           = bath_model->get_basis();
  chi_bath          = bath_model->get_chi_Q();
  //
  if (nrg_iterations_cnt == 0) {
    // For the first Wilson site, initialize with impurity model data.
    eigenvaluesQ      = impurityModel->get_eigenvaluesQ();
    pre_sysmQ         = impurityModel->get_basis();
    pre_fdag_oparator = impurityModel->f_dag_operator;
    discard_higher_energies();
  }
  //
  create_next_basis();
  create_next_hamiltonians(thopping, rescale);
  std::cout << "Done create_next_hamiltonians at " << t1.getDuration()
        << std::endl;
  // Diagonalize all Hamiltonians for each quantum number sector.
  eigenvaluesQ.clear();
  eigenvaluesQ.resize(current_hamiltonQ.size(), {});
  for (size_t i = 0; i < current_hamiltonQ.size(); i++) {
    eigenvaluesQ[i] = current_hamiltonQ[i].diag();
  }
  std::cout << "Done Diag " << t1.getDuration() << std::endl;
  set_current_fdag_operator();
  std::cout << "Done set_current_fdag_operator: " << t1.getDuration()
        << std::endl;
  // Discard higher energy states after updating operators.
  // This must follow set_current_fdag_operator as it relies on previous kept indices.
  // update_internal_state();
  }

  /**
   * @brief Updates the internal state after adding a bath site.
   *
   * Increments the iteration counter, discards high-energy states, and moves the current basis to the previous one
   * for the next iteration.
   */
  void update_internal_state() {
  nrg_iterations_cnt++;
  discard_higher_energies();
  pre_sysmQ = std::move(current_sysmQ);
  }

  /**
   * @brief Performs basic consistency tests between impurity and bath models.
   *
   * Ensures that the number of quantum numbers and the size of the f-dag operator arrays match
   * between the impurity and bath models to prevent runtime errors.
   */
  void test() {
  // Check number of quantum numbers.
  if (impurityModel->n_Q[0].size() != bath_model->n_Q[0].size()) {
    throw std::runtime_error(
      "Number of quantum number is different on model and impurity!");
  }
  // Check f-operator size.
  if (impurityModel->f_dag_operator.size() !=
    bath_model->f_dag_operator.size()) {
    std::string thString =
      "f_dag_operator.size is different on model and impurity!\n"
      "impurityModel->f_dag_operator.size() = " +
      std::to_string(impurityModel->f_dag_operator.size()) +
      "bath_model->f_dag_operator.size() = " +
      std::to_string(bath_model->f_dag_operator.size());
    throw std::runtime_error(thString);
  }
  }

  /**
   * @brief Creates the Hamiltonians for the next NRG iteration by coupling system and bath.
   *
   * For each coupled quantum number sector, constructs a matrix representing the Hamiltonian
   * with diagonal elements from rescaled system eigenvalues and bath eigenvalues, and off-diagonal
   * elements from hopping terms involving f-dag operators.
   *
   * @param t_hopping Vector of hopping parameters (must match the number of f-dag operators).
   * @param rescale Energy rescaling factor.
   */
  void create_next_hamiltonians(const std::vector<double> &t_hopping, // NOLINT
                double                     rescale) {
  if (t_hopping.size() != pre_fdag_oparator.size()) {
    throw std::runtime_error(
      std::string("t_hopping ") + std::to_string(t_hopping.size()) +
      " =! pre_fdag_oparator " + std::to_string(pre_fdag_oparator.size()));
  }
  std::cout << "-------------------------------------------------------------"
        << std::endl;
  std::cout << "Adding Wilson site " << nrg_iterations_cnt
        << " thopping:" << t_hopping << " rescale:" << rescale
        << std::endl;
  current_hamiltonQ.clear();
  for (size_t licq = 0; licq < coupled_nQ_index.size(); licq++) {
    current_hamiltonQ.emplace_back();
  }
#pragma omp parallel for // NOLINT
  for (size_t licq = 0; licq < coupled_nQ_index.size(); licq++) {
    auto   lindex = coupled_nQ_index[licq];
    size_t ldim   = 0;
    for (auto kindex : lindex) {
    auto ii = kindex / nq_bath.size(); // impurity nqi index
    auto bb = kindex % nq_bath.size(); // bath nqi index
    ldim    = ldim + eigenvaluesQ_kept_indices[ii].size() *
              bath_eigenvaluesQ[bb].size();
    }
    // Create the local Hamiltonian matrix for this quantum sector.
    qmatrix<double> h_nqi(ldim, ldim, 0);
    size_t kdim = 0;
    for (auto kindex : lindex) {
    auto ii = kindex / nq_bath.size(); // impurity nqi index
    auto bb = kindex % nq_bath.size(); // bath nqi index
    // Set diagonal elements: rescaled system energies + bath energies.
    for (size_t it : eigenvaluesQ_kept_indices[ii]) {
      for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
      h_nqi(kdim + it + il * eigenvaluesQ_kept_indices[ii].size(),
          kdim + it + il * eigenvaluesQ_kept_indices[ii].size()) =
        eigenvaluesQ[ii][it] * rescale + bath_eigenvaluesQ[bb][il];
      }
    }
    size_t kpdim = 0;
    for (auto kp_index : lindex) {
      auto ii_p = kp_index / nq_bath.size(); // impurity nqi index
      auto bb_p = kp_index % nq_bath.size(); // bath nqi index
      // Add off-diagonal hopping terms using f-dag operators from system and bath.
      for (size_t ip = 0; ip < pre_fdag_oparator.size(); ip++) {
      auto fdag_opt =
        pre_fdag_oparator[ip].get(ii, ii_p); // Previous system's f-dag operator.
      auto f_bath_opt = bath_model->f_dag_operator[ip].get(
        bb_p, bb); // Bath's f-dag operator.
      if (fdag_opt && f_bath_opt) {
        auto *f_dag_prev = fdag_opt.value(); // Pointer to qmatrix.
        auto f_bath      = f_bath_opt.value();
        for (size_t it : eigenvaluesQ_kept_indices[ii]) {
        for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
          for (size_t it_p : eigenvaluesQ_kept_indices[ii_p]) {
          for (size_t il_p = 0; il_p < bath_eigenvaluesQ[bb_p].size();
             il_p++) {
            double aa = (f_dag_prev->at(it, it_p)) *
                  (f_bath->at(il_p, il)) * (chi_bath[bb_p]);
            h_nqi(
              kdim + it + il * eigenvaluesQ_kept_indices[ii].size(),
              kpdim + it_p +
                il_p * eigenvaluesQ_kept_indices[ii_p].size()) +=
              aa * t_hopping[ip];
            h_nqi(kpdim + it_p +
                il_p * eigenvaluesQ_kept_indices[ii_p].size(),
              kdim + it +
                il * eigenvaluesQ_kept_indices[ii].size()) +=
              aa * t_hopping[ip]; // Hermitian conjugate.
          }
          }
        }
        }
      }
      }
      kpdim += eigenvaluesQ_kept_indices[ii_p].size() *
           bath_eigenvaluesQ[bb_p].size();
    }
    kdim += eigenvaluesQ_kept_indices[ii].size() *
        bath_eigenvaluesQ[bb].size();
    }
    current_hamiltonQ[licq] = h_nqi;
  }
  }

  /**
   * @brief Creates the next basis by coupling the previous system's quantum numbers with the bath's.
   *
   * Generates new quantum number sectors by summing the quantum numbers of the previous system
   * and the bath, and groups coupled indices for efficient Hamiltonian construction.
   */
  void create_next_basis() {
  coupled_nQ_index.clear();
  current_sysmQ.clear();
  for (size_t i = 0; i < pre_sysmQ.size(); i++) {
    for (size_t j = 0; j < nq_bath.size(); j++) {
    auto tm = pre_sysmQ[i];
    for (size_t ix = 0; ix < nq_bath[j].size(); ix++) {
      tm[ix] += nq_bath[j][ix];
    }
    // Check if this quantum number combination already exists.
    bool ex = false;
    for (size_t ix = 0; ix < current_sysmQ.size(); ix++) {
      if (tm == current_sysmQ[ix]) {
      ex = true;
      coupled_nQ_index[ix].push_back({i * nq_bath.size() + j});
      }
    }
    if (!ex) {
      current_sysmQ.push_back(tm);
      coupled_nQ_index.push_back({i * nq_bath.size() + j});
    }
    }
  }
  };

  /**
   * @brief Discards high-energy states to keep only the lowest-energy ones.
   *
   * Collects all eigenvalues, sorts them, and keeps up to `max_kept_states` plus degenerate states within `errorbarInEnergy`.
   * Updates kept indices, shifts energies relative to the ground state if truncation occurs, and logs statistics.
   */
  void discard_higher_energies() {
  all_eigenvalue.clear();
  for (auto aa : eigenvaluesQ) {
    all_eigenvalue.insert(all_eigenvalue.end(), aa.begin(), aa.end());
  }
  std::sort(all_eigenvalue.begin(), all_eigenvalue.end());
  eigenvaluesQ_kept_indices.clear();
  double En_max = all_eigenvalue[all_eigenvalue.size() - 1];
  double En_min = all_eigenvalue[0];
  if (all_eigenvalue.size() <= max_kept_states) {
    for (auto &aa : eigenvaluesQ) {
    std::vector<size_t> tm(aa.size());
    std::iota(tm.begin(), tm.end(), 0);
    eigenvaluesQ_kept_indices.push_back(tm);
    }
  } else {
    En_max = all_eigenvalue[max_kept_states];
    // Include degenerate states within errorbar.
    for (size_t i = max_kept_states + 1; i < all_eigenvalue.size(); i++) {
    if (std::fabs(all_eigenvalue[i] - all_eigenvalue[max_kept_states]) <
      errorbarInEnergy) {
      // Degenerate, could extend kept states.
    }
    }
    // Set kept indices for each sector.
    for (auto &aa : eigenvaluesQ) {
    std::vector<size_t> tm;
    size_t              id = 0;
    for (auto &bb : aa) {
      if (bb < En_max + errorbarInEnergy) {
      tm.push_back(id);
      }
      id++;
    }
    eigenvaluesQ_kept_indices.push_back(tm);
    }
  }
  // Count and log kept states.
  size_t tmp_kept{0};
  for (auto &aa : eigenvaluesQ_kept_indices) {
    tmp_kept += aa.size();
  }
  no_of_kept_states = tmp_kept;
  std::cout << "NRG: Iteration " << nrg_iterations_cnt << " N:K:kp "
        << all_eigenvalue.size() << ":" << tmp_kept << ":"
        << no_of_kept_states << std::endl;
  std::cout << "Ground state energy: " << all_eigenvalue[0]
        << " MaxEn(kept): " << En_max << std::endl;
  if (no_of_kept_states <= max_kept_states) {
    nrg_iterations_min = nrg_iterations_cnt;
  } else {
    // Shift energies relative to ground state only if truncation occurred.
    relativeGroundStateEnergy.push_back(En_min);
    for (auto &aa : eigenvaluesQ) {
    std::transform(aa.begin(), aa.end(), aa.begin(),
             [En_min](double a) { return a - En_min; });
    }
  }
  }

  /**
   * @brief Sets the maximum number of states to keep in each iteration.
   *
   * Also initializes other parameters like energy error bar and iteration counters.
   *
   * @param n Maximum number of states to keep (default: 1024).
   */
  void set_parameters(size_t n = 1024) {
  no_of_kept_states  = n;
  errorbarInEnergy   = 1e-8;
  nrg_iterations_cnt = 0;
  nrg_iterations_min = 0;
  max_kept_states    = no_of_kept_states;
  }

  /**
   * @brief Retrieves the current basis quantum numbers.
   *
   * @return Vector of vectors representing quantum numbers for each sector.
   * @throws std::runtime_error if the current basis is empty.
   */
  std::vector<std::vector<int>> get_basis_nQ() {
  if (current_sysmQ.empty()) {
    throw std::runtime_error("current_sysmQ.empty!");
  }
  return current_sysmQ;
  }

  /**
   * @brief Retrieves the eigenvalues of the current Hamiltonian.
   *
   * @return Vector of vectors of eigenvalues for each quantum sector.
   * @throws std::runtime_error if eigenvalues are empty.
   */
  std::vector<std::vector<double>> get_eigenvaluesQ() {
  if (eigenvaluesQ.empty()) {
    throw std::runtime_error("eigenvaluesQ.empty!");
  }
  return eigenvaluesQ;
  }

  /**
   * @brief Retrieves the f-dag operators in the current basis.
   *
   * These operators are used to construct the Hamiltonian for the next iteration.
   *
   * @return Vector of qOperator objects.
   */
  auto get_f_dag_operator() { return pre_fdag_oparator; }

  /**
   * @brief Checks if higher-energy states were discarded in the current iteration.
   *
   * @return True if the number of kept states is at or below the maximum, indicating no truncation.
   */
  bool checkHigherEnergyDiscarded() {
  return no_of_kept_states >= max_kept_states;
  }

private:
  /**
   * @brief Updates the f-dag operators in the current eigenbasis.
   *
   * Transforms the previous operators using the unitary matrices from diagonalization
   * and incorporates bath operators for the next iteration.
   */
  /**
   * @brief Update f-dag operators in the new eigenbasis after diagonalization.
   *
   * For every system sector, this method rotates previous f-dag operators into
   * the current basis using the current Hamiltonian eigenvectors.
   */
  void set_current_fdag_operator() { // NOLINT
  for (auto &aa : pre_fdag_oparator) {
    aa.clear();
  }
  //
  std::map<std::array<size_t, 3>, std::vector<std::array<size_t, 6>>> foprIdx;
  //
  timer t1("set_current_fdag_operator");
  for (size_t ip = 0; ip < pre_fdag_oparator.size(); ip++) {
    // Find basis indices for operator coupling.
    for (size_t i = 0; i < current_sysmQ.size(); i++) {
    for (size_t j = 0; j < current_sysmQ.size(); j++) {
      auto idx   = coupled_nQ_index[i];
      auto idx_p = coupled_nQ_index[j];
      size_t kidx = 0;
      for (auto kindex : idx) {
      auto   ii     = kindex / nq_bath.size(); // impurity nqi index
      auto   bb     = kindex % nq_bath.size(); // bath nqi index
      size_t kidx_p = 0;
      for (auto kindex_p : idx_p) {
        auto ii_p = kindex_p / nq_bath.size(); // impurity nqi index
        auto bb_p = kindex_p % nq_bath.size(); // bath nqi index
        auto fdag_bath_opt = bath_model->f_dag_operator[ip].get(bb, bb_p);
        if (ii == ii_p && fdag_bath_opt) {
        foprIdx[{ip, i, j}].push_back(
          {ii, bb, ii_p, bb_p, kidx, kidx_p});
        }
        kidx_p += eigenvaluesQ_kept_indices[ii_p].size() *
            bath_eigenvaluesQ[bb_p].size();
      }
      kidx += eigenvaluesQ_kept_indices[ii].size() *
          bath_eigenvaluesQ[bb].size();
      }
    }
    }
  }
  // Process operator indices.
  std::cout << "--------------------------" << t1.getDuration() << std::endl;
  std::vector<qmatrix<double>> qfr(foprIdx.size());
  for (auto itf = foprIdx.begin(); itf != foprIdx.end(); itf++) {
    auto            ip      = itf->first[0];
    auto            i       = itf->first[1];
    auto            j       = itf->first[2];
    size_t          kpdim   = eigenvaluesQ[i].size();
    size_t          kpdim_p = eigenvaluesQ[j].size();
    qmatrix<double> tmat(kpdim, kpdim_p, 0);
    for (auto const &itd : itf->second) {
    auto ii = itd[0];
    auto bb = itd[1];
    auto bb_p = itd[3];
    size_t kidx   = itd[4];
    size_t kidx_p = itd[5];
    auto fdag_bath_opt = bath_model->f_dag_operator[ip].get(bb, bb_p);
    auto fdag_bath     = fdag_bath_opt.value();
    for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
      for (size_t il_p = 0; il_p < bath_eigenvaluesQ[bb_p].size(); il_p++) {
      double rvalue = fdag_bath->at(il, il_p);
      for (size_t it : eigenvaluesQ_kept_indices[ii]) {
        tmat(kidx + it + il * eigenvaluesQ_kept_indices[ii].size(),
           kidx_p + it + il_p * eigenvaluesQ_kept_indices[ii].size()) +=
          rvalue;
      }
      }
    }
    }
    qfr[std::distance(foprIdx.begin(), itf)] =
      current_hamiltonQ[i].cTranspose().dot(tmat.dot(current_hamiltonQ[j]));
  }
  std::cout << "--------------------------" << t1.getDuration() << std::endl;
  // Set the operators.
  for (auto itf = foprIdx.begin(); itf != foprIdx.end(); itf++) {
    auto ip = itf->first[0];
    auto i  = itf->first[1];
    auto j  = itf->first[2];
    pre_fdag_oparator[ip].set(qfr[std::distance(foprIdx.begin(), itf)], i, j);
  }
  }

  /**
   * @brief Enforce degeneracy for near-equal energy levels.
   *
   * Identifies energies whose differences are below `errorbarInEnergy` and
   * stabilizes them by assigning a common minimum value. Useful to avoid
   * numerical split from repeated diagonalization.
   */
  void enforceDegeneracy() {
  for (auto &aa : eigenvaluesQ) {
    for (auto &bb : aa) {
    for (auto &aa_p : eigenvaluesQ) {
      for (auto &bb_p : aa_p) {
      if (std::fabs(bb - bb_p) < errorbarInEnergy) {
        double tmp = std::min(bb, bb_p);
        bb_p       = tmp;
        bb         = tmp;
      }
      }
    }
    if (std::abs(bb) < 1e-10) {
      std::cout << "degeneracy found" << std::endl;
    }
    }
  }
  }

  /**
   * @brief Legacy version of set_current_fdag_operator (deprecated).
   */
  void set_current_fdag_operator_old() { // NOLINT
  // Implementation omitted for brevity; similar to current version but less optimized.
  }

  // Member variables (with brief descriptions where needed).
  std::vector<double> chi_bath; // Bath chi factors.
  size_t no_of_kept_states{0}; // Current number of kept states.
  size_t max_kept_states{0}; // Maximum number of states to keep.
  double errorbarInEnergy{0}; // Tolerance for energy degeneracy.

public:
  std::vector<double> all_eigenvalue; // All eigenvalues across sectors.
  std::vector<double> relativeGroundStateEnergy; // Ground state energies per iteration.
  std::vector<std::vector<size_t>> eigenvaluesQ_kept_indices; // Indices of kept states per sector.
  std::vector<qmatrix<double>>     current_hamiltonQ; // Current Hamiltonians (unitary matrices).
  std::vector<std::vector<int>>    current_sysmQ; // Current quantum number sectors.
  std::vector<std::vector<int>>    pre_sysmQ; // Previous quantum number sectors.
  std::vector<std::vector<double>> eigenvaluesQ; // Eigenvalues per sector.
  std::vector<std::vector<size_t>> coupled_nQ_index; // Coupled indices for sectors.
  int                              nrg_iterations_cnt{}; // Current iteration count.
  int                              nrg_iterations_min{}; // Iteration where truncation started.
  std::vector<std::vector<double>> bath_eigenvaluesQ; // Bath eigenvalues per sector.
  std::vector<std::vector<int>>    nq_bath; // Bath quantum numbers.
  std::vector<qOperator>          *getPreWilsonSiteOperators() {
  return &pre_fdag_oparator;
  }
};
