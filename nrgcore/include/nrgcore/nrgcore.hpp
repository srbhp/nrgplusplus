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
// Take two hamiltonian as a template parameter
// Create a new hamiltonian
/**
 * @brief We use this class  solve the NRG problem
 * of a bath and an Impurity. This class take two template parameter
 * one for the impurity and one for the bath. The bath and impurity
 * are the should have the same quantum number.
 *
 * The bath and impurity Hamiltonian is diagonalized in the block basis
 * of the quantum numbers. These eigenvales are zero for a metallic bath
 * and non-zero for a superconductor bath.
 *
 * We assume that the
 * impurity and bath are connected by a few interactions terms.
 *
 * \f[ H_{inter} = \sum_{i} (\lambda_i a_{i}^{\dagger}b_{i} + h.c ) \f]
 * where \f$ a_{i}^{\dagger} \f$ is the creation operator of the impurity
 * and \f$ b_{i}^{\dagger} \f$ is the creation operator of the bath.
 * The \f$ \lambda_i \f$ are the coupling constants.
 *
 *
 *
 *
 *
 * The number of operator for the bath and impurity should be the same.
 * If this is the case (i.e., Anderson Model) then we can pass the impurity
 * class and bath class to the nrgcore class as a  template parameter. If this
 * is not the case (i.e., Kondo Model )then  we create the impurity class out of
 * the impurity site and the first Wilson site. Then we pass this combined
 * impurity class and bath class to the nrgcore class.
 *
 *
 * @tparam im_type Type of Impurity class.
 * @tparam bath_type Type of bath class.
 * @param im_hamilt Impurity class  Hamiltonian.
 * @param bt_hamilt Bath class Hamiltonian.
 * @return [TODO:return]
 */
template <typename im_type,   // Impurity type
          typename bath_type> // bath type
class nrgcore {
  im_type   *impurityModel;
  bath_type *bath_model;
  //
  std::vector<qOperator> pre_fdag_oparator;

public:
  /**
   * @brief [TODO:description]
   *
   * @param im_hamilt [TODO:parameter]
   * @param bt_hamilt [TODO:parameter]
   */
  nrgcore(im_type &im_hamilt, bath_type &bt_hamilt)
      : impurityModel(&im_hamilt), bath_model(&bt_hamilt),
        chi_bath(bath_model->get_chi_Q()),
        bath_eigenvaluesQ(bath_model->get_eigenvaluesQ()),
        nq_bath(bath_model->get_basis()) {
    // create the staring basis states ....
    // TODO(sp): bath_eigenvaluesQ may be different on each wilson site
    set_parameters(); // set the default parameters
    test();
  }
  /**
   * @brief This function is called to add a bath site. This is
   * done for each iteration of the simulation. This function
   * create a full Hamiltonian from the impurity and bath class
   * and diagonalize it. The eigenvalues are stored in eigenvaluesQ.
   *
   * @param thopping This a array for the \f \lambda_i \f
   * @param rescale This is the rescaling factor  \f \sqrt{\Lambda} \f
   */
  void add_bath_site(const std::vector<double> &thopping, double rescale) {
    timer t1("add_bath_site " + std::to_string(nrg_iterations_cnt));
    // copy the bath model stuff. This important for superconductor bath
    // this was not if we had used pointer.
    bath_eigenvaluesQ = bath_model->get_eigenvaluesQ();
    nq_bath           = bath_model->get_basis();
    chi_bath          = bath_model->get_chi_Q();
    //
    if (nrg_iterations_cnt == 0) {
      // This is  only for the first wilson site
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
    // solve all the hamiltonians
    eigenvaluesQ.clear(); // TODO(sp): save values See if one needs
    eigenvaluesQ.resize(current_hamiltonQ.size(), {});
    for (size_t i = 0; i < current_hamiltonQ.size(); i++) {
      // std::cout << current_hamiltonQ[i].getcolumn() << " "
      //          << current_hamiltonQ[i].getrow() << std::endl
      // if (current_sysmQ[i] == std::vector{1, 1})
      //  std::cout << "Nqi" << current_sysmQ[i] << "current_hamiltonQ"
      //            << current_hamiltonQ[i];
      eigenvaluesQ[i] = current_hamiltonQ[i].diag();
      // std::cout << "eigenvaluesQ" << eigenvaluesQ[i];
    }
    std::cout << "Done Diag " << t1.getDuration() << std::endl;
    set_current_fdag_operator();
    std::cout << "Done set_current_fdag_operator: " << t1.getDuration()
              << std::endl;
    // discard higher energy states.
    // This should be done after set_current_fdag_operator
    // as previous eigenvaluesQ_kept_indices is needed
    // for set_current_fdag_operator
    // update_internal_state();
  }
  /**
   * @brief This function is called to discard the higher energy state and
   * update some  the internal state of the nrgcore class. This function
   * should be called after add_bath_site. If we need to update some
   * bath or impurity operators (i.e., f_dag_operator) then we should call
   * this function after the update is done.
   *
   */
  void update_internal_state() {
    nrg_iterations_cnt++;
    discard_higher_energies();
    pre_sysmQ = std::move(current_sysmQ);
  }
  /**
   * @brief [TODO: We have added few basic test for the impurity and bath class
   * here.]
   */
  void test() {
    // Number of quantum number
    if (impurityModel->n_Q[0].size() != bath_model->n_Q[0].size()) {
      throw std::runtime_error(
          "Number of quantum number is different on model and impurity!");
    }
    // foperator size
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
      // TODO(sp): Not all the eigenvalues are always available
      // create the local nqi hamiltonian
      qmatrix<double> h_nqi(ldim, ldim, 0);
      // std::cout << "-------------------------------------";
      // std::cout << "ldim " << ldim << " lindex " << lindex << std::endl;
      size_t kdim = 0;
      for (auto kindex : lindex) {
        auto ii = kindex / nq_bath.size(); // impurity nqi index
        auto bb = kindex % nq_bath.size(); // bath nqi index
        // set diagonal matrix elements
        for (size_t it : eigenvaluesQ_kept_indices[ii]) {
          for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
            // this assumes that the diagonal and off-diagonal elements of the
            // impurity and bath site hamiltonians are diagonal;
            h_nqi(kdim + it + il * eigenvaluesQ_kept_indices[ii].size(),
                  kdim + it + il * eigenvaluesQ_kept_indices[ii].size()) =
                eigenvaluesQ[ii][it] * rescale + bath_eigenvaluesQ[bb][il];
          }
        }
        size_t kpdim = 0;
        for (auto kp_index : lindex) {
          auto ii_p = kp_index / nq_bath.size(); // impurity nqi index
          auto bb_p = kp_index % nq_bath.size(); // bath nqi index
          // get f_operators of the bath and
          // f_dag_opr_list  of the previous
          // site.
          // This is optional, mean sometime the operator
          // does not exist i.e.,{0}
          for (size_t ip = 0; ip < pre_fdag_oparator.size(); ip++) {
            auto fdag_opt =
                pre_fdag_oparator[ip].get(ii, ii_p); // previous f f_operators
            auto f_bath_opt = bath_model->f_dag_operator[ip].get(
                bb_p, bb); // get the f_dag operator
            // auto f_bath_opt = bath_model->f_operator[ip].get(bb, bb_p);
            if (fdag_opt && f_bath_opt) {
              auto *f_dag_prev =
                  fdag_opt.value(); // This is a pointer to qmatrix
              auto f_bath = f_bath_opt.value();
              for (size_t it : eigenvaluesQ_kept_indices[ii]) {
                for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
                  for (size_t it_p : eigenvaluesQ_kept_indices[ii_p]) {
                    for (size_t il_p = 0; il_p < bath_eigenvaluesQ[bb_p].size();
                         il_p++) {
                      // sum over all the f operator
                      double aa = 0;
                      // std::cout << it << it_p << il << il_p << std::endl;
                      aa += (f_dag_prev->at(it, it_p)) //
                            * (f_bath->at(il_p, il))   //
                            * (chi_bath[bb_p]);
                      h_nqi(
                          kdim + it + il * eigenvaluesQ_kept_indices[ii].size(),
                          kpdim + it_p +
                              il_p * eigenvaluesQ_kept_indices[ii_p].size()) +=
                          aa * t_hopping[ip];
                      h_nqi(kpdim + it_p +
                                il_p * eigenvaluesQ_kept_indices[ii_p].size(),
                            kdim + it +
                                il * eigenvaluesQ_kept_indices[ii].size()) +=
                          aa * t_hopping[ip];
                    }
                  }
                  // end of iner loop
                }
              }
            }
            // end off local nqi index
          }
          kpdim = kpdim + eigenvaluesQ_kept_indices[ii_p].size() *
                              bath_eigenvaluesQ[bb_p].size();
        }
        kdim = kdim + eigenvaluesQ_kept_indices[ii].size() *
                          bath_eigenvaluesQ[bb].size();
      }
      current_hamiltonQ[licq] = h_nqi;
    }
  }
  void create_next_basis() {
    coupled_nQ_index.clear();
    current_sysmQ.clear();
    for (size_t i = 0; i < pre_sysmQ.size(); i++) {
      for (size_t j = 0; j < nq_bath.size(); j++) {
        auto tm = pre_sysmQ[i];
        for (size_t ix = 0; ix < nq_bath[j].size(); ix++) {
          tm[ix] += nq_bath[j][ix];
        }
        // check if tm is already exist or not
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
  // function to create in the nrgcore
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
      // size_t tdim{max_kept_states};
      for (size_t i = max_kept_states + 1; i < all_eigenvalue.size(); i++) {
        if (std::fabs(all_eigenvalue[i] - all_eigenvalue[max_kept_states]) <
            errorbarInEnergy) {
          // std::cout << "En_max" << En_max << " " << all_eigenvalue[i]
          //         << std::endl;
          // tdim++;
        }
      }
      // Set kept indices
      // check if higher energy discardtion is needed
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
    // Save the energy eigenvalues
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
    // std::cout << "eigenvaluesQ_kept_indices" << eigenvaluesQ_kept_indices;
    if (no_of_kept_states <= max_kept_states) {
      nrg_iterations_min = nrg_iterations_cnt;
    }
    // Set energy values relative to the ground state.
    else { // Shift energy of the ground state energy only if the higher
      // energy state energy discarded
      relativeGroundStateEnergy.push_back(En_min);
      for (auto &aa : eigenvaluesQ) {
        std::transform(aa.begin(), aa.end(), aa.begin(),
                       [En_min](double a) { return a - En_min; });
        // for (auto &bb : aa) {
        //     bb = bb - En_min;
        // }
      }
    }
  }
  /**
   * @brief Set the maximum numbers of states to be kept states. The
   * actual number of states is determined by adding few more degenarete
   * states.
   *
   * @param n: maximum number of states to be kept
   */
  void set_parameters(size_t n = 1024) {
    no_of_kept_states  = n;
    errorbarInEnergy   = 1e-8;
    nrg_iterations_cnt = 0;
    nrg_iterations_min = 0;
    max_kept_states    = no_of_kept_states;
  }
  // variables
  /**
   * @brief Get the current Basis vector.
   *
   * @return std::vector<std::vector<int>>
   */
  std::vector<std::vector<int>> get_basis_nQ() {
    if (current_sysmQ.empty()) {
      throw std::runtime_error("current_sysmQ.empty!");
    }
    return current_sysmQ; // because current  stuffs are move
  }
  /**
   * @brief This function returns the eigenvalues of the Hamiltonian
   * in the current basis.
   *
   * @return  std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> get_eigenvaluesQ() {
    if (eigenvaluesQ.empty()) {
      throw std::runtime_error("eigenvaluesQ.empty!");
    }
    return eigenvaluesQ;
  }
  /**
   * @brief Returns the \f f^{\dagger} \f operator in the current basis
   * current wilson site. This `qOperator` is used to construct the Hamiltonian
   * for the next iteration.
   *
   * @return std::vector<qOperator>
   */
  auto get_f_dag_operator() { return pre_fdag_oparator; }
  /**
   * @brief Check whether any states are discarded in the current iteration.
   *
   * @return bool
   */
  bool checkHigherEnergyDiscarded() {
    return no_of_kept_states >= max_kept_states;
  }

private:
  void set_current_fdag_operator() { // NOLINT
    for (auto &aa : pre_fdag_oparator) {
      aa.clear();
    }
    //
    std::map<std::array<size_t, 3>, std::vector<std::array<size_t, 6>>> foprIdx;
    //
    timer t1("set_current_fdag_operator");
    for (size_t ip = 0; ip < pre_fdag_oparator.size(); ip++) {
      // find the i,j index of the basis
      for (size_t i = 0; i < current_sysmQ.size(); i++) {
        for (size_t j = 0; j < current_sysmQ.size(); j++) {
          // TODO(sp): Not all the eigenvalues are always available
          auto idx   = coupled_nQ_index[i];
          auto idx_p = coupled_nQ_index[j];
          // std::cout << "--------------------------";
          // bool checkValue = false;
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
                // checkValue = true;
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
    // start the foprIdx Iteration
    std::cout << "--------------------------" << t1.getDuration() << std::endl;
    std::vector<qmatrix<double>> qfr(foprIdx.size());
    // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
    // #pragma omp parallel for
    for (auto itf = foprIdx.begin(); itf != foprIdx.end(); itf++) {
      //
      //
      auto            ip      = itf->first[0];
      auto            i       = itf->first[1];
      auto            j       = itf->first[2];
      size_t          kpdim   = eigenvaluesQ[i].size();
      size_t          kpdim_p = eigenvaluesQ[j].size();
      qmatrix<double> tmat(kpdim, kpdim_p, 0);
      // value
      //
      for (auto const &itd : itf->second) {
        auto ii = itd[0];
        auto bb = itd[1];
        // auto ii_p = itd[2];
        auto bb_p = itd[3];
        //
        size_t kidx   = itd[4];
        size_t kidx_p = itd[5];
        //
        // auto fdag_bath = bath_model->get_fdag_operator(bb, bb_p);
        // create previous bath id matrix
        // std::cout << "bb:bb_p" << bb << ":" << bb_p << std::endl;
        auto fdag_bath_opt = bath_model->f_dag_operator[ip].get(bb, bb_p);
        auto fdag_bath     = fdag_bath_opt.value();
        // TODO(sp): check the pragma openmp effect for the multi orbital
        //
        // #pragma omp parallel for collapse(2)
        for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
          for (size_t il_p = 0; il_p < bath_eigenvaluesQ[bb_p].size(); il_p++) {
            double rvalue = fdag_bath->at(il, il_p);
            for (size_t it : eigenvaluesQ_kept_indices[ii]) {
              tmat(kidx + it + il * eigenvaluesQ_kept_indices[ii].size(),
                   kidx_p + it + il_p * eigenvaluesQ_kept_indices[ii].size()) +=
                  rvalue;
              // fdag_bath->at(il, il_p);
            }
          }
        }
      } // set the operator
      // pre_fdag_oparator[ip].set(current_hamiltonQ[i].cTranspose().dot(tmat.dot(
      //                               current_hamiltonQ[j])), //
      //                               UnitaryTransform
      //                           i, j);
      qfr[std::distance(foprIdx.begin(), itf)] =
          current_hamiltonQ[i].cTranspose().dot(tmat.dot(current_hamiltonQ[j]));
      // end of ip loop
    }
    std::cout << "--------------------------" << t1.getDuration() << std::endl;
    // set the foparator
    for (auto itf = foprIdx.begin(); itf != foprIdx.end(); itf++) {
      //
      //
      auto ip = itf->first[0];
      auto i  = itf->first[1];
      auto j  = itf->first[2];
      pre_fdag_oparator[ip].set(qfr[std::distance(foprIdx.begin(), itf)], i, j);
    }
  }
  void enforceDegeneracy() {
    // enforce the degenarecy of the energy levels
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
  void set_current_fdag_operator_old() { // NOLINT
    for (auto &aa : pre_fdag_oparator) {
      aa.clear();
    }
    //
    for (size_t i = 0; i < current_sysmQ.size(); i++) {
      for (size_t j = 0; j < current_sysmQ.size(); j++) {
        // TODO(sp): Not all the eigenvalues are always available
        size_t kpdim   = eigenvaluesQ[i].size();
        size_t kpdim_p = eigenvaluesQ[j].size();
        auto   idx     = coupled_nQ_index[i];
        auto   idx_p   = coupled_nQ_index[j];
        // std::cout << "--------------------------";
        // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
        std::vector<qmatrix<double>> c_fopr(pre_fdag_oparator.size(),
                                            qmatrix<double>(kpdim, kpdim_p, 0));
        std::vector<bool>            hasValue(pre_fdag_oparator.size(), false);
        size_t                       kidx = 0;
        for (auto kindex : idx) {
          auto   ii     = kindex / nq_bath.size(); // impurity nqi index
          auto   bb     = kindex % nq_bath.size(); // bath nqi index
          size_t kidx_p = 0;
          for (auto kindex_p : idx_p) {
            auto ii_p = kindex_p / nq_bath.size(); // impurity nqi index
            auto bb_p = kindex_p % nq_bath.size(); // bath nqi index
            if (ii == ii_p) {
              // auto fdag_bath = bath_model->get_fdag_operator(bb, bb_p);
              // create previous bath id matrix
              for (size_t ip = 0; ip < pre_fdag_oparator.size(); ip++) {
                auto fdag_bath_opt =
                    bath_model->f_dag_operator[ip].get(bb, bb_p);
                if (fdag_bath_opt) {
                  // std::cout << "bb:bb_p" << bb << ":" << bb_p << std::endl;
                  hasValue[ip]   = true;
                  auto fdag_bath = fdag_bath_opt.value();
                  for (size_t il = 0; il < bath_eigenvaluesQ[bb].size(); il++) {
                    for (size_t il_p = 0; il_p < bath_eigenvaluesQ[bb_p].size();
                         il_p++) {
                      for (size_t it : eigenvaluesQ_kept_indices[ii]) {
                        c_fopr[ip](
                            kidx + it +
                                il * eigenvaluesQ_kept_indices[ii].size(),
                            kidx_p + it +
                                il_p * eigenvaluesQ_kept_indices[ii].size()) +=
                            fdag_bath->at(il, il_p);
                      }
                    }
                  }
                }
              }
            }
            kidx_p += eigenvaluesQ_kept_indices[ii_p].size() *
                      bath_eigenvaluesQ[bb_p].size();
          }
          kidx += eigenvaluesQ_kept_indices[ii].size() *
                  bath_eigenvaluesQ[bb].size();
        }
        // set the matrix elements
        // end of lm loop
        // End of matrix generation.
        // Rotate the c operator in the eigen basis
        for (size_t iv = 0; iv < hasValue.size(); iv++) {
          if (hasValue[iv]) {
            pre_fdag_oparator[iv].set(
                current_hamiltonQ[i].cTranspose().dot(
                    c_fopr[iv].dot(current_hamiltonQ[j])), // UnitaryTransform
                i, j);
          }
        } // End of loop save
          // std::cout << "c_fopr" << c_fopr;
      }
    }
  }
  // Store and set previous hamiltonian
  // indices of of coupled system and bath n_Q i.e.,
  // which subspaces are connected
  // get the bath stuff so that we dont have to create them
  // in every site
  std::vector<double> chi_bath;
  // Construct the nrgcore
  // from two hamiltonian type
  size_t no_of_kept_states{0};
  size_t max_kept_states{0};
  double errorbarInEnergy{0};

public:
  std::vector<double> all_eigenvalue;
  std::vector<double> relativeGroundStateEnergy;
  // No need save this for backward iteration.
  std::vector<std::vector<size_t>> eigenvaluesQ_kept_indices;
  std::vector<qmatrix<double>>     current_hamiltonQ; // next hamiltonians
  std::vector<std::vector<int>>    current_sysmQ;     // next symmetries
  std::vector<std::vector<int>>    pre_sysmQ;         // previous symmetries
  std::vector<std::vector<double>> eigenvaluesQ;      // Eigenvalues
  std::vector<std::vector<size_t>> coupled_nQ_index;
  // These are not update on each wilson site
  int                              nrg_iterations_cnt{};
  int                              nrg_iterations_min{};
  std::vector<std::vector<double>> bath_eigenvaluesQ;
  std::vector<std::vector<int>>    nq_bath;
  std::vector<qOperator>          *getPreWilsonSiteOperators() {
    return &pre_fdag_oparator;
  }
};
