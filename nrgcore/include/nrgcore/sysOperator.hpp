#pragma once
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include "utils/timer.hpp"
#include <cstddef>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <tuple>
#include <vector>

/**
 * @brief Updates the system operators for the NRG object.
 *
 * This function iterates over the quantum numbers and updates the system
 * operators based on the coupled quantum indices and bath eigenvalues.
 *
 * @tparam nrgcore_type The type of the NRG core object.
 * @param nrg_object Pointer to the NRG core object.
 * @param systemo_oparator_nQ Vector of system operators to be updated.
 */
template <typename nrgcore_type> // nrgcore_type is a type of
void update_system_operator(nrgcore_type           *nrg_object, // NOLINT
                            std::vector<qOperator> &systemo_oparator_nQ) {
  size_t                 f_operators_size = systemo_oparator_nQ.size();
  std::vector<qOperator> updated_system_operator(systemo_oparator_nQ.size(),
                                                 qOperator());
  // calc th3e cummelative dimension
  // std::cout << "foperator size: " << f_operators_size
  //           << "Sysq Size :: " << nrg_object->current_sysmQ.size() <<
  //           std::endl;
  std::vector<std::vector<size_t>> cummulative_dim(
      nrg_object->current_sysmQ.size(), std::vector<size_t>{});
  for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
    size_t kidx = 0;
    for (auto kindex : nrg_object->coupled_nQ_index[i]) {
      cummulative_dim[i].push_back(kidx);
      auto ii = kindex / nrg_object->nq_bath.size(); // impurity nqi index
      auto bb = kindex % nrg_object->nq_bath.size(); // bath nqi index
      kidx += nrg_object->eigenvaluesQ_kept_indices[ii].size() *
              nrg_object->bath_eigenvaluesQ[bb].size();
    }
  }
  // #pragma omp parallel for collapse(2) schedule(runtime)
  // // Does notimprove
  // improve
  for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
    for (size_t j = 0; j < nrg_object->current_sysmQ.size(); j++) {
      // mkl_set_num_threads_local(1);
      size_t kpdim   = nrg_object->eigenvaluesQ[i].size();
      size_t kpdim_p = nrg_object->eigenvaluesQ[j].size();
      // std::cout << "--------------------------";
      // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
      qmatrix<>              tmat(kpdim, kpdim_p, 0);
      std::vector<qmatrix<>> c_fopr(f_operators_size, tmat);
      std::vector<bool>      hasValue(f_operators_size, false);
      // TODO(sp): This
      for (size_t cdx = 0; cdx < nrg_object->coupled_nQ_index[i].size();
           cdx++) {
        for (size_t cdx_p = 0; cdx_p < nrg_object->coupled_nQ_index[j].size();
             cdx_p++) {
          auto kidx     = cummulative_dim[i][cdx];
          auto kindex   = nrg_object->coupled_nQ_index[i][cdx];
          auto ii       = kindex / nrg_object->nq_bath.size(); // impurity nqi
          auto bb       = kindex % nrg_object->nq_bath.size(); // bath nqi
          auto kidx_p   = cummulative_dim[j][cdx_p];
          auto kindex_p = nrg_object->coupled_nQ_index[j][cdx_p];
          auto ii_p =
              kindex_p / nrg_object->nq_bath.size(); // impurity nqi index
          auto bb_p = kindex_p % nrg_object->nq_bath.size(); // bath nqi index
          if (bb == bb_p) {
            // #pragma omp parallel for // DOES NOT improve
            for (size_t ip = 0; ip < f_operators_size; ip++) {
              auto sys_opr_opt = systemo_oparator_nQ[ip].get(ii, ii_p);
              if (sys_opr_opt) {
                auto sys_opr = sys_opr_opt.value();
                hasValue[ip] = true;
                // create previous bath id matrix
                for (size_t il = 0;
                     il < nrg_object->bath_eigenvaluesQ[bb].size(); il++) {
                  // We may parallelthe next two loop
                  for (size_t it : nrg_object->eigenvaluesQ_kept_indices[ii]) {
                    for (size_t it_p :
                         nrg_object->eigenvaluesQ_kept_indices[ii_p]) {
                      c_fopr[ip](
                          kidx + it +
                              nrg_object->eigenvaluesQ_kept_indices[ii].size() *
                                  il,
                          kidx_p + it_p +
                              nrg_object->eigenvaluesQ_kept_indices[ii_p]
                                      .size() *
                                  il) += sys_opr->at(it, it_p);
                    }
                  }
                }
              }
            }
          }
        }
      }
      // set the matrix elements
      // end of lm loop
      // End of matrix generation.
      // Rotate the c operator in the eigen basis
      for (size_t iv = 0; iv < hasValue.size(); iv++) {
        if (hasValue[iv]) {
          updated_system_operator[iv].set(
              nrg_object->current_hamiltonQ[i].cTranspose().dot(c_fopr[iv].dot(
                  nrg_object->current_hamiltonQ[j])), // UnitaryTransform
              i, j);
        }
      }
      // Rotate the c operator in the eigen basis
      // std::cout << "c_fopr" << c_fopr;
    }
  }
  // move the operator
  systemo_oparator_nQ = updated_system_operator;
  // std::cout << "Done ! update_system_operatorQ" << std::endl;
} // End of  update_system_operatorQ
