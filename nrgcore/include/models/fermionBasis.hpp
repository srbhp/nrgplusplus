
#pragma once
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <nrgcore/qsymmetry.hpp>
#include <ostream>
#include <set>
/**
 * @class fermionBasis
 * @brief This is the main basis class for the electrons or fermions.
 * Every other fermion model class is derived from this class.
 * We normally dont inherit from this class. Just use as a local object.
 * See for example spinhalf.
 *
 */
class fermionBasis {
public:
  /**
   * @brief Symmetries used used in the proble
   * chargeOnly: Only charge of the system is conserved.
   * spinOnly : Only \f$S_z\f$ component of the total spin is coserved. This is
   * use for the superconducting Bath.
   * chargeAndSpin: Both charge and spin is conserved
   */
  enum modelSymmetry {
    chargeOnly, // Only charge of the system is conserved.
    spinOnly,
    chargeAndSpin
  }; // Add a key for the exact diagonisation
  /**
   * @brief \f$ f^\dagger \f$ of the system written using symmetries of the
   * problem.
   */
  std::vector<qOperator> f_dag_operator; // f_dag operator
  /**
   * @brief Quantization numbers used to write various system operators
   */
  std::vector<std::vector<double>> fnParticle;
  /**
   * @brief
   *
   * @param dof : Fermion degree of freedom
   * @param l : Symmetries of the system
   */
  fermionBasis(size_t dof, modelSymmetry l) {
    // create the basis in nstates x nstates
    createFermionBasis(dof); // create the basis in nstates x nstates
    // create the switch statement for l
    switch (l) {
    case chargeOnly:
      std::cout << "chargeOnly" << std::endl;
      create_QuantumNChargeonly(); // create the quantum numbers
      break;
    case spinOnly:
      std::cout << "spinOnly" << std::endl;
      {
        std::vector<size_t> tm1;
        std::vector<size_t> tm2;
        for (size_t j = 0; j < dof; j++) {
          if (j % 2 == 0) { // Spin up
            tm1.push_back(j);
          } else { // Spin down
            tm2.push_back(j);
          }
        }
        std::vector<std::vector<size_t>> SpinSz = {tm1, tm2};
        create_QuantumSpinOnly(SpinSz); // create the quantum numbers
      }                                 //
      break;
    case chargeAndSpin:
      std::cout << "chargeAndSpin" << std::endl;
      create_QuantumNspinCharge();
      break;
    }
    create_Block_structure(); // create the block structure
    set_f_dag_operators();    // set the fermion operators
    DebugView = false;        //
  }
  /**
   * @brief Get the full Matrix of the `f` operators.
   *
   * @return
   */
  [[nodiscard]] auto get_raw_f_dag_operator() const { return fermionOprMat; }
  /**
   * @brief Get `f` operators of the site in terms of the `qmatrix` notations
   *
   * @return
   */
  [[nodiscard]] auto get_f_dag_operator() const { return f_dag_operator; }
  void               create_QuantumNChargeonly() {
    // Here We assume that
    // charge of the individual spin is conserved
    nQuantum.clear();
    for (size_t i = 0; i < fnParticle[0].size(); i++) {
      int tm1{0};
      for (auto &j : fnParticle) {
        tm1 += static_cast<int>(j[i]);
        // std::cout << fnpr[j][i] << " ";
      }
      nQuantum.push_back({tm1});
      // std::cout << ":Q:" << std::vector({tm1, tm2});
      // std::cout << std::endl;
    }
  }
  void create_QuantumNspinCharge() {
    // Here We assume that
    // charge of the individual spin is conserved
    nQuantum.clear();
    for (size_t i = 0; i < fnParticle[0].size(); i++) {
      int tm1{0};
      int tm2{0};
      for (size_t j = 0; j < fnParticle.size(); j++) {
        if (j % 2 == 0) { // Spin up
          tm1 += static_cast<int>(fnParticle[j][i]);
        } else { // Spin down
          tm2 += static_cast<int>(fnParticle[j][i]);
        }
        // std::cout << fnpr[j][i] << " ";
      }
      nQuantum.push_back({tm1, tm2});
      // std::cout << ":Q:" << std::vector({tm1, tm2});
      // std::cout << std::endl;
    }
  }
  void create_QuantumSpinOnly(std::vector<std::vector<size_t>> &qsymmetry) {
    // qsymmetry defines the symmetries of propblem
    // qsymmetry.size = 2 . If there is four fermion flavours the qsymmetry is
    // (Charge of individual spin is conserved)
    // i.e  qsymmetry = {{0,2},{1,3}}
    // check for empty
    if (qsymmetry.empty()) {
      std::cout << "Warning: qsymmetry is empty" << std::endl;
      // create odd even pair for spin up and down.
      std::vector<size_t> tm1;
      std::vector<size_t> tm2;
      for (size_t j = 0; j < fnParticle.size(); j++) {
        if (j % 2 == 0) { // Spin up
          tm1.push_back(j);
        } else { // Spin down
          tm2.push_back(j);
        }
      }
      qsymmetry = std::vector<std::vector<size_t>>({tm1, tm2});
    }
    if (qsymmetry.size() != 2) {
      throw std::runtime_error("qsymmetry.size() != 2");
    }
    nQuantum.clear();
    for (size_t i = 0; i < fnParticle[0].size(); i++) {
      std::vector<int> tmVec;
      {
        size_t iq = 0;
        int    tm1{0}; // up spin
        int    tm2{0}; // down spin
        for (auto jq : qsymmetry[iq]) {
          tm1 += static_cast<int>(fnParticle[jq][i]);
        }
        iq = 1;
        for (auto jq : qsymmetry[iq]) {
          tm2 += static_cast<int>(fnParticle[jq][i]);
        }
        tmVec.push_back(tm1 - tm2);
      }
      nQuantum.push_back(tmVec);
      // std::cout << ":Q:" << std::vector({tm1, tm2});
      // std::cout << std::endl;
    }
  }
  void createQNumbers(const std::vector<std::vector<size_t>> &qsymmetry) {
    // qsymmetry defines the symmetries of propblem
    // qsymmetry.size = 2 for a system with charge and conserved
    // (Charge of individual spin is conserved)
    // i.e  qsymmetry = {{0,2},{1,3}}
    nQuantum.clear();
    for (size_t i = 0; i < fnParticle[0].size(); i++) {
      std::vector<int> tmVec;
      for (const auto &iq : qsymmetry) {
        int tm1{0};
        for (auto jq : iq) {
          tm1 += static_cast<int>(fnParticle[jq][i]);
        }
        tmVec.push_back(tm1);
      }
      nQuantum.push_back(tmVec);
      // std::cout << ":Q:" << std::vector({tm1, tm2});
      // std::cout << std::endl;
    }
  }
  void create_Block_structure() {
    nQBlocks.clear();
    unique_Qnumbers.clear();
    // unique_Indices of the quantum numbers
    std::vector<size_t> unique_Indices;
    {
      std::set<std::vector<int>> tm;
      for (size_t i = 0; i < fnParticle[0].size(); i++) {
        if (tm.insert(nQuantum[i]).second) {
          unique_Indices.push_back(i);
        }
      }
    }
    // sub blocks indices where Quantum numbers are same.
    for (auto i : unique_Indices) {
      unique_Qnumbers.push_back(nQuantum[i]);
      std::vector<size_t> tm1{i};
      for (size_t j = i; j < fnParticle[0].size(); j++) {
        if (nQuantum[i] == nQuantum[j]) {
          if (i != j) {
            tm1.push_back(j);
          }
        }
      }
      if (DebugView) {
        std::cout << "Coupled : " << nQuantum[i] << " : ";
        for (auto j : tm1) {
          for (auto &jx : fnParticle) {
            std::cout << jx[j] << " ";
          }
          std::cout << std::endl;
        }
      }
      // << tm1 << std::endl;
      nQBlocks.push_back(tm1);
    }
    // std::cout << "nQBlocks : " << nQBlocks.size() << std::endl;
    //
  }
  auto get_unique_Qnumbers() {
    if (unique_Qnumbers.empty()) {
      std::cout << "Warning: unique_Qnumbers is empty !" << std::endl;
    }
    return unique_Qnumbers;
  }
  std::vector<qOperator>
  get_block_operators(const std::vector<qmatrix<>> &sys_operators) {
    // This function can be used to calculate any other
    // operators which are a combinations of the f operators
    std::vector<qOperator> block_operators(sys_operators.size(), qOperator());
    for (size_t i = 0; i < nQBlocks.size(); i++) {
      for (size_t j = 0; j < nQBlocks.size(); j++) {
        size_t tdim   = nQBlocks[i].size();
        size_t tdim_p = nQBlocks[j].size();
        for (size_t ipr = 0; ipr < sys_operators.size(); ipr++) {
          qmatrix<> fmatrix(tdim, tdim_p, 0);
          double    tvalue{0};
          for (size_t ik = 0; ik < nQBlocks[i].size(); ik++) {
            for (size_t ik_p = 0; ik_p < nQBlocks[j].size(); ik_p++) {
              fmatrix(ik, ik_p) =
                  sys_operators[ipr](nQBlocks[i][ik], nQBlocks[j][ik_p]);
              tvalue += std::fabs(
                  sys_operators[ipr](nQBlocks[i][ik], nQBlocks[j][ik_p]));
            }
          }
          if (tvalue > 0) {
            block_operators[ipr].set(fmatrix, i, j);
          }
        }
      }
    } // End of operator
    return block_operators;
  }
  qOperator get_block_Hamiltonian(const qmatrix<double> &sys_operators) {
    // This function can be used to calculate any other
    // check the Hamiltonian is Harmitian i.e H == H.T
    {
      double tvalue = 0;
      for (size_t i = 0; i < sys_operators.getcolumn(); i++) {
        for (size_t j = i + 1; j < sys_operators.getrow(); j++) {
          tvalue += std::fabs(sys_operators(i, j) - sys_operators(j, i));
        }
      }
      if (tvalue > 1e-10) {
        throw std::runtime_error("Hamiltonian is not Hermitian\n ");
      }
    }
    // operators which are a combinations of the f operators
    //   std::cout << "&sys_operators" << sys_operators;
    qOperator block_operators;
    for (size_t i = 0; i < nQBlocks.size(); i++) {
      size_t    tdim = nQBlocks[i].size();
      qmatrix<> fmatrix(tdim, tdim, 0);
      for (size_t ik = 0; ik < nQBlocks[i].size(); ik++) {
        for (size_t ik_p = 0; ik_p < nQBlocks[i].size(); ik_p++) {
          fmatrix(ik, ik_p) = sys_operators(nQBlocks[i][ik], nQBlocks[i][ik_p]);
        }
      }
      block_operators.set(fmatrix, i, i);
    } //
      // test all the offdiagonal elements are zero
    double tvalue = 0;
    for (size_t i = 0; i < nQBlocks.size(); i++) {
      for (size_t j = i + 1; j < nQBlocks.size(); j++) {
        for (size_t ik = 0; ik < nQBlocks[i].size(); ik++) {
          for (size_t ik_p = 0; ik_p < nQBlocks[j].size(); ik_p++) {
            tvalue += sys_operators(nQBlocks[i][ik], nQBlocks[j][ik_p]);
          }
        }
      }
    }
    if (std::abs(tvalue) > 1e-10) {
      throw std::runtime_error("Warning: Hamiltonian is not block diagonal \n");
    }
    //
    // / ///////////////////////////////////////
    return block_operators;
  }
  // End of operator
  [[nodiscard]] auto get_basis() const { return fnParticle; }
  void               setDebugMode(bool debug = true) { DebugView = debug; }
  void               createFermionBasis(size_t ldof) {
    qmatrix<double> fdag({0, 0, 1, 0});
    qmatrix<double> sigz({1, 0, 0, -1});
    qmatrix<double> id2({1, 0, 0, 1});
    fermionOprMat.clear();
    for (size_t i = 0; i < ldof; i++) {
      qmatrix<> prev(1, 1, 1);
      for (size_t j = 0; j < ldof - i - 1; j++) {
        prev = id2.krDot(prev);
      }
      qmatrix<> nxt(1, 1, 1);
      for (size_t j = 0; j < i; j++) {
        nxt = sigz.krDot(nxt);
      }
      fermionOprMat.push_back(prev.krDot(fdag.krDot(nxt)));
    }
    fnParticle.clear();
    //  std::cout << fup << fdw;
    for (size_t i = 0; i < ldof; i++) {
      fnParticle.push_back(
          (fermionOprMat[i].dot(fermionOprMat[i].cTranspose())).getdiagonal());
    }
    //
  }
  std::vector<qmatrix<>> fermionOprMat;
  void                   set_f_dag_operators() {
    f_dag_operator = get_block_operators(fermionOprMat);
  }

private:
  bool DebugView{false};
  // foperator in the full basis
  // no of particles operator
  // basis states -- and quantum numbers
  // n_up n_down as a quntum number
  std::vector<std::vector<int>> nQuantum;
  std::vector<std::vector<int>> unique_Qnumbers;
  // Blocked no of particles i.e.,
  std::vector<std::vector<size_t>> nQBlocks;
};
