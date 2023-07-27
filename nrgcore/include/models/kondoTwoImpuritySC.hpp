#pragma once
#include "models/fermionBasis.hpp"
#include "nrgcore/qOperator.hpp"
#include "nrgcore/qsymmetry.hpp"
#include "utils/qmatrix.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>
class kondoTwoImpuritySC : public fermionBasis {
  /** This class is for a single orbital with spin up and down
   * f operator. SIAM can be made entirely from this class.
   *
   *
   */
public:
  /**
   * @brief Construct a new kondoTwoImpuritySC.hpp object
   *
   * @param JKondo: J value for the Kondo interaction
   * @param spinS: Spin value i.e,. 1/2 or 3/2. Integer Spin may not work
   */
  //    ########################################
  explicit kondoTwoImpuritySC(const std::map<std::string, double> &params) {
    double JKondo1    = params.at("Jkondo1");
    double JKondo2    = params.at("Jkondo2");
    double Jrkky      = params.at("J_RKKY");
    double spinS      = params.at("spinS");
    double localDelta = params.at("Delta_sc");
    auto   spinDim    = static_cast<size_t>(2. * spinS + 1.);
    // TODO: check spinS is multiple of 2.
    if (std::fabs(spinDim - (2. * spinS + 1.)) > 0
        //
        //  or zero
        || abs(spinS) == 0) {
      throw std::invalid_argument(
          "spinS must be a half integer or integer. spinS:" +
          std::to_string(spinS));
    }
    std::cout << "spinDim: " << spinDim << std::endl;
    qmatrix<double> spinXz(spinDim, spinDim, 0);
    qmatrix<double> spinXplus(spinDim, spinDim, 0);
    // Spin Z
    for (size_t i = 0; i < spinDim; i++) {
      double m     = -spinS + i * 1.0;
      spinXz(i, i) = m;
      if (i != (spinDim - 1)) {
        std::cout << "spinDim: " << i << " " << spinDim << std::endl;
        spinXplus(i, i + 1) = std::sqrt(spinS * (spinS + 1.) - m * (m + 1));
      }
    }
    // std::cout << "spinSz: \n" << spinSz << std::endl;
    // std::cout << "spinSplus: \n" << spinSplus << std::endl;
    auto spinXminus = spinXplus.cTranspose();
    auto spinXx     = (spinXplus + spinXminus) * 0.5;
    auto spinXy     = (spinXplus - spinXminus) * 0.5; // -i is omitted here.
    auto spinXId    = spinXz.id();
    // ###########################################
    auto spin1Sz = spinXz.krDot(spinXId);
    auto spin1Sx = spinXx.krDot(spinXId);
    auto spin1Sy = spinXy.krDot(spinXId);
    // second spin
    auto spin2Sz = spinXId.krDot(spinXz);
    auto spin2Sx = spinXId.krDot(spinXx);
    auto spin2Sy = spinXId.krDot(spinXy);
    // ########################################
    // std::cout << "spinSx: " << spinSx << std::endl;
    // std::cout << "spinSy: " << spinSy << std::endl;
    // std::cout << "spinSz: " << spinSz << std::endl;
    // ########################################
    // create basis
    qmatrix<> fdag = {0, 0, 1, 0};
    qmatrix<> sigz = {1, 0, 0, -1};
    qmatrix<> id2  = {1, 0, 0, 1};
    // Create c_up and c_Down operator
    // Addting and additional operators should be done in the same way.
    auto c_up_dag   = fdag.krDot(id2);
    auto c_Down_dag = sigz.krDot(fdag);
    // Create Wilson Site spin operators
    auto wSpinx = (c_up_dag.dot(c_Down_dag.cTranspose()) +
                   c_Down_dag.dot(c_up_dag.cTranspose())) *
                  0.5;
    auto wSpiny = (c_Down_dag.dot(c_up_dag.cTranspose()) -
                   c_up_dag.dot(c_Down_dag.cTranspose())) *
                  0.5; // -i is omitted here.
    auto wSpinz = (c_up_dag.dot(c_up_dag.cTranspose()) -
                   c_Down_dag.dot(c_Down_dag.cTranspose())) *
                  0.5;
    auto wNtotal =
        (c_up_dag.dot(c_up_dag.cTranspose()) +
         c_Down_dag.dot(
             c_Down_dag.cTranspose())); // Create Hamiltonian spin operator
    // ########################################################
    // std::cout << "wSpinx: " << wSpinx << std::endl;
    // std::cout << "wSpiny: " << wSpiny << std::endl;
    // std::cout << "wSpinz: " << wSpinz << std::endl;
    // ############################################
    auto Hamiltonian =
        // first spin coupling
        ((spin1Sz.krDot(wSpinz) + spin1Sx.krDot(wSpinx) -
          spin1Sy.krDot(wSpiny)) *
         JKondo1) +
        // second spin coupling
        ((spin2Sz.krDot(wSpinz) + spin2Sx.krDot(wSpinx) //
          - spin2Sy.krDot(wSpiny)) * // Imaginary part is taken care here
         JKondo2)
        // rkky interaction between the two impurities
        + (spin1Sz.dot(spin2Sz) + spin1Sx.dot(spin2Sx) - spin1Sy.dot(spin2Sy))
                  .krDot(wSpinz.id()) *
              Jrkky
        // SC term
        + spin1Sz.id().krDot(
              (c_up_dag.dot(c_Down_dag) +
               c_Down_dag.cTranspose().dot(c_Down_dag.cTranspose()))) *
              localDelta //
        ;
    // End
    // std::cout << "spinSz: " << spinSz << std::endl;
    auto qspin1Sz = spin1Sz.krDot(wSpinx.id());
    auto qspin2Sz = spin2Sz.krDot(wSpinx.id());
    // std::cout << "spinSz: " << spinSz << std::endl;
    auto qwSpinz  = spin1Sx.id().krDot(wSpinz);
    auto qwNtotal = spin1Sx.id().krDot(wNtotal);
    // std::cout << "wSpinz: " << wSpinz * 2 << std::endl;
    c_up_dag   = spin1Sx.id().krDot(c_up_dag);
    c_Down_dag = spin1Sx.id().krDot(c_Down_dag);
    fnParticle.clear();
    // Number of particles
    fnParticle.push_back(((qspin1Sz.id() + qspin1Sz * 2.) * 0.5).getdiagonal());
    fnParticle.push_back(((qspin1Sz.id() - qspin1Sz * 2.) * 0.5).getdiagonal());
    //
    fnParticle.push_back(((qspin1Sz.id() + qspin2Sz * 2.) * 0.5).getdiagonal());
    fnParticle.push_back(((qspin1Sz.id() - qspin2Sz * 2.) * 0.5).getdiagonal());
    //
    fnParticle.push_back(((qwNtotal + qwSpinz * 2.) * 0.5).getdiagonal());
    fnParticle.push_back(((qwNtotal - qwSpinz * 2.) * 0.5).getdiagonal());
    //
    //
    std::cout << "fnParticle: " << fnParticle << std::endl;
    // set the symmetries of the system
    // create_QuantumNspinCharge();
    // spin is only quantum number that is conserved
    std::vector<size_t> tm1;
    std::vector<size_t> tm2;
    for (size_t j = 0; j < fnParticle.size(); j++) {
      if (j % 2 == 0) { // Spin up
        tm1.push_back(j);
      } else { // Spin down
        tm2.push_back(j);
      }
    }
    std::vector<std::vector<size_t>> spinIdx = {tm1, tm2};
    create_QuantumSpinOnly(spinIdx);
    // create the quantum numbers
    // create the quantum numbers
    create_Block_structure();
    // ####################################################################
    n_Q = get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    //
    // set foperator
    auto h_blocked = get_block_Hamiltonian(Hamiltonian);
    // std::cout << "h_blocked: " << h_blocked << std::endl;
    // std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
    }
    // TODO: rotate the f operator
    // ####################################################################
    f_dag_operator = get_block_operators({c_up_dag, c_Down_dag});
    std::cout << "f_dag_operators: " << f_dag_operator.size() << std::endl;
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
    f_dag_operator = topr;
  }
  //    ######################################
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; }
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }
  [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }
  //
  std::vector<std::vector<double>> eigenvalues_Q;
  std::vector<double>              chi_Q;
  std::vector<std::vector<int>>    n_Q;
};
