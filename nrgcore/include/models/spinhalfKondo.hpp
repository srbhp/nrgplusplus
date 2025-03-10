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
class spinhalfKondo {
  /** This class is for a single orbital with spin up and down
   * f operator. SIAM can be made entirely from this class.
   *
   *
   */
public:
  /**
   * @brief Construct a new spinhalfKondo object
   *
   * @param JKondo: J value for the Kondo interaction
   * @param spinS: Spin value i.e,. 1/2 or 3/2. Integer Spin may not work
   */
  spinhalfKondo(double JKondo, double spinS) {
    auto nstates = static_cast<size_t>(
        (2. * spinS + 1) * 4); // std::pow(2, dof); // Dimension of the system
    auto dof = static_cast<size_t>(
        (2. * spinS + 1) + 2); // std::pow(2, dof); // Dimension of the system
    std::cout << "nstates: " << nstates << "\n"
              << "dof: " << dof << "\n";
    createBasis(JKondo, spinS); // create the basis in nstates x nstates
  }
  [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; }
  [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
    return eigenvalues_Q;
  }
  [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }
  //
  std::vector<std::vector<double>> eigenvalues_Q;
  std::vector<double>              chi_Q;
  std::vector<qOperator>           f_dag_operator;
  std::vector<std::vector<int>>    n_Q;
  //    ########################################
private:
  void createBasis(double JKondo, double spinS) { // NOLINT
    // TODO(sp): check spinS is multiple of 2.
    auto spinDim = static_cast<size_t>(2. * spinS + 1.);
    std::cout << "spinDim: " << spinDim << std::endl;
    qmatrix<double> spinSz(spinDim, spinDim, 0);
    qmatrix<double> spinSplus(spinDim, spinDim, 0);
    // Spin Z
    for (size_t i = 0; i < spinDim; i++) {
      double m     = -spinS + static_cast<double>(i) * 1.0;
      spinSz(i, i) = m;
      if (i != (spinDim - 1)) {
        std::cout << "spinDim: " << i << " " << spinDim << std::endl;
        spinSplus(i, i + 1) = std::sqrt(spinS * (spinS + 1.) - m * (m + 1));
      }
    }
    std::cout << "spinSz: \n" << spinSz << std::endl;
    std::cout << "spinSplus: \n" << spinSplus << std::endl;
    auto spinSminus = spinSplus.cTranspose();
    auto spinSx     = (spinSplus + spinSminus) * 0.5;
    auto spinSy     = (spinSplus - spinSminus) * 0.5; // -i is omitted here.
    // ########################################
    std::cout << "spinSx: " << spinSx << std::endl;
    std::cout << "spinSy: " << spinSy << std::endl;
    std::cout << "spinSz: " << spinSz << std::endl;
    // ########################################
    // create basis
    fermionBasis fnBasis(2, fermionBasis::chargeAndSpin);
    // Create c_up and c_Down operator
    // Addting and additional operators should be done in the same way.
    auto c_up_dag   = fnBasis.fermionOprMat[0];
    auto c_Down_dag = fnBasis.fermionOprMat[1];
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
    std::cout << "wSpinx: " << wSpinx << std::endl;
    std::cout << "wSpiny: " << wSpiny << std::endl;
    std::cout << "wSpinz: " << wSpinz << std::endl;
    // ############################################
    auto Hamiltonian = (spinSz.krDot(wSpinz) + spinSx.krDot(wSpinx) //
                        - spinSy.krDot(wSpiny)) *
                       // Imaginary part is taken care here
                       JKondo;
    // End
    fnBasis.fnParticle.clear();
    // std::cout << "spinSz: " << spinSz << std::endl;
    spinSz = qmatrix(spinSz.krDot(wSpinx.id()));
    // std::cout << "spinSz: " << spinSz << std::endl;
    wSpinz  = spinSx.id().krDot(wSpinz);
    wNtotal = spinSx.id().krDot(wNtotal);
    std::cout << "wSpinz: " << wSpinz * 2 << std::endl;
    c_up_dag   = spinSx.id().krDot(c_up_dag);
    c_Down_dag = spinSx.id().krDot(c_Down_dag);
    // Number of particles
    fnBasis.fnParticle.push_back(
        ((spinSz.id() + spinSz * 2.) * 0.5).getdiagonal());
    fnBasis.fnParticle.push_back(
        ((spinSz.id() - spinSz * 2.) * 0.5).getdiagonal());
    fnBasis.fnParticle.push_back(((wNtotal + wSpinz * 2.) * 0.5).getdiagonal());
    fnBasis.fnParticle.push_back(((wNtotal - wSpinz * 2.) * 0.5).getdiagonal());
    //
    std::cout << "fnParticle: " << fnBasis.fnParticle << std::endl;
    //
    fnBasis.create_QuantumNspinCharge();
    fnBasis.create_Block_structure();
    //####################################################################
    n_Q = fnBasis.get_unique_Qnumbers();
    // set chi_Q
    chi_Q.clear();
    for (auto ai : n_Q) {
      double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
      chi_Q.push_back(std::pow(-1., t_charge));
    }
    //
    // set foperator
    auto h_blocked = fnBasis.get_block_Hamiltonian(Hamiltonian);
    std::cout << "h_blocked: " << h_blocked << std::endl;
    std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
    // Diagonalize the hamilton
    eigenvalues_Q.clear();
    eigenvalues_Q.resize(n_Q.size(), {});
    for (size_t i = 0; i < n_Q.size(); i++) {
      eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
      std::cout << "eigenvalues_Q[i]: " << eigenvalues_Q[i] << std::endl;
    }
    // TODO(sp): rotate the f operator
    //####################################################################
    f_dag_operator = fnBasis.get_block_operators({c_up_dag, c_Down_dag});
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
};
