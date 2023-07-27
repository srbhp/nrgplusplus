
.. _program_listing_file_nrgcore_include_models_kondoSC.hpp:

Program Listing for File kondoSC.hpp
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_models_kondoSC.hpp>` (``nrgcore/include/models/kondoSC.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   class kondoSC {
     fermionBasis localSCbabsis;
   public:
     kondoSC(double JKondo, double spinS, double localDelta)
         : localSCbabsis(2, fermionBasis::spinOnly) {
       createBasis(JKondo, spinS,
                   localDelta); // create the basis in nstates x nstates
     }
     [[nodiscard]] std::vector<std::vector<int>> get_basis() const { return n_Q; }
     [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
       return eigenvalues_Q;
     }
     [[nodiscard]] std::vector<double> get_chi_Q() const { return chi_Q; }
     //
     std::vector<std::vector<double>> eigenvalues_Q;
     std::vector<double>              chi_Q;
     std::vector<std::vector<int>>    n_Q;
     //    ########################################
   private:
     void createBasis(double JKondo, double spinS, double localDelta) {
       // TODO(saurabh): check spinS is multiple of 2.
       auto spinDim = static_cast<size_t>(2. * spinS + 1.);
       std::cout << "spinDim: " << spinDim << std::endl;
       qmatrix<double> spinSz(spinDim, spinDim, 0);
       qmatrix<double> spinSplus(spinDim, spinDim, 0);
       // Spin Z
       for (size_t i = 0; i < spinDim; i++) {
         double m     = -spinS + i * 1.0;
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
       std::cout << "wSpinx: " << wSpinx << std::endl;
       std::cout << "wSpiny: " << wSpiny << std::endl;
       std::cout << "wSpinz: " << wSpinz << std::endl;
       // ############################################
       auto Hamiltonian =
           (spinSz.krDot(wSpinz) + spinSx.krDot(wSpinx) //
            - spinSy.krDot(wSpiny)) * // Imaginary part is taken care here
               JKondo
           // SC term
           + spinSz.id().krDot(
                 (c_up_dag.dot(c_Down_dag) +
                  c_Down_dag.cTranspose().dot(c_Down_dag.cTranspose()))) *
                 localDelta;
       // End
       // std::cout << "spinSz: " << spinSz << std::endl;
       spinSz = qmatrix(spinSz.krDot(wSpinx.id()));
       // std::cout << "spinSz: " << spinSz << std::endl;
       wSpinz  = spinSx.id().krDot(wSpinz);
       wNtotal = spinSx.id().krDot(wNtotal);
       std::cout << "wSpinz: " << wSpinz * 2 << std::endl;
       c_up_dag   = spinSx.id().krDot(c_up_dag);
       c_Down_dag = spinSx.id().krDot(c_Down_dag);
       fnParticle.clear();
       // Number of particles
       fnParticle.push_back(((spinSz.id() + spinSz * 2.) * 0.5).getdiagonal());
       fnParticle.push_back(((spinSz.id() - spinSz * 2.) * 0.5).getdiagonal());
       fnParticle.push_back(((wNtotal + wSpinz * 2.) * 0.5).getdiagonal());
       fnParticle.push_back(((wNtotal - wSpinz * 2.) * 0.5).getdiagonal());
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
       std::vector<std::vector<size_t>> SpinSz = {tm1, tm2};
       create_QuantumSpinOnly(SpinSz);
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
       std::cout << "h_blocked: " << h_blocked << std::endl;
       std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
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
   };
