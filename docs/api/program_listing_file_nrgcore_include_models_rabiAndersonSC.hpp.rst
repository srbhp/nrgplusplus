
.. _program_listing_file_nrgcore_include_models_rabiAndersonSC.hpp:

Program Listing for File rabiAndersonSC.hpp
===========================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_models_rabiAndersonSC.hpp>` (``nrgcore/include/models/rabiAndersonSC.hpp``)

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
   class rabiAndersonSC : public fermionBasis {
     // This class is the same as the rabi-Anderson Model.
     // Except the fact  that We now have Superconductor
     // in the system. So Spin is now only conserved quatity
     // and charge isn't conserved any more.
     // CHECK: The class for the Wilson chaini has to be consistent.
   public:
     explicit rabiAndersonSC(const std::map<std::string, double> &params) {
       createBasis(params); // create the basis in nstates x nstates
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
     void createBasis(const std::map<std::string, double> &params) {
       double gammaZero    = params.at("gammaZero");
       double epsilonTrion = params.at("epsilonTrion");
       double omega        = params.at("Omega");
       double epsilonD     = params.at("epsilonImpurity");
       double UCoulumb     = params.at("UColoumbImpurity");
       double UTrion       = params.at("UColoumbTrion");
       std::cout << "=======Models Parameter=========" << std::endl;
       for (const auto &para : params) {
         std::cout << para.first << ": " << para.second << "\n";
       }
       std::cout << "================================" << std::endl;
       createFermionBasis();
       std::cout << "FermionBasis Size" << fermionOprMat.size() << "\n";
       // trion particle
       auto trionHole = fermionOprMat[4].cTranspose().dot(fermionOprMat[4]);
       auto nDown     = fermionOprMat[3].dot(fermionOprMat[3].cTranspose());
       auto nUp       = fermionOprMat[2].dot(fermionOprMat[2].cTranspose());
       // ########################################################
       auto Hamiltonian =
           // Impurity Onsite Energy
           (nDown + nUp) * epsilonD +
           nDown.dot(nUp) * UCoulumb
           // Trion Onsite Energy
           + (trionHole)*epsilonTrion +               // Trion onsite
           (nDown + nUp).dot(trionHole) * (-UTrion) + //
           //          Wilson Site coupling
           (fermionOprMat[0].dot(fermionOprMat[2].cTranspose()) +
            fermionOprMat[2].dot(fermionOprMat[0].cTranspose()) +
            fermionOprMat[1].dot(fermionOprMat[3].cTranspose()) +
            fermionOprMat[3].dot(fermionOprMat[1].cTranspose())) *
               gammaZero
           // Trion coupling
           + (fermionOprMat[2].dot(fermionOprMat[4].cTranspose()) +
              fermionOprMat[4].dot(fermionOprMat[2].cTranspose())) *
                 omega;
       // End
       //
       //
       // create_QuantumNspinCharge();
       // Main change happen here because of SC.
       // Spin (S_z) only conserved quantity.
       std::vector<std::vector<size_t>> spinIndex({{0, 2, 4}, // UpSpin
                                                   {1, 3}}    // Down Spin
       );
       create_QuantumSpinOnly(spinIndex);
       create_Block_structure();
       // Only S_z is quatized
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
       //    std::cout << "h_blocked: " << h_blocked << std::endl;
       //    std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
       // Diagonalize the hamilton
       eigenvalues_Q.clear();
       eigenvalues_Q.resize(n_Q.size(), {});
       for (size_t i = 0; i < n_Q.size(); i++) {
         eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
       }
       std::cout << "Eigenvalues: " << eigenvalues_Q << std::endl;
       // TODO: rotate the f operator
       // ####################################################################
       f_dag_operator = get_block_operators({fermionOprMat[0], fermionOprMat[1]});
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
