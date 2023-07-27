
.. _program_listing_file_nrgcore_include_models_spinhalfSC.hpp:

Program Listing for File spinhalfSC.hpp
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_models_spinhalfSC.hpp>` (``nrgcore/include/models/spinhalfSC.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "models/fermionBasis.hpp"
   #include "nrgcore/qOperator.hpp"
   #include "nrgcore/qsymmetry.hpp"
   #include "utils/qmatrix.hpp"
   #include <algorithm>
   #include <cmath>
   #include <iostream>
   #include <iterator>
   #include <map>
   #include <numeric>
   #include <optional>
   #include <string>
   #include <tuple>
   #include <vector>
   class spinhalfSC {
     std::vector<qmatrix<double>> f_dag_raw;
     fermionBasis                 localSCbabsis;
   
   public:
     spinhalfSC(double teps, double tUint, double tmag) // NOLINT
         : localSCbabsis(2, fermionBasis::spinOnly) {
       double Uint           = tUint;
       double epsilon_d      = teps;
       double magnetic_field = tmag;
       //
       // create_basis();
       // set_foperator();
       // set_chi_Q();
       // create_hamiltonian();
       // localSCbabsis = fermionBasis(2, fermionBasis::spinOnly);
       // Number of fermion channels/spins
       f_dag_operator = localSCbabsis.get_f_dag_operator();
       f_dag_raw      = localSCbabsis.get_raw_f_dag_operator();
       // n_Q
       n_Q = localSCbabsis.get_unique_Qnumbers();
       // set chi_Q
       chi_Q.clear();
       for (auto ai : n_Q) {
         double t_charge = std::accumulate(ai.begin(), ai.end(), 0);
         chi_Q.push_back(std::pow(-1., t_charge));
       }
       // std::cout << "chi_Q: " << chi_Q << std::endl;
       // set up the Hamiltonian
       auto n_up   = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
       auto n_down = f_dag_raw[1].dot(f_dag_raw[1].cTranspose());
       // std::cout << "n_up" << n_up << "n_down" << n_down;
       auto H = (n_up + n_down) * epsilon_d        // onsite energy
                + (n_up - n_down) * magnetic_field // Magnetic Field
                + (n_up.dot(n_down)) * Uint;       // Columb Energy
       // get the hamiltonian  in the blocked structure
       auto h_blocked = localSCbabsis.get_block_Hamiltonian(H);
       //  h_blocked.display();
       // Diagonalize the hamilton
       eigenvalues_Q.clear();
       eigenvalues_Q.resize(n_Q.size(), {});
       for (size_t i = 0; i < n_Q.size(); i++) {
         eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
       }
       // TODO(sp): rotate the f operator
       // ####################################################################
       f_dag_operator =
           localSCbabsis.get_block_operators({f_dag_raw[0], f_dag_raw[1]});
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
       // End of the constructor
     }
     //
     void addSCDelta(double delta) {
       auto H = (f_dag_raw[0].cTranspose().dot(f_dag_raw[1].cTranspose()) +
                 f_dag_raw[1].dot(f_dag_raw[0])) *
                delta;
       auto h_blocked = localSCbabsis.get_block_Hamiltonian(H);
       //  h_blocked.display();
       // Diagonalize the hamilton
       eigenvalues_Q.clear();
       eigenvalues_Q.resize(n_Q.size(), {});
       for (size_t i = 0; i < n_Q.size(); i++) {
         eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
       }
       // TODO(sp): rotate the f operator
       // ####################################################################
       f_dag_operator =
           localSCbabsis.get_block_operators({f_dag_raw[0], f_dag_raw[1]});
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
     [[nodiscard]] std::vector<std::vector<int>> get_basis() const {
       return n_Q;
     }
     [[nodiscard]] std::vector<std::vector<double>> get_eigenvaluesQ() const {
       return eigenvalues_Q;
     }
     [[nodiscard]] std::vector<double> get_chi_Q() const {
       return chi_Q;
     }
     // protected:
     // functions
     void set_chi_Q();
     //
     std::vector<qOperator>           f_dag_operator;
     std::vector<std::vector<double>> eigenvalues_Q;
     std::vector<double>              chi_Q;
     std::vector<std::vector<int>>    n_Q;
   };
