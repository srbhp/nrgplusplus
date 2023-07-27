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
class rabiSingleParticle : public fermionBasis {
    /** This class is for a single orbital with spin up and down
     * f operator. SIAM can be made entirely from this class.
     *
     *
     */
public:
    rabiSingleParticle(double JKondo, double omega) {
        // One additional state for the trion state
        nstates = 3 * 4; // std::pow(2, dof); // Dimension of the system
        dof     = 3      // for spin + trion
                  + 2;       // std::pow(2, dof); // Dimension of the system
        std::cout << "nstates: " << nstates << "\n"
                  << "dof: " << dof << "\n";
        createBasis(JKondo, omega); // create the basis in nstates x nstates
    }
    std::vector<std::vector<int>>    get_basis() {
        return n_Q;
    }
    std::vector<std::vector<double>> get_eigenvaluesQ() {
        return eigenvalues_Q;
    }
    std::vector<double>              get_chi_Q() {
        return chi_Q;
    }
    //
    std::vector<std::vector<double>> eigenvalues_Q;
    std::vector<double>              chi_Q;
    std::vector<std::vector<int>>    n_Q;
    //    ########################################
private:
    void createBasis(double JKondo, double omega) {
        // TODO: check spinS is multiple of 2.
        size_t spinDim = 3;
        std::cout << "spinDim: " << spinDim << std::endl;
        qmatrix<double> spinSz(spinDim, spinDim, 0);
        qmatrix<double> spinSx(spinDim, spinDim, 0);
        qmatrix<double> spinSy(spinDim, spinDim, 0);
        qmatrix<double> hOmega(spinDim, spinDim, 0);
        qmatrix<double> trParticle(spinDim, spinDim, 0); // trion particle Number
        qmatrix<double> sParticle(spinDim, spinDim, 0);  // trion particle Number
        // Set Matrices
        spinSz(1, 1) = -0.5;
        spinSz(2, 2) = 0.5;
        //
        spinSx(1, 2) = 0.5;
        spinSx(2, 1) = 0.5;
        //
        spinSy(1, 2) = 0.5;
        spinSy(2, 1) = -0.5;
        //
        hOmega(0, 1) = omega;
        hOmega(1, 0) = omega;
        //
        trParticle(0, 0) = 1.0;
        //
        sParticle(1, 1) = 1.0;
        sParticle(2, 2) = 1.0;
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
            JKondo +
            hOmega.krDot(wSpinz.id()); // Hamiltonian
        // End
        fnParticle.clear();
        // std::cout << "spinSz: " << spinSz << std::endl;
        spinSz     = qmatrix(spinSz.krDot(wSpinx.id()));
        trParticle = qmatrix(trParticle.krDot(wSpinx.id()));
        sParticle  = qmatrix(sParticle.krDot(wSpinx.id()));
        // std::cout << "spinSz: " << spinSz << std::endl;
        wSpinz  = spinSx.id().krDot(wSpinz);
        wNtotal = spinSx.id().krDot(wNtotal);
        std::cout << "wSpinz: " << wSpinz * 2 << std::endl;
        c_up_dag   = spinSx.id().krDot(c_up_dag);
        c_Down_dag = spinSx.id().krDot(c_Down_dag);
        // Number of particles
        fnParticle.push_back(trParticle.getdiagonal());
        fnParticle.push_back(((sParticle + spinSz * 2.) * 0.5).getdiagonal());
        fnParticle.push_back(((sParticle - spinSz * 2.) * 0.5).getdiagonal());
        fnParticle.push_back(((wNtotal + wSpinz * 2.) * 0.5).getdiagonal());
        fnParticle.push_back(((wNtotal - wSpinz * 2.) * 0.5).getdiagonal());
        //
        std::cout << "fnParticle: " << fnParticle << std::endl;
        //
        createQNumbers({{1, 3}, {0, 2, 4}});
        create_Block_structure();
        //####################################################################
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
        //####################################################################
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
