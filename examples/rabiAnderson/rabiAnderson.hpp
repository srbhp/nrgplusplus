
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
class rabiAnderson : public fermionBasis {
    // See PHYSICAL REVIEW B 101, 085110 (2020)
    /** This class is for a single orbital with spin up and down
     * f operator. SIAM can be made entirely from this class.
     *
     *
     */
public:
    /**
     * @brief [TODO:description]
     *
     * @param UCoulumb [TODO:description]
     * @param gamma [TODO:description]
     * @param omega [TODO:description]
     */
    rabiAnderson(const std::map<std::string, double> &params) {
        dof = 1 + // for the trion State (up)
              2 + // for the electron state
              2;  // for the wilson site
        // Dimension of the system
        nstates = std::pow(2, dof);
        // Dimension of the system
        std::cout << "nstates: " << nstates << "\n"
                  << "dof: " << dof << "\n";
        createBasis(params); // create the basis in nstates x nstates
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
    void createBasis(const std::map<std::string, double> &params) {
        double gammaZero    = params.at("gammaZero");
        double epsilonTrion = params.at("epsilonTrion");
        double omega        = params.at("Omega");
        double epsilonD     = params.at("epsilonImpurity");
        double UCoulumb     = params.at("UColoumbImpurity");
        double UTrion       = params.at("UColoumbTrion");
        // if (std::fabs(omega) < 1e-20) {
        //     UTrion = 0;
        // }
        //
        createFermionBasis();
        std::cout << "FermionBasis Size" << fermionOprMat.size() << "\n";
        // trion particle
        auto trParticle = fermionOprMat[4].dot(fermionOprMat[4].cTranspose());
        auto nDown      = fermionOprMat[3].dot(fermionOprMat[3].cTranspose());
        auto nUp        = fermionOprMat[2].dot(fermionOprMat[2].cTranspose());
        // ########################################################
        auto Hamiltonian =
            // Impurity Onsite Energy
            (nDown + nUp) * epsilonD +
            nDown.dot(nUp) * UCoulumb
            // Trion Onsite Energy
            + (trParticle)*epsilonTrion +               // Trion onsite
            (nDown + nUp).dot(trParticle) * (-UTrion) + //
            //          Wilson Site coupling
            (fermionOprMat[0].dot(fermionOprMat[2].cTranspose()) +
             fermionOprMat[2].dot(fermionOprMat[0].cTranspose()) +
             fermionOprMat[1].dot(fermionOprMat[3].cTranspose()) +
             fermionOprMat[3].dot(fermionOprMat[1].cTranspose())) *
            gammaZero
            // Trion coupling
            + ((fermionOprMat[3].cTranspose().dot(fermionOprMat[4].cTranspose()) +
                fermionOprMat[4].dot(fermionOprMat[3])) *
               omega);
        // End
        //
        //
        create_QuantumNspinCharge();
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
        //    std::cout << "h_blocked: " << h_blocked << std::endl;
        //    std::cout << "Hamiltonian: " << Hamiltonian << std::endl;
        // Diagonalize the hamilton
        eigenvalues_Q.clear();
        eigenvalues_Q.resize(n_Q.size(), {});
        for (size_t i = 0; i < n_Q.size(); i++) {
            eigenvalues_Q[i] = (h_blocked.get(i, i)).value()->diag();
        }
        //    std::cout << "Eigenvalues: " << eigenvalues_Q << std::endl;
        // TODO: rotate the f operator
        //####################################################################
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
