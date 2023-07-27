
.. _program_listing_file_nrgcore_include_dynamics_fdmSpectrum.hpp:

Program Listing for File fdmSpectrum.hpp
========================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_dynamics_fdmSpectrum.hpp>` (``nrgcore/include/dynamics/fdmSpectrum.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "nrgcore/qOperator.hpp"
   #include "utils/qmatrix.hpp"
   #include "utils/timer.hpp"
   #include <algorithm> // std::min_element
   #include <cmath>
   #include <iostream>
   #include <iterator>
   #include <map>
   #include <numeric>
   #include <optional>
   #include <tuple>
   #include <vector>
   template <typename nrgcore_type> // nrgcore_type is a type of
   // We don't want to inherit the nrgcore_type.
   // Just want to use few of its objects as a pointer
   // We don't want to reference to any specific model here.
   //
   // Warning: this class DOES NOT take care of the fermion sign.
   // You need to provide the operator with  sign
   // TODO(sp): : kBt != 0 and a = b^\dag
   class fdmSpectrum {
     nrgcore_type *nrg_object;
     double        kBT{0}; // Temperature of nrg system i.e., in FDM formalism
   public:
     // fdmSpectrum() {}
     explicit fdmSpectrum(nrgcore_type *t_nrg_object // nrgcore_type
     ) {
       setup(t_nrg_object);
     }
     void setup(nrgcore_type *t_nrg_object) {
       lastiteration           = true;
       nrg_object              = t_nrg_object;
       globalGroundStateEnergy = 0; //= nrg_object->all_eigenvalue[0];
     }
     void calcSpectrum(double energyScale) {
       energyRescale = energyScale;
       // Clear the operator
       setCurrentIndex();
       setRhoZero();
       rhoDotOperators();
       setReduceDensityMatrix();
       lastiteration = false;
     }
     void setReduceDensityMatrix() {
       // timer t1("setReduceDensityMatrix");
       // Set reducedRho
       reducedRho.clear();
       for (size_t i = 0; i < nrg_object->pre_sysmQ.size(); i++) {
         reducedRho.push_back( //
             qmatrix<>(nrg_object->eigenvaluesQ_kept_indices[i].size(),
                       nrg_object->eigenvaluesQ_kept_indices[i].size(), 0));
       }
       // Rotate the eigen basis
       // #pragma omp parallel for
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         // U. rhoZero . U.T : TODO: check
         rhoZero[i] = nrg_object->current_hamiltonQ[i].dot(
             rhoZero[i].dot(nrg_object->current_hamiltonQ[i].cTranspose()));
       }
       // Set reducedRho
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         size_t kidx = 0;
         for (auto kindex : nrg_object->coupled_nQ_index[i]) {
           auto ii = kindex / nrg_object->nq_bath.size(); // impurity nqi index
           auto bb = kindex % nrg_object->nq_bath.size(); // bath nqi index
           // create previous bath id matrix
           // TODO(sp): Add parallelization here
           for (size_t it : nrg_object->eigenvaluesQ_kept_indices[ii]) {
             for (size_t it_p : nrg_object->eigenvaluesQ_kept_indices[ii]) {
               double aa{0};
               for (size_t il = 0; il < nrg_object->bath_eigenvaluesQ[bb].size();
                    il++) {
                 aa += rhoZero[i].at(
                     kidx + it +
                         nrg_object->eigenvaluesQ_kept_indices[ii].size() * il,
                     kidx + it_p +
                         nrg_object->eigenvaluesQ_kept_indices[ii].size() * il);
               }
               reducedRho[ii].at(it, it_p) += aa;
             }
           }
           kidx += nrg_object->eigenvaluesQ_kept_indices[ii].size() *
                   nrg_object->bath_eigenvaluesQ[bb].size();
         }
         // End of matrix generation.
       }
     }
     void setLocalPartitionFunction() {
       localGroundStateEnergy = 0;
       localPartitionFunction = 0; // Ground state degenarecy
       for (auto &aa : nrg_object->eigenvaluesQ) {
         if (aa.size() != 0) {
           double result          = *std::min_element(aa.begin(), aa.end());
           localGroundStateEnergy = std::min(localGroundStateEnergy, result);
         }
       }
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
           double energy =
               std::fabs(nrg_object->eigenvaluesQ[i][ie] - localGroundStateEnergy);
           if (energy < energyErrorBar) {
             localPartitionFunction += 1.;
           }
         }
       } //
       // std::cout << "localGroundStateEnergy" << localGroundStateEnergy
       //           << " localPartitionFunction: " << localPartitionFunction
       //           << std::endl;
       BoltzmannFactor = nrg_object->eigenvaluesQ;
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
           double energy =
               std::fabs(nrg_object->eigenvaluesQ[i][ie] - localGroundStateEnergy);
           if (energy < energyErrorBar) {
             BoltzmannFactor[i][ie] = 1. / localPartitionFunction;
           } else {
             BoltzmannFactor[i][ie] = 0;
           }
         }
       } //
     }
     std::vector<std::vector<double>> BoltzmannFactor;
     void                             setRhoZero() {
       if (lastiteration) {
         setLocalPartitionFunction();
       }
       double rhoTrace{0};
       // Calc partition function of the shell ::
       vecPartitions.push_back(localPartitionFunction);
       rhoZero.clear();
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         size_t kpdim = nrg_object->eigenvaluesQ[i].size();
         // std::cout << "--------------------------";
         // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
         qmatrix<> tmat(kpdim, kpdim, 0);
         // Discarded states
         if (lastiteration) { // Only for the laast iterationT= 0
           for (size_t ie = 0; ie < nrg_object->eigenvaluesQ[i].size(); ie++) {
             tmat(ie, ie) = BoltzmannFactor[i][ie];
           }
         }
         // Kept states
         if (!reducedRho.empty()) { // Condition for the last Wilson site
           //  Just Override the kept states
           for (auto ik : currentKeptIndex[i]) {
             for (auto ikp : currentKeptIndex[i]) {
               tmat(ik, ikp) = reducedRho[i](ik, ikp);
             }
           }
         }
         rhoTrace += tmat.trace();
         // Save the matrix
         rhoZero.push_back(tmat);
       }
       std::cout << "rhoTrace: " << rhoTrace << std::endl;
       // move the operator
     }                        // End of  update_system_operatorQ
     void rhoDotOperators() { // NOLINT
       double specSum = 0.0;
       //
       for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
         for (size_t j = 0; j < nrg_object->current_sysmQ.size(); j++) {
           size_t kpdim   = nrg_object->eigenvaluesQ[i].size();
           size_t kpdim_p = nrg_object->eigenvaluesQ[j].size();
           // std::cout << "--------------------------";
           // std::cout << "idx" << idx << " idx_p" << idx_p << std::endl;
           // TODO(sp): This
           for (size_t ip = 0; ip < bOperator->size(); ip++) {
             auto sys_opr_opt = (*bOperator)[ip].get(i, j);
             if (sys_opr_opt) {
               auto     *sys_opr = sys_opr_opt.value();
               qmatrix<> aMatrix;
               if (aOperator == nullptr) {
                 aMatrix = sys_opr_opt.value()->cTranspose();
               } // Else this should be just boperator
               // set kept-kept part operator
               for (auto iv : currentKeptIndex[i]) {
                 for (auto iv_p : currentKeptIndex[j]) {
                   aMatrix.at(iv_p, iv) = 0;
                 }
               }
               auto rhoA = sys_opr->dot(rhoZero[j]);
               auto ARho = rhoZero[i].dot(*sys_opr);
               for (size_t iv = 0; iv < kpdim; iv++) {
                 for (size_t iv_p = 0; iv_p < kpdim_p; iv_p++) {
                   double aa{0};
                   double bbv{0};
                   // for (size_t iv_pp = 0; iv_pp < kpdim_p; iv_pp++) {
                   //  aa += sys_opr->at(iv, iv_pp) * rhoZero[j](iv_pp, iv_p);
                   //}
                   // for (size_t iv_pp = 0; iv_pp < kpdim; iv_pp++) {
                   //  bbv += rhoZero[i](iv, iv_pp) * sys_opr->at(iv_pp, iv_p);
                   //}
                   aa          = rhoA(iv, iv_p) * aMatrix.at(iv_p, iv);
                   bbv         = ARho(iv, iv_p) * aMatrix.at(iv_p, iv);
                   int tmindex = int(
                       std::log(1 +
                                (std::fabs(energyRescale *
                                           (nrg_object->eigenvaluesQ[i][iv] -
                                            nrg_object->eigenvaluesQ[j][iv_p])) /
                                 minEnergy)) *
                       delE);
                   if (tmindex < energyPts && tmindex >= 0) {
                     if (std::fabs(aa) > spWeightErrorBar) {
                       positiveWeight[ip][tmindex] += std::fabs(aa);
                       specSum += std::fabs(aa);
                       // std::cout << "Raw:  " << tmindex << " " << std::fabs(aa)
                       //         << " " << std::fabs(bbv) << std::endl;
                     }
                     if (std::fabs(bbv) > spWeightErrorBar) {
                       negativeWeight[ip][tmindex] += std::fabs(bbv);
                       specSum += std::fabs(bbv);
                     }
                   }
                 }
               }
             }
           }
           // all the matrices are set
         }
       }
       std::cout << nrg_object->nrg_iterations_cnt << "specSum: " << specSum
                 << " Scale: " << energyRescale << std::endl;
       // set the matrix elements
       // end of lm loop
       // End of matrix generation.
       // Rotate the c operator in the eigen basis
     }
     //
     void setTemperature(double at) { kBT = at; }
     void setCurrentIndex() {
       if (currentKeptIndex.empty()) {
         for (size_t i = 0; i < nrg_object->current_sysmQ.size(); i++) {
           currentKeptIndex.emplace_back();
         }
       } else {
         currentKeptIndex = previoudKeptIndex;
       }
       previoudKeptIndex = nrg_object->eigenvaluesQ_kept_indices;
     }
     void setOperator(std::vector<qOperator> *bopr,
                      std::vector<qOperator> *aopr = nullptr) {
       aOperator = aopr;
       bOperator = bopr;
       for (size_t i = 0; i < bOperator->size(); i++) {
         positiveWeight.emplace_back(energyPts, 0);
         negativeWeight.emplace_back(energyPts, 0);
       }
     }
     template <typename filetype> void saveFinalData(filetype *pfile) {
       std::vector<double> energyPoints(energyPts, 0);
       for (int i = 0; i < energyPts; i++) {
         energyPoints[i] = (minEnergy * std::exp(i * 1. / delE));
       }
       std::string hstr = "GreenFn";
       pfile->write(energyPoints, hstr + "EnergyPoints");
       pfile->write(positiveWeight, hstr + "PositiveWeight");
       pfile->write(negativeWeight, hstr + "NegativeWeight");
       for (size_t i = 0; i < bOperator->size(); i++) {
         positiveWeight[i].clear();
         negativeWeight[i].clear();
       }
       // pfile->write(vecPartitions, "vecPartitions");
     }
   
   private:
     double localGroundStateEnergy{0};
     double localPartitionFunction{0}; // Ground state degenarecy
     std::vector<std::vector<size_t>> previoudKeptIndex;
     double spWeightErrorBar{1e-20}; // We dont care for the lower value
     double globalGroundStateEnergy{0};
     std::vector<qmatrix<>> reducedRho;
     std::vector<double>    vecPartitions;
     // rho.B and B.rho
     std::vector<std::vector<double>> positiveWeight;
     std::vector<std::vector<double>> negativeWeight;
     // value
     int    energyPts = 100000;
     double maxEnergy = 10; // One decade more
     double minEnergy = 1e-10;
     double delE      = (energyPts - 1.0) / (std::log(maxEnergy / minEnergy));
   
   public: // Give access for openchain class
     std::vector<std::vector<size_t>> currentKeptIndex;
     std::vector<qmatrix<>>           rhoZero;
     std::vector<qOperator>          *aOperator{};
     std::vector<qOperator>          *bOperator{};
     bool                             lastiteration{true};
     double                           energyRescale{1};
     double                           energyErrorBar{1e-5};
   };
