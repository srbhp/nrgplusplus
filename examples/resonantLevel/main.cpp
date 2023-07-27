/* This is the main function
 * Which does all the NRG calculations,
 * It will also call the function
 * for the backward iteration
 */
#include "dynamics/fdmSpectrum.hpp"
#include "models/resonantTwoLead.hpp"
#include "models/spinlessTwoLead.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "nrgcore/sysOperator.hpp"
#include "openChain/currentShellCouplingOpenChain.hpp"
#include "openChain/staticExpectationOpenChain.hpp"
#include <iostream>
#include <string>
double LAMBDA = 3.00; // Dont do this
double hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                   (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
// int main(int argc, char *argv[]) {
int main() {
  std::srand(std::time(0));
  timer mtime("Total time : ");
  // Parameter ############################################
  size_t nMax          = 6; // Number of NRG iteration
  int    minIterations = 0; // Number of NRG iteration after which higher energy
  // states are discarded
  double bandWidth   = 1.0;              // Overall energy scale.
  double U_int       = bandWidth * 0.05; // default 0.25 std::stod(argv[1]) *
  double GAMMA       = 0.01 * bandWidth;
  double GAMMA_left  = GAMMA;
  double GAMMA_right = GAMMA;
  double fc          = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  double V_left      = std::sqrt(2.0 * fc * GAMMA_left / std::acos(-1.));
  double V_right     = std::sqrt(2.0 * fc * GAMMA_right / std::acos(-1.));
  double eps         = -U_int * 0.50;
  // Parameter ############################################
  h5stream::h5stream rfile(
      std::string("resultsRLM") + //
                                  // std::to_string(U_int)       //
      +".h5");
  // #########################################
  // Set up impurity and bath
  //  spinhalf impurity(eps, U_int);
  //  // Enlarge the no of fermion on the impurity
  //  // i.e., just copy the spin up and down operator
  //  impurity.f_dag_operator.push_back(impurity.f_dag_operator[0]);
  //  impurity.f_dag_operator.push_back(impurity.f_dag_operator[1]);
  //  twoChannel bath_model;
  // spinhalf bath_model(0, 0); // For Siam
  // resonantTwoLead
  resonantTwoLead  impurity(eps, V_left, V_right, U_int);
  spinnlessTwoLead bath_model; // For Siam
  // ------------------------Create operators
  // Ask the model to provide it .
  nrgcore siamTc(impurity, bath_model);
  siamTc.set_parameters(128); // set max number of states to be kept
  // first site.
  // This is consistentent with
  // Bulla's RMP
  std::vector<qOperator> nUpDownOperator{impurity.impurityNparticle};
  currentOpenChain       fullOpenCalc(&siamTc,          // NRG objeect
                                      &nUpDownOperator, // NRG local operator
                                      LAMBDA);
  // Hopping element for each operator
  std::vector<double> thop;
  for (size_t in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    thop = {hopping(in), hopping(in), hopping(in), hopping(in)};
    siamTc.add_bath_site(thop, rescale);
    update_system_operator(&siamTc, nUpDownOperator);
    if (!siamTc.checkHigherEnergyDiscarded()) {
      minIterations++;
    }
    // bath operators are always saved
    fullOpenCalc.saveFullNRGState();
    // Update all the previous bath operator
    // Every site has bath operators
    // openChainTC.saveBathOperators();
    //
    // End of openChainData
    siamTc.update_internal_state();
    // Bulla's RMP
    // std::cout << "Eigenvalues: " << siamTc.all_eigenvalue << std::endl;
    // Save the eigenvalue of the current iteration
    rfile.write(siamTc.all_eigenvalue, "Eigenvalues" + std::to_string(in));
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Backward iteration starts now.
  // BOOM:: BOOM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%% Backward Iteration %%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // Start ofd the fdmSpectrum
  // Voltage  ############################################
  std::vector<double> voltageArray;
  double              vMax    = 1; // Max Voltage can be 2* bandwidth
  int                 tPoints = 20;
  //   ####### Linear scale
  for (int i = -tPoints + 1; i < tPoints; i++) {
    voltageArray.push_back(i * vMax / (1. * tPoints));
  }
  //  double vMin    = 1e-6;
  //  double delTime = (tPoints - 1.0) / (std::log(vMax / vMin));
  //  // Log scale
  //  for (int i = -tPoints + 1; i < tPoints; i++) {
  //    if (i == 0) { // Mid point of the voltage is zero
  //      voltageArray.push_back(0);
  //    } else {
  //      voltageArray.push_back(sign(i) * vMin *
  //                             std::exp(std::fabs(i) * 1. / delTime));
  //    }
  //  }
  // #####################################################
  fullOpenCalc.setVolltageArray(voltageArray);
  fullOpenCalc.openChainFullCalculation();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%Second Backward Iteration %%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  //  for (int in = nMax; in > minIterations; in--) {
  //    // We dont need load the nrg data for the
  //    // last iteration
  //    std::cout << "--Started Backward iteration: " << in
  //              << " Mtr: " << minIterations << std::endl;
  //    siamTc.nrg_iterations_cnt--;
  //    // fullOpenCalc.calcSpectrum(std::pow(LAMBDA, -(in - 1.0) / 2.0));
  //    fullOpenCalc.openChainFullCalculation(std::pow(LAMBDA, -(in - 1.0)
  //    / 2.0)); std::cout << "##########################" << std::endl;
  //  }
  fullOpenCalc.saveOpenChainFinalData(&rfile);
  //
  // Dont do this
  rfile.close();
  return 0;
}
