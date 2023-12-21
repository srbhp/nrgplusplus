#include "dynamics/fdmSpectrum.hpp"
#include "models/resonantTwoLead.hpp"
#include "models/spinlessTwoLead.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "nrgcore/sysOperator.hpp"
#include <iostream>
#include <string>
const double LAMBDA = 3.00; // Dont do this
double       hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                         (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
// int main(int argc, char *argv[]) {
int main() {
  std::srand(std::time(nullptr));
  timer mtime("Total time : ");
  // Parameter ############################################
  size_t nMax = 6; // Number of NRG iteration
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
  // currentOpenChain       fullOpenCalc(&siamTc,          // NRG objeect
  //                                     &nUpDownOperator, // NRG local operator
  //                                     LAMBDA);
  // // Hopping element for each operator
  std::vector<double> thop;
  for (size_t in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    thop = {hopping(in), hopping(in), hopping(in), hopping(in)};
    siamTc.add_bath_site(thop, rescale);
    update_system_operator(&siamTc, nUpDownOperator);
    siamTc.update_internal_state();
    // Bulla's RMP
    // std::cout << "Eigenvalues: " << siamTc.all_eigenvalue << std::endl;
    // Save the eigenvalue of the current iteration
    rfile.write(siamTc.all_eigenvalue, "Eigenvalues" + std::to_string(in));
  }
  rfile.close();
  return 0;
}
