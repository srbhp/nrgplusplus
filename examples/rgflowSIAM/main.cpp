#include "models/spinhalf.hpp"
#include "nrgcore/nrgcore.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
double LAMBDA = 2.50; // Dont do this
double hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                   (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  h5stream::h5stream rfile("resultSIAM.h5");
  timer              mtime("Total time : ");
  int                nMax{41}; // Number of NRG iteration
  double             U_int = 0.20;
  double             GAMMA = 0.0100;
  // double fc = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  // double V = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  double   V   = std::sqrt(2.0 * GAMMA / std::acos(-1.));
  double   eps = -U_int * 0.50;
  spinhalf impurity(eps, U_int);
  spinhalf bathModel(0, 0); // set parameters
  std::cout << "Done !" << std::endl;
  nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
  siam.set_parameters(1024);       // set max number of states to be kept
  siam.add_bath_site({V, V}, 1.0); // first site. This is consistentent with
  siam.update_internal_state();
  // Bulla's RMP
  std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
  for (int in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    siam.add_bath_site({hopping(in), hopping(in)}, rescale);
    // std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
    // Update operators if needed here. then ->
    siam.update_internal_state();
    rfile.write(
        siam.all_eigenvalue, // Save the eigenvalue of the current iteration
        "Eigenvalues" + std::to_string(in));
  }
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
