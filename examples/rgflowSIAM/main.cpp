#include "models/spinhalf.hpp"
#include "nrgcore/nrgcore.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
double hopping(int site, double lmbd) {
  // renormalized Hopping Parameters
  return 0.5 * (1.0 + 1.0 / lmbd) * (1. - std::pow(lmbd, -site - 1)) /
         std::sqrt((1.0 - std::pow(lmbd, -2. * site - 1)) *
                   (1.0 - std::pow(lmbd, -2. * site - 3)));
}
int main() {
  // file where outputput will be wriiten
  h5stream::h5stream rfile("resultSIAM.h5");
  const double       LAMBDA = 2.50; // NRG Lambda
  timer              mtime("Total time : ");
  int                nMax{41}; // Number of NRG iteration
  double             U_int = 0.20;
  double             GAMMA = 0.0100;
  double fc  = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  double V   = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  double eps = -U_int * 0.50;
  // Impurity Model class
  spinhalf impurity(eps, U_int);
  //  Bath site class with zero onsite and Coloumb Energy
  spinhalf bathModel(0, 0);
  // Build the NRG class and add the First site.
  nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
  siam.set_parameters(1024);       // set max number of states to be kept
  siam.add_bath_site({V, V}, 1.0); // first site.
  siam.update_internal_state();
  // Iteratively Add all the site.
  std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
  for (int in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    siam.add_bath_site({hopping(in, LAMBDA), hopping(in, LAMBDA)}, rescale);
    siam.update_internal_state();
    // Save the eigenvalue of the current iteration
    rfile.write(siam.all_eigenvalue, "Eigenvalues" + std::to_string(in));
  }
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
