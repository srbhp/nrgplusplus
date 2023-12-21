#include "models/spinhalf.hpp"
#include "models/spinhalfKondo.hpp"
#include "nrgcore/nrgcore.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
const double LAMBDA = 2.0; // Dont do this
double       hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                         (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  h5stream::h5stream rfile("results.h5");
  timer              mtime("Total time : ");
  int                nMax{41}; // Number of NRG iteration
  // double             GAMMA = 0.0100;
  // double fc = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  // double V = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  double        Jkondo = 0.250;
  spinhalfKondo impurity(Jkondo, 0.5);
  spinhalf      bathModel(0, 0); // set parameters
  std::cout << "Done !" << std::endl;
  nrgcore kondoModel(impurity, bathModel);
  std::cout << "f_dag_operators: " << impurity.f_dag_operator.size()
            << std::endl;
  kondoModel.set_parameters(1024); // set max number of states to be kept
  std::cout << "Eigenvalues: " << kondoModel.all_eigenvalue << std::endl;
  for (int in = 0; in < nMax; in++) {
    // Rescaling Factor is different for the Kondo Model
    double rescale = std::sqrt(LAMBDA);
    kondoModel.add_bath_site({hopping(in), hopping(in)}, rescale);
    // Update operators if needed here. then ->
    kondoModel.update_internal_state();
    // std::cout << "Eigenvalues: " << kondoModel.all_eigenvalue << std::endl;
    rfile.write(
        kondoModel.all_eigenvalue,
        "Eigenvalues" +
            std::to_string(in)); // Save the eigenvalue of the current iteration
  }
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
