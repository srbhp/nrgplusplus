#include "nrgcore/nrgcore.hpp"
#include "nrgcore/spinhalf.hpp"
#include "nrgcore/twoChannel.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
#include <ostream>
double LAMBDA = 2.50; // Dont do this
double hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                   (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  h5stream::h5stream rfile("resultsTc.h5");
  timer mtime("Total time : ");
  size_t nMax{40}; // Number of NRG iteration
  double U_int = 0.20;
  double GAMMA = 0.0100;
  double fc = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  double V = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  std::sqrt(2.0 * GAMMA / std::acos(-1.));
  double eps = -U_int * 0.50;
  spinhalf impurity(eps, U_int);
  // Enlarge the no of fermion on the impurity

  impurity.f_dag_operator.push_back(impurity.f_dag_operator[0]);
  impurity.f_dag_operator.push_back(impurity.f_dag_operator[1]);
  twoChannel bath_model;
  // --------------- - --------
  // std::cout << "impurity.f_dag_operator " << std::endl;
  // for (auto aa : impurity.f_dag_operator) {
  //  std::cout << "-------" << std::endl;
  //  aa.display();
  //}
  // std::cout << "impurity.f_dag_operator " << std::endl;
  // for (auto aa : bath_model.f_dag_operator) {
  //  std::cout << "-------" << std::endl;
  //  aa.display();
  //}
  // -----------------------------------------------------NRG
  nrgcore<spinhalf, twoChannel> siam(impurity, bath_model);
  siam.set_parameters(100); // set max number of states to be kept
                            // first site. This is consistentent with
                            // Bulla's RMP
  auto thop =
      std::vector<double>{V, V, V, V}; // Hopping element for each operator

  siam.add_bath_site(thop, 1.0);

  siam.update_internal_state();
  std::cout << "Eigenvalues: " << siam.all_eigenvalue.size() << " |"
            << siam.all_eigenvalue << std::endl;
  for (size_t in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    thop = {hopping(in), hopping(in), hopping(in), hopping(in)};
    siam.add_bath_site(thop, rescale);
    siam.update_internal_state();
    // Bulla's RMP
    // std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
    // Save the eigenvalue of the current iteration
    rfile.write("Eigenvalues" + std::to_string(in), siam.all_eigenvalue);
  }
  rfile.close();
  // create operators that connect to the environment
  return 0;
}