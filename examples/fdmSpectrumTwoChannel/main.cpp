#include "dynamics/fdmSpectrum.hpp"
#include "models/fermionBasis.hpp"
#include "models/spinhalf.hpp"
#include "models/twoChannel.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "nrgcore/sysOperator.hpp"
#include <iostream>
#include <string>
double LAMBDA = 3.00; // Dont do this
double hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                   (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  timer mtime("Total time : ");
  // Parameter ############################################
  size_t nMax{40};         // Number of NRG iteration
  size_t minIterations{0}; // Number of NRG iteration after which higher energy
  // states are discarded
  double bandWidth   = 1.0; // Overall energy scale.
  double U_int       = 0.250 * bandWidth;
  double GAMMA       = 0.0100 * bandWidth;
  double GAMMA_left  = 1.0 * GAMMA;
  double GAMMA_right = 1. * GAMMA;
  double fc          = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  double V_left      = std::sqrt(2.0 * fc * GAMMA_left / std::acos(-1.));
  double V_right     = std::sqrt(2.0 * fc * GAMMA_right / std::acos(-1.));
  double eps         = -U_int * 0.50;
  // Parameter ############################################
  h5stream::h5stream rfile(std::string("resultsTcRT") +
                           std::to_string(std::log(GAMMA_right / GAMMA_left)) +
                           ".h5");
  spinhalf           impurity(eps, U_int);
  // Enlarge the no of fermion on the impurity
  // i.e., just copy the spin up and down operator
  impurity.f_dag_operator.push_back(impurity.f_dag_operator[0]);
  impurity.f_dag_operator.push_back(impurity.f_dag_operator[1]);
  twoChannel bath_model;
  // spinhalf bath_model(0, 0); // For Siam
  // ------------------------Create operators
  // Ask the model to provide it .
  nrgcore siamTc(impurity, bath_model);
  siamTc.set_parameters(1024); // set max number of states to be kept
  // first site. This is consistentent with
  // Bulla's RMP
  std::vector<qOperator> dUpDownDagger; // Up and Down particle number
  {
    fermionBasis spinhalfBasis(2); // Number of fermion channels/spins
    auto         f_dag_raw = spinhalfBasis.get_raw_f_dag_operator();
    // set f_operator
    // Total Number of particle
    auto ntotal = f_dag_raw[0].dot(f_dag_raw[0].cTranspose());
    for (size_t i = 1; i < f_dag_raw.size(); i++) {
      ntotal = ntotal + f_dag_raw[i].dot(f_dag_raw[i].cTranspose());
    }
    auto chi = ntotal;
    for (size_t i = 0; i < chi.getrow(); i++) {
      chi(i, i) = std::pow(-1, ntotal(i, i));
    }
    // d. \Chi
    for (auto &aa : f_dag_raw) {
      aa = aa.dot(chi);
    }
    dUpDownDagger = spinhalfBasis.get_block_operators(f_dag_raw);
  }
  NrgData rawData(&siamTc);
  // Hopping element for each operator
  auto thop = std::vector<double>{V_left, V_left, V_right, V_right};
  siamTc.add_bath_site(thop, 1.0);
  update_system_operator(&siamTc, dUpDownDagger);
  siamTc.update_internal_state();
  for (size_t in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    thop = {hopping(in), hopping(in), hopping(in), hopping(in)};
    siamTc.add_bath_site(thop, rescale);
    update_system_operator(&siamTc, dUpDownDagger);
    if (siamTc.checkHigherEnergyDiscarded()) {
      rawData.saveCurrentData();
      rawData.saveqOperator(&dUpDownDagger, "nparticle");
    } else {
      minIterations++;
    }
    // Update all the previious bath operator
    siamTc.update_internal_state();
    // Bulla's RMP
    // std::cout << "Eigenvalues: " << siamTc.all_eigenvalue << std::endl;
    // Save the eigenvalue of the current iteration
    rfile.write(siamTc.all_eigenvalue, "Eigenvalues" + std::to_string(in));
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Backword iteration starts now.
  // BOOM:: BOOM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%% Backward Iteration %%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // Load data
  // Start ofd the fdmSpectrum
  fdmSpectrum fdmTc(&siamTc);
  fdmTc.setOperator(&dUpDownDagger);
  for (size_t in = nMax; in > minIterations; in--) {
    // We dont need load the nrg data for the
    // last iteration
    std::cout << "--Started Backward iteration: " << in
              << " Mtr: " << minIterations << std::endl;
    siamTc.nrg_iterations_cnt--;
    rawData.loadCurrentData();
    rawData.loadqOperator(&dUpDownDagger, "nparticle");
    fdmTc.calcSpectrum(std::pow(LAMBDA, -(in - 1.0) / 2.0));
  }
  fdmTc.saveFinalData(&rfile);
  //
  //
  // Dont do this
  rfile.close();
  return 0;
}