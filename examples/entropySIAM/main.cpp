#include "models/spinhalf.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "static/fdmThermodynamics.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
#include <nrgcore/sysOperator.hpp>
double hopping(int site, double LAMBDA) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                   (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  double LAMBDA = 2.50; // Dont do this
  timer  mtime("Total time : ");
  int    nMax{11};         // Number of NRG iteration
  size_t minIterations{0}; // minimum number of iteration
  double U_int = 0.20;
  double GAMMA = 0.0100;
  // double fc = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  // double V = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  double V   = std::sqrt(2.0 * GAMMA / std::acos(-1.));
  double eps = -U_int * 0.40; // Warning: Epsilon is not half of the U
  // next two are for the calcThermodynamics of the bare band without the
  // impurity
  // spinhalf           impurity(0, 0);
  // h5stream::h5stream  rfile("resultSIAMZero.h5");
  h5stream::h5stream rfile("resultSIAMU" + std::to_string(U_int) + ".h5");
  spinhalf           impurity(eps, U_int);
  auto               doubleOcc = impurity.getDoubleOcc();
  // bath
  spinhalf bathModel(0, 0); // set parameters
  std::cout << "Done !" << std::endl;
  nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
  siam.set_parameters(1024);       // set max number of states to be kept
  siam.add_bath_site({V, V}, 1.0); // first site. This is consistentent with
                                   //
  update_system_operator(&siam, doubleOcc);
  //
  siam.update_internal_state();
  NrgData rawData(&siam);
  // Bulla's RMP
  std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
  for (int in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    siam.add_bath_site({hopping(in, LAMBDA), hopping(in, LAMBDA)}, rescale);
    // std::cout << "Eigenvalues: " << siam.all_eigenvalue << std::endl;
    // Update operators if needed here. then ->
    update_system_operator(&siam, doubleOcc);
    if (siam.checkHigherEnergyDiscarded()) {
      rawData.saveCurrentData();
      rawData.saveqOperator(&doubleOcc, "/doubleOcc" + std::to_string(in));
    } else {
      minIterations++;
    }
    siam.update_internal_state();
    // rfile.write(
    //     siam.all_eigenvalue, // Save the eigenvalue of the current iteration
    //     "Eigenvalues" + std::to_string(in));
  }
  // Start of the Backward Iteration
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Backward iteration starts now.
  // BOOM:: BOOM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%% Backward Iteration %%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // set temperatureArray
  double              vMax    = 2; // Max Voltage can be 2* bandwidth
  double              vMin    = 1e-10;
  int                 tPoints = 200;
  double              delTime = (tPoints - 1.0) / (std::log(vMax / vMin));
  std::vector<double> temperatureArray(tPoints);
  // Log scale
  for (int i = 0; i < tPoints; i++) {
    temperatureArray[i] = vMin * std::exp(std::fabs(i) * 1. / delTime);
  }
  // Load data
  // Start ofd the fdmSpectrum
  fdmThermodynamics fdmTc(&siam, temperatureArray);
  fdmTc.setSystemOperator(&doubleOcc);
  for (size_t in = nMax; in > minIterations; in--) {
    // We dont need load the nrg data for the
    // last iteration
    double enScale =
        std::pow(LAMBDA, -(1.0 * static_cast<double>(in) - 1.0) / 2.0);
    std::cout << "--Started Backward iteration: " << in
              << " Mtr: " << minIterations << " enScale:" << enScale
              << std::endl;
    siam.nrg_iterations_cnt--;
    rawData.loadCurrentData();
    rawData.loadqOperator(&doubleOcc, "/doubleOcc" + std::to_string(in - 1));
    fdmTc.calcThermodynamics(enScale);
    // fdmTc.calcThermodynamics(std::pow(LAMBDA, -(in - 1.0) / 2.0));
  }
  fdmTc.saveFinalData(&rfile);
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
