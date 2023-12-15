#include "models/spinhalf.hpp"
#include "models/spinhalfKondo.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "static/fdmThermodynamics.hpp"
#include "utils/h5stream.hpp"
#include <iostream>
const double LAMBDA = 2.0; // Dont do this
double       hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                         (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  timer  mtime("Total time : ");
  size_t minIterations = 1;
  // double             GAMMA = 0.0100;
  // double fc = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  // double V = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  //
  spinhalf           impurity(0, 0);
  size_t             nMax{41}; // Number of NRG iteration
  h5stream::h5stream rfile("resultsZero.h5");
  // size_t             nMax{42}; // Number of NRG iteration
  // spinhalfKondo      impurity(Jkondo, 0.5);
  // h5stream::h5stream rfile("results.h5");
  spinhalf bathModel(0, 0); // set parameters
  // std::cout << "Done !" << std::endl;
  nrgcore kondoModel(impurity, bathModel);
  // Start calc
  std::cout << "f_dag_operators: " << impurity.f_dag_operator.size()
            << std::endl;
  kondoModel.set_parameters(1024); // set max number of states to be kept
  // siam.add_bath_site({V, V}, 1.0); // first site. This is consistentent with
  // siam.update_internal_state();
  //
  // Bulla's RMP
  //  h5stream::h5stream rfile("resultsOmega" + std::to_string(Omega / Tk) +
  //  ".h5");
  NrgData rawData(&kondoModel);
  for (size_t in = 0; in < nMax; in++) {
    // Rescaling Factor is different for the Kondo Model
    double rescale = std::sqrt(LAMBDA);
    kondoModel.add_bath_site({hopping(in), hopping(in)}, rescale);
    //
    if (kondoModel.checkHigherEnergyDiscarded()) {
      rawData.saveCurrentData();
    } else {
      minIterations++;
    }
    // Update operators if needed here. then ->
    kondoModel.update_internal_state();
    //    std::cout << "Eigenvalues: " << kondoModel.all_eigenvalue <<
    //    std::endl;
    rfile.write(
        kondoModel.all_eigenvalue,
        "Eigenvalues" +
            std::to_string(in)); // Save the eigenvalue of the current iteration
  }
  // Start of the Backward Iteration
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Backword iteration starts now.
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
  fdmThermodynamics fdmTc(&kondoModel, temperatureArray);
  for (size_t in = nMax; in > minIterations; in--) {
    // We dont need load the nrg data for the
    // last iteration
    double enScale = std::pow(LAMBDA, -(1. * in - 1.0) / 2.0);
    std::cout << "--Started Backward iteration: " << in
              << " Mtr: " << minIterations << " enScale:" << enScale
              << std::endl;
    kondoModel.nrg_iterations_cnt--;
    rawData.loadCurrentData();
    fdmTc.calcThermodynamics(enScale);
    // fdmTc.calcThermodynamics(std::pow(LAMBDA, -(in - 1.0) / 2.0));
  }
  fdmTc.saveFinalData(&rfile);
  // Dont do this
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
