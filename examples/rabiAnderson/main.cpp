#include "fdmThermodynamics.hpp"
#include "models/spinhalf.hpp"
#include "models/spinhalfKondo.hpp"
#include "nrgcore/nrgData.hpp"
#include "nrgcore/nrgcore.hpp"
#include "rabiAnderson.hpp"
#include "rabiSingleParticle.hpp"
#include "rabiSpinless.hpp"
#include "utils/h5stream.hpp"
#include <iomanip>
#include <iostream>
const double LAMBDA = 2.70; // Dont do this
double       hopping(int site) {
  return 0.5 * (1.0 + 1.0 / LAMBDA) * (1. - std::pow(LAMBDA, -site - 1)) /
         std::sqrt((1.0 - std::pow(LAMBDA, -2. * site - 1)) *
                         (1.0 - std::pow(LAMBDA, -2. * site - 3)));
}
int main() {
  timer  mtime("Total time : ");
  size_t minIterations = 1;
  size_t nMax          = 62; // Number of NRG iteration
  double GAMMA         = 0.0080;
  double fc            = 0.5 * std::log(LAMBDA) * (1. + LAMBDA) / (LAMBDA - 1.);
  double V             = std::sqrt(2.0 * fc * GAMMA / std::acos(-1.));
  double UColoumb      = 0.10;
  // Trion
  double Tk           = 1.5e-4;
  double Omega        = 0.000 * Tk;
  double epsilonTrion = 0.14893549596507; // 0.148935587;
  // double epsilonTrion = 0.14893549596505; // 0.148935587;
  double uTrion = UColoumb;
  //
  // spinhalfKondo      impurity(Jkondo, 0.5);
  // rabiSingleParticle impurity(Jkondo, Omega);
  // rabiSpinless       impurity(Jkondo, Omega);
  std::map<std::string, double> params;
  params["UColoumbImpurity"] = UColoumb;
  params["epsilonImpurity"]  = -UColoumb * 0.5;
  params["Omega"]            = Omega;
  params["epsilonTrion"]     = epsilonTrion;
  params["UColoumbTrion"]    = uTrion;
  params["gammaZero"]        = V;
  rabiAnderson impurity(params);
  // h5stream::h5stream rfile("results" + std::to_string(Omega) + ".h5");
  h5stream::h5stream rfile("results.h5");
  spinhalf           bathModel(0, 0); // set parameters
  // std::cout << "Done !" << std::endl;
  nrgcore nrgModel(impurity, bathModel);
  // Start calc
  std::cout << "f_dag_operators: " << impurity.f_dag_operator.size()
            << std::endl;
  nrgModel.set_parameters(1024); // set max number of states to be kept 1800
  NrgData rawData(&nrgModel);
  for (size_t in = 0; in < nMax; in++) {
    // Rescaling Factor is different for the Kondo Model
    // double rescale = std::sqrt(LAMBDA);
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    nrgModel.add_bath_site({hopping(in), hopping(in)}, rescale);
    //
    if (nrgModel.checkHigherEnergyDiscarded()) {
      rawData.saveCurrentData();
    } else {
      minIterations++;
    }
    // Update operators if needed here. then ->
    nrgModel.update_internal_state();
    // std::cout << "Eigenvalues: " << nrgModel.all_eigenvalue <<
    // std::endl;
    // Save the eigenvalue of the current iteration
    rfile.write(nrgModel.all_eigenvalue, "Eigenvalues" + std::to_string(in));
  }
  // calculate the ground state
  auto        grEnergy  = nrgModel.relativeGroundStateEnergy;
  long double absEnergy = 0;
  size_t      ic        = 0;
  for (size_t in = minIterations; ic < nMax - 2; in++) {
    // We dont need load the nrg data for the
    // last iteration
    double enScale = std::pow(LAMBDA, -(1.0 * in - 2.0) / 2.0);
    // std::cout << "--Started Backward iteration: " << in
    //           << " Mtr: " << minIterations << " enScale:" << enScale
    //           << std::endl;
    absEnergy += enScale * grEnergy.at(ic);
    ic++;
  }
  std::cout << grEnergy.size() << " : " << ic << std::endl;
  std::cout << "Absolute Energy: " << std::setprecision(16) << absEnergy
            << std::endl;
  // Start of the Backward Iteration
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Backward iteration starts now.
  // BOOM:: BOOM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // std::cout << "%%% Backward Iteration %%%%%%" << std::endl;
  // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // // set temperatureArray
  // double              vMax    = 2; // Max Voltage can be 2* bandwidth
  // double              vMin    = 1e-10;
  // int                 tPoints = 200;
  // double              delTime = (tPoints - 1.0) / (std::log(vMax / vMin));
  // std::vector<double> temperatureArray(tPoints);
  // // Log scale
  // for (int i = 0; i < tPoints; i++) {
  //     temperatureArray[i] = vMin * std::exp(std::fabs(i) * 1. / delTime);
  // }
  // // Load data
  // // Start ofd the fdmSpectrum
  // fdmThermodynamics fdmTc(&nrgModel, temperatureArray);
  // for (size_t in = nMax; in > minIterations; in--) {
  //     // We dont need load the nrg data for the
  //     // last iteration
  //     double enScale = std::pow(LAMBDA, -(1. * in - 1.0) / 2.0);
  //     std::cout << "--Started Backward iteration: " << in
  //               << " Mtr: " << minIterations << " enScale:" << enScale
  //               << std::endl;
  //     nrgModel.nrg_iterations_cnt--;
  //     rawData.loadCurrentData();
  //     fdmTc.calcThermodynamics(enScale);
  //     // fdmTc.calcThermodynamics(std::pow(LAMBDA, -(in - 1.0) / 2.0));
  // }
  // fdmTc.saveFinalData(&rfile);
  // Dont do this
  rfile.close();
  // create operators that connect to the environment
  return 0;
}
