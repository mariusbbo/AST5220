#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h		     =  0.67;
  double OmegaB      = 0.05;
  double OmegaCDM	 = 0.267;
  double OmegaK      = 0.0;
  double Neff = 0; // 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp			 = 0.0; //0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  //cosmo.create_Hp_spline();
  //cosmo.solve_t_x();
  cosmo.solve();
  cosmo.info();

  // Run tests to check if results are valid
  //cosmo.tests("times.txt");

  // Output background evolution quantities
  //cosmo.output("cosmology.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();
  //rec.create_Xe_of_x_saha_and_T_spline();

  // Output recombination quantities
  //rec.output("recombination.txt");
  //rec.output_times("rec_times.txt");

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  std::cout << "finished" << "\n";
  pert.info();

  // Output perturbation quantities
  pert.output(0.1 / Constants.Mpc, "perturbations_k0.1.txt");
  pert.output(0.01 / Constants.Mpc, "perturbations_k0.01.txt");
  pert.output(0.001 / Constants.Mpc, "perturbations_k0.001.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
