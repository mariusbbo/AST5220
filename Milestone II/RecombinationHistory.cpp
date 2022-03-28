#include"RecombinationHistory.h"
#include <math.h>

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  // Set up x-array and make arrays to store X_e(x) and n_e(x)
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector log_ne_arr(npts_rec_arrays);

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  // Calculate recombination history
  for (int i = 0; i < npts_rec_arrays; i++){

      // Compute Xe with the Saha equation
      auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

      // Current values of the electron fraction and number density
      const double Xe_current = Xe_ne_data.first;
      const double ne_current = Xe_ne_data.second;

      // Check if we are in the Saha regime
      // If we are, add current values of Xe and ne to arrays
      if (Xe_current > Xe_saha_limit) {
          Xe_arr[i] = Xe_current;
          log_ne_arr[i] = log(ne_current);
      }

      else {
         
          // Create array with redshift values to check redshift when Saha stop
          Vector z_array(npts_rec_arrays);
          double a;
          for (int k = 0; k < npts_rec_arrays; k++) { 
              z_array[k] = exp(-x_array[k]) - 1;
          }
          // Check at which values of x and z the code switches from Saha to Peebles
          // std::cout << "x_last_saha " << x_array[i] << "\n";
          // std::cout << "z_last_saha " << z_array[i] << "\n";

          // Temporary array with x-values used when solving Peebles equation
          // The array starts at the x-value when Saha stopped
          Vector x_array_new = Utils::linspace(x_array[i], x_end, 1000);


          // The Peebles ODE equation
          ODESolver peebles_Xe_ode;
          ODEFunction dXedx = [&](double x, const double* Xe, double* dXedx) {
              return rhs_peebles_ode(x, Xe, dXedx);
          };

          // Set initial conditions
          double Xe_ini = Xe_arr[i-1]; // Last Xe value above 0.99
          Vector Xe_ic{ Xe_ini };

          // Solve the ODE
          peebles_Xe_ode.solve(dXedx, x_array_new, Xe_ic);

          // Get solution
          auto Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);
          
          // Create spline with solution from Peebles
          Xe_spline_peebles.create(x_array_new, Xe_arr_peebles, "Xe_temp");

          // Fetch cosmological parameters
          const double OmegaB = cosmo->get_OmegaB(0);
          const double H0 = cosmo->get_H0();

          // Critical density today and nuber density of baryons
          double rho_c0 = 3 * pow(H0, 2) / (8 * M_PI * Constants.G);
          double nb = 0.0;

          // Add values of Xe and ne from Peebles to the arrays where the Saha values are stored
          int j = i;
          while (j < npts_rec_arrays) {
              Xe_arr[j] = Xe_spline_peebles(x_array[j]);
              nb = OmegaB * rho_c0 / (Constants.m_H * exp(3 * x_array[j]));
              log_ne_arr[j] = log(Xe_spline_peebles(x_array[j]) * nb);
              j += 1;
          }

          // Create splines for Xe and log(ne)
          Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
          log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

          // Stop loop after Peebles equation has been solved and 
          // splines have been created for Xe and log(ne)
          break;
      }      
  }
  
  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  
  // Current value of a
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double H0          = cosmo->H_of_x(0);
  const double TCMB        = cosmo->get_TCMB(0);

  // Critical density, baryon number density and baryon temperature
  double rho_c0 = 3 * pow(H0, 2) / (8 * M_PI * G);
  double nb = OmegaB * rho_c0 / (m_H * pow(a, 3));
  double Tb = TCMB / a;

  // Compute Xe and ne from the Saha equation
  double C = 1 / (nb * pow(hbar, 3)) * pow(m_e * k_b * Tb / (2 * M_PI), 3.0 / 2.0) * exp(-epsilon_0 / (k_b * Tb));  
  double Xe = 2 / (1 + sqrt(1 + 4/C));
  double ne = Xe * nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB(0);
  const double OmegaR      = cosmo->get_OmegaR(0);
  const double H0          = cosmo->H_of_x(0);
  const double TCMB        = cosmo->get_TCMB(0);
  
  // Critical density
  const double rho_c0 = 3 * pow(H0, 2) / (8 * M_PI * G);

  // Baryon temperature
  double T_b = TCMB / a;

  // dXedx
  double H = cosmo->H_of_x(x);
  double phi2 = 0.448 * log(epsilon_0 / (k_b * T_b));
  double alpha = m_e * sqrt(3 * sigma_T / (8 * M_PI));
  double alpha2 = 64 * M_PI * pow(alpha, 2) * sqrt(epsilon_0) * phi2 * c / (sqrt(27 * M_PI * T_b * k_b) * pow(m_e, 2));
  double beta = alpha2 / pow(hbar, 3) * pow(m_e * k_b * T_b / (2 * M_PI), 3.0 / 2.0) * exp(-epsilon_0 / (k_b * T_b));
  double beta2 = alpha2 / pow(hbar, 3) * pow(m_e * k_b * T_b / (2 * M_PI), 3.0 / 2.0) * exp(-epsilon_0 / (4 * k_b * T_b));
  double nH = (1 - Yp) * 3 * pow(H0, 2) * OmegaB / (8 * M_PI * G * m_H * pow(a, 3));
  double n1s = (1 - X_e) * nH;
  double lambda_alpha = H * pow(3 * epsilon_0 / (c * hbar), 3) / (pow(8 * M_PI, 2) * n1s);
  double Cr = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2);
  
  double dXe_dx =  Cr / H * (beta * (1 - X_e) - nH * alpha2 * pow(X_e, 2));
  dXedx[0] = dXe_dx;

  return GSL_SUCCESS;

}


//=====================================================================
// Solve for the optical depth tau and compute the visibility function
//=====================================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-array, integrate backwards in time
  const int npts = 2000;
  Vector x_array_temp = Utils::linspace(x_end, x_start, npts);

  // The ODE dtau/dx
  ODEFunction dtaudx = [&](double x, const double* tau, double* dtaudx) {
      // Set the derivative for photon optical depth
      dtaudx[0] = -Constants.c * Constants.sigma_T * exp(ne_of_x(x)) / cosmo->H_of_x(x);
      return GSL_SUCCESS;
  };

    // Set initial conditions
    double tau_ini = 0; // Optical depth today
    Vector tau_ic{ tau_ini };

    // Solve the ODE
    ODESolver ode;
    ode.solve(dtaudx, x_array_temp, tau_ic);

    // Get the solution
    auto tau_array_temp = ode.get_data_by_component(0);

    // Arrays for storing the values of x and tau in the right order
    Vector x_array(npts);
    Vector tau_array(npts);
    Vector dtaudx_arr(npts);

    // Reverse the order of the arrays for x and tau because the DE for 
    // tau has been solved backwards in time
    // Compute dtaudx
    for (int i = 0; i < npts; i++) {
        tau_array[i] = tau_array_temp[npts - i - 1];
        x_array[i] = x_array_temp[npts - i - 1];
        dtaudx_arr[i] = -Constants.c * exp(ne_of_x(x_array[i])) * Constants.sigma_T / cosmo->H_of_x(x_array[i]);
    }

    // Create splines
    tau_of_x_spline.create(x_array, tau_array, "tau");
    dtaudx_spline.create(x_array, dtaudx_arr, "dtaudx");
  
    // Visibility function
    Vector g_tilde_array(npts);
    Vector g_tilde_array_temp(npts);

    // Compute sum of g_tilde values for normalisation
    double sum = 0.0;
    for (int i = 0; i < npts; i++) {
        g_tilde_array_temp[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
        sum += g_tilde_array_temp[i];
    }

    // Normalise g_tilde
    for (int i = 0; i < npts; i++) {
        g_tilde_array[i] = g_tilde_array_temp[i] / sum;
    }

    // Create spline for g_tilde
    g_tilde_of_x_spline.create(x_array, g_tilde_array, "g");


  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
    return dtaudx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
    return dtaudx_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return log_ne_of_x_spline(x);
}

//================================================================================
// Create splines for Xe computed only with the Saha equation and the temperature
//================================================================================
void RecombinationHistory::create_Xe_of_x_saha_and_T_spline() {
    
    // Set up x-array and create arrays for storing Xe and T
    Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
    Vector Xe_array_saha(npts_rec_arrays);
    Vector T_array(npts_rec_arrays);
    Vector kB_T_array(npts_rec_arrays);

    // Temperature of CMB today
    double TCMB = cosmo->get_TCMB(0);

    // Compute Xe from the Saha equation
    // Compute temperature
    for (int i = 0; i < npts_rec_arrays; i++) {
        Xe_array_saha[i] = electron_fraction_from_saha_equation(x_array[i]).first;
        T_array[i] = TCMB * exp(-x_array[i]);
        kB_T_array[i] = Constants.k_b * TCMB * exp(-x_array[i]);
    }

    // Create splines for Xe computed only with the Saha equation and T
    Xe_of_x_saha_spline.create(x_array, Xe_array_saha, "Xe_saha");
    T_of_x_spline.create(x_array, T_array, "T");
    kB_T_of_x_spline.create(x_array, kB_T_array, "kB_T");
}

double RecombinationHistory::Xe_of_x_saha(double x) const {
    return Xe_of_x_saha_spline(x);
}

double RecombinationHistory::T_of_x(double x) const {
    return T_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//==========================================================
// Print and output x and z for recombination and decoupling
//==========================================================
void RecombinationHistory::output_times(const std::string filename_times) const {
    double x_rec            = Utils::binary_search_for_value(Xe_of_x_spline, 0.5);
    double x_rec_saha       = Utils::binary_search_for_value(Xe_of_x_saha_spline, 0.5);
    double x_decoupling     = Utils::binary_search_for_value(tau_of_x_spline, 1.0);
    double z_rec            = exp(-x_rec) - 1;
    double z_rec_saha       = exp(-x_rec_saha) - 1;
    double z_decoupling     = exp(-x_decoupling) - 1;
    double Xe_freeze_out    = Xe_of_x(0);
    double x_T_equal_eps    = Utils::binary_search_for_value(kB_T_of_x_spline, Constants.epsilon_0);
    double z_T_equal_eps    = exp(-x_T_equal_eps) - 1;
    double x_dtaudx_equal_1 = Utils::binary_search_for_value(dtaudx_spline, -1.0);
    double z_dtaudx_equal_1 = exp(-x_dtaudx_equal_1) - 1;

    //std::cout << "x_decoupling             = " << x_decoupling      << "\n";
    //std::cout << "x_recombination          = " << x_rec             << "\n";
    //std::cout << "x_recombination_saha     = " << x_rec_saha        << "\n";
    //std::cout << "z_decoupling             = " << z_decoupling      << "\n";
    //std::cout << "z_recombination          = " << z_rec             << "\n";
    //std::cout << "z_recombination_saha     = " << z_rec_saha        << "\n";
    //std::cout << "Freeze out abundance: Xe = " << Xe_freeze_out     << "\n";
    //std::cout << "k_BT = epsilon0:       x = " << x_T_equal_eps     << "\n";
    //std::cout << "k_BT = epsilon0:       z = " << z_T_equal_eps     << "\n";

    std::ofstream write(filename_times.c_str());
    write << x_decoupling << " " << x_rec << " " << x_rec_saha << " " << x_T_equal_eps << " " "\n";
    write << z_decoupling << " " << z_rec << " " << z_rec_saha << " " << z_T_equal_eps << " " "\n";
    write << Xe_freeze_out  << "\n";    
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << Xe_of_x_saha(x)      << " ";
    fp << exp(ne_of_x(x))      << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << T_of_x(x)            << " ";

    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

