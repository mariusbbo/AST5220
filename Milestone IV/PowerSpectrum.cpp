    #include"PowerSpectrum.h"

// Constructors
PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

// Function for running all computations
void PowerSpectrum::solve(){

  // Create vectors with values of the wavenumber k
  // The vector with logarithimically spaced values is used when computing Cell
  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector log_k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));

  // Create splines for the spherical Bessel function
  generate_bessel_function_splines();
  
  // Perform line-of-sight integration and create splines for theta_ell
  line_of_sight_integration(k_array);
  
  // Compute the CMB power spectrum and create splines for Cell
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

// Function for computing and creating splines for spherical Bessel function
void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  // Set the number of values used when computing the spherical Bessel function
  double arg_start = k_max * cosmo->eta_of_x(0.0);          
  double dx = 2 * M_PI / 16;                                
  int n_x = arg_start / dx;                                 
  Vector arg_array = Utils::linspace(0, arg_start, n_x);    

  //std::cout << "dx " << dx << "\n";
  //std::cout << "nx " << n_x << "\n";
  //std::cout << "k_max*eta0 " << arg_start << "\n";

  // Compute and store the spherical Bessel function for each ell
  #pragma omp parallel for
  for(int i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    Vector j_ell_array(n_x);
    for (int j = 0; j < n_x; j++) {
        j_ell_array[j] = Utils::j_ell(ell, arg_array[j]);
    }
    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(arg_array, j_ell_array, "j_ell");
  }
  Utils::EndTiming("besselspline");
}


// Function for performing the line-of-sight integration for one source function
Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector& k_array,
    std::function<double(double, double)>& source_function) {
    Utils::StartTiming("lineofsight");

    // Make storage for the results
    Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

    // Set the number of x-values used in the integration
    const int n = 500;
    double dx = 2 * M_PI / n;
    int n_x = -x_start / dx;

    //std::cout << "x_start " << x_start << "\n";
    //std::cout << "dx " << dx << "\n";
    //std::cout << "n_x " << n_x << "\n";

    Vector x_array = Utils::linspace(x_start, 0, n_x);
    double eta0 = cosmo->eta_of_x(0.0);

    // For each ell, perform the line-of-sight integration for all k
    #pragma omp parallel for
    for (int il = 0; il < ells.size(); il++) {

        #pragma omp parallel for
        for (int ik = 0; ik < n_k; ik++) {
            double k = k_array[ik];
            double theta_l = 0.0;

            // Perform line-of-sight integration using the trapezoidal rule
            for (int ix = 0; ix < n_x-1; ix++){
                double arg_i = k * (eta0 - cosmo->eta_of_x(x_array[ix]));                             // Argument for Bessel function time x_i
                double arg_ip = k * (eta0 - cosmo->eta_of_x(x_array[ix+1]));                          // Argument for Bessel function time x_(i+1)
                double theta_i = (source_function(x_array[ix], k) * j_ell_splines[il](arg_i));        // Theta at time x_i
                double theta_ip = (source_function(x_array[ix+1], k) * j_ell_splines[il](arg_ip));    // Theta at time x_(i+1)
                theta_l += (theta_i + theta_ip) * dx;
            }
            // Store the result for Source_ell(k) in results[ell][ik]
            result[il][ik] = theta_l / 2;
        }
    }
    std::cout << "\n";
    Utils::EndTiming("lineofsight");
    return result;
}

// Function for doing the line-of-sight integration for all source functions
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  // Number of ells
  const int nells = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  // Function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Performs line-of-sight integration for each source function
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline result and store it in thetaT_ell_of_k_spline
  for (int i = 0; i < nells; i++) {
      thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i], "theta_ell");
  }

}

// Function for computing the CMB power spectrum
Vector PowerSpectrum::solve_for_cell(
    Vector & k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  
  // Number of ells
  const int nells = ells.size();

  // Vector for storing Cell
  Vector results(nells);

  // Compute and store Cell
  #pragma omp parallel for
  for (int il = 0; il < nells; il++) {
      double Cell = 0.0;
      for (int ik = 0; ik < n_k - 1; ik++) {
          double dk = k_array[ik + 1] - k_array[ik]; // Compute stepsize
          double int_i = primordial_power_spectrum(k_array[ik]) * pow(thetaT_ell_of_k_spline[il](k_array[ik]), 2) / k_array[ik];
          double int_ip = primordial_power_spectrum(k_array[ik+1]) * pow(thetaT_ell_of_k_spline[il](k_array[ik+1]), 2) / k_array[ik+1];
          Cell += (int_i + int_ip) * dk;
      }
      results[il] = 4 * M_PI * Cell / 2;
   }
  return results;
}

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc, n_s - 1.0);
}

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{ 
  double delta_M = 2.0 / 3.0 * pow(Constants.c * k_mpc / cosmo->H_of_x(0), 2) 
                 * pert->get_Psi(x, k_mpc) * exp(x) / (cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0));

  return pow(delta_M, 2) * primordial_power_spectrum(k_mpc);
}

double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

// Output Cells to file
void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(0.0), 2);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    fp << cell_TT_spline( ell ) << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

// Output matter power spectrum to file
void PowerSpectrum::output_matter_pk(std::string filename) const{
    std::ofstream fp(filename.c_str());

    // Find the value of k corresponding to the particle horizon at matter-radiation equality and output this value to file
    double a_eq = (cosmo->get_OmegaR(0.0) + cosmo->get_OmegaNu(0.0)) / (cosmo->get_OmegaB(0.0) + cosmo->get_OmegaCDM(0.0));
    double z_eq = 1 / a_eq - 1;
    double x_eq = log(a_eq);
    double k_eq = cosmo->Hp_of_x(x_eq);
    fp << k_eq * Constants.Mpc / cosmo->get_h() / Constants.c << "\n";

    double eta0 = cosmo->eta_of_x(0.0);
    Vector ell_array = { 6, 100, 200, 500, 1000 };

    // Compute the value of k at which the Bessel function peaks for different ells
    for (int l : ell_array) {
        fp << l / eta0 * Constants.Mpc / cosmo->get_h() << " ";
    }
    fp << "\n";

    const int n_k = 10000;
    auto k_array = Utils::linspace(k_min, k_max, n_k);
    auto print_data = [&](const double k) {
        fp << k * Constants.Mpc / cosmo->get_h() << " ";
        fp << get_matter_power_spectrum(0, k) * 2 * pow(M_PI, 2) / pow(k * Constants.Mpc / cosmo->get_h(), 3) << " ";
        fp << Constants.c * k / cosmo->H_of_x(0.0) << " ";
        fp << thetaT_ell_of_k_spline[4](k) << " ";
        fp << thetaT_ell_of_k_spline[19](k) << " ";
        fp << thetaT_ell_of_k_spline[24](k) << " ";
        fp << thetaT_ell_of_k_spline[32](k) << " ";
        fp << thetaT_ell_of_k_spline[42](k) << " ";
        fp << "\n";
    };
    std::for_each(k_array.begin(), k_array.end(), print_data);
}

// Output spherical Bessel function to file
void PowerSpectrum::output_bessel(std::string filename) const {
    std::ofstream fp(filename.c_str());
    auto arg_array = Utils::linspace(k_max * cosmo->eta_of_x(0.0), 0, 20000);
    auto print_data = [&](const double arg) {
        fp << arg << " ";
        fp << j_ell_splines[0](arg) << " ";
        fp << Utils::j_ell(ells[0], arg) << " ";
        fp << j_ell_splines[10](arg) << " ";
        fp << Utils::j_ell(ells[10], arg) << " ";
        fp << j_ell_splines[25](arg) << " ";
        fp << Utils::j_ell(ells[25], arg) << " ";
        fp << j_ell_splines[40](arg) << " ";
        fp << Utils::j_ell(ells[40], arg) << " ";
        fp << j_ell_splines[60](arg) << " ";
        fp << Utils::j_ell(ells[60], arg) << " ";
        fp << "\n";
    };
    std::for_each(arg_array.begin(), arg_array.end(), print_data);
}



