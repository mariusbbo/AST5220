#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB,
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
    H0 = Constants.H0_over_h * h; // Value Hubble constant today
    OmegaR = (2 * std::pow(3.141592, 2) * pow((Constants.k_b * TCMB), 4) * 8 * 3.14159 * Constants.G /
             (30 * pow(Constants.hbar, 3) * pow(Constants.c, 5) * 3 * pow(H0, 2)) ); // Energy density photons
    OmegaNu = Neff * 7 / 8 * pow((4.0 / 11.0), (4.0 / 3.0)) * OmegaR; // Energy density neutrinos
    OmegaLambda = 1 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu); // Energy density dark energy
}

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");

  // Create array with x-values
  Vector x_array = Utils::linspace(x_start, x_end, 1000);

  // The ODE for deta/dx = c/Hp
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Set initial conditions
  double eta_ini = Constants.c / Hp_of_x(x_array[0]);
  Vector eta_ic{ eta_ini };

  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  // Get the solution
  auto eta_array = ode.get_data_by_component(0);

  // Create spline
  eta_of_x_spline.create(x_array, eta_array, "Eta");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
    return H0 * sqrt((OmegaB + OmegaCDM) * exp(-3 * x) + (OmegaR + OmegaNu) * exp(-4 * x) + OmegaK * exp(-2 * x) + OmegaLambda);
}

double BackgroundCosmology::Hp_of_x(double x) const {
    return exp(x) * H_of_x(x);
}

void BackgroundCosmology::create_Hp_spline(){
    // Create array with x-values
    int npts = 1000;
    Vector x_array = Utils::linspace(x_start, x_end, npts);

    // Create array with Hp values
    Vector Hp_array(npts);
    for (int i = 0; i < npts; i++) {
        double x = x_array[i];
        Hp_array[i] = H0 * sqrt((OmegaB + OmegaCDM) * exp(-x) + (OmegaR + OmegaNu) * exp(-2 * x) + OmegaK + OmegaLambda * exp(2 * x));
    }

    // Create spline
    Hp_of_x_spline.create(x_array, Hp_array, "Hp");
}

double BackgroundCosmology::Hp_of_x_new(double x) const{
    return Hp_of_x_spline(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
    return Hp_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  return Hp_of_x_spline.deriv_xx(x);
}

double BackgroundCosmology::dHpdx_analytical(double x) const{
    double u = (OmegaB + OmegaCDM) * exp(-x) + (OmegaR+OmegaNu) * exp(-2 * x) + OmegaK + OmegaLambda * exp(2 * x);
    double v = -(OmegaB + OmegaCDM) * exp(-x) - 2 * (OmegaR+ OmegaNu) * exp(-2 * x) + 2 * OmegaLambda * exp(2 * x);
    return 0.5 * H0 * v / sqrt(u);
}

double BackgroundCosmology::ddHpddx_analytical(double x) const{
    double u = (OmegaB + OmegaCDM) * exp(-x) + (OmegaR+OmegaNu) * exp(-2 * x) + OmegaK + OmegaLambda * exp(2 * x);
    double v = -(OmegaB + OmegaCDM) * exp(-x) - 2 * (OmegaR + OmegaNu) * exp(-2 * x) + 2 * OmegaLambda * exp(2 * x);
    double w = (OmegaB + OmegaCDM) * exp(-x) + 4 * (OmegaR + OmegaNu) * exp(-2 * x) + 4 * OmegaLambda * exp(2 * x);
    return 0.5 * H0 * (w / sqrt(u) - 0.5 * v * v / (sqrt(u) * u));
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  return (OmegaB * pow(get_H0(), 2)) / (exp(3 * x) * pow(H_of_x(x), 2));;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  return (OmegaR * pow(get_H0(), 2)) / (exp(4 * x) * pow(H_of_x(x), 2));;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  return (OmegaNu * pow(get_H0(), 2)) / (exp(4 * x) * pow(H_of_x(x), 2));;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  return (OmegaCDM * pow(get_H0(), 2)) / (exp(3 * x) * pow(H_of_x(x), 2));;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  return (OmegaLambda * pow(get_H0(), 2)) / pow(H_of_x(x), 2);

}

double BackgroundCosmology::get_OmegaK(double x) const{ 
    if (x == 0.0) return OmegaK;
    return (OmegaK * pow(get_H0(), 2)) / (exp(2 * x) * pow(H_of_x(x), 2));
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::detadx(double x) const {
    return eta_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

void BackgroundCosmology::solve_t_x() {
    // Create array with x-values
    Vector x_array = Utils::linspace(x_start, x_end, 1000);

    // The ODE for dt/dx = 1/H
    ODEFunction dtdx = [&](double x, const double* t, double* dtdx) {
        dtdx[0] = 1 / H_of_x(x);
        return GSL_SUCCESS;
    };

    // Set initial conditions
    double t_ini = 1 / (2 * H_of_x(x_array[0]));
    Vector t_ic{ t_ini };

    // Solve the ODE
    ODESolver ode;
    ode.solve(dtdx, x_array, t_ic);

    // Get the solution
    auto t_array = ode.get_data_by_component(0);

    // Create splines
    t_of_x_spline.create(x_array, t_array, "t");
    x_of_t_spline.create(t_array, x_array, "x");
}

double BackgroundCosmology::t_of_x(double x) const{
    return t_of_x_spline(x);
}

double BackgroundCosmology::x_of_t(double x) const {
    return x_of_t_spline(t_of_x(x));
}

double BackgroundCosmology::chi_of_x(double x) const{
    return eta_of_x(0.0) - eta_of_x(x);
}

double BackgroundCosmology::comoving_distance(double x) const{
    double chi = chi_of_x(x);
    double parenthesis_value = sqrt(abs(OmegaK)) * H0 * chi / Constants.c;
    if (OmegaK < 0)
        return chi * sin(parenthesis_value) / parenthesis_value;
    else if (OmegaK > 0)
        return chi * sinh(parenthesis_value) / parenthesis_value;
    return chi;
}

double BackgroundCosmology::angular_distance(double x) const{
    return exp(x) * comoving_distance(x);
}

double BackgroundCosmology::luminosity_distance(double x) const{ 
    return angular_distance(x) / exp(2 * x);
}

// Function for testing the code
void BackgroundCosmology::tests(const std::string filename_times) const {
    
    // Create array with x-values
    Vector x_array = Utils::linspace(x_start, x_end, 1000);

    // Compute time of last scattering
    double z_ls = 1100;
    double a_ls = 1 / (z_ls + 1);
    double x_ls = log(a_ls);
    double t_ls = t_of_x(x_ls);
    std::cout << "Time of last scattering: " << t_ls / (365 * 24 * 3600) << " years" << "\n";

    double OmegaM = OmegaB + OmegaCDM;
    double OmegaRad = OmegaR + OmegaNu;

    // Check if: 
    // - sum of Omegas = 1
    // - Hp/(e^xH) = 1 
    // - detadx*Hp/c = 1
    for (double x : x_array) {
        //std::cout << "Sum omegas: " << get_OmegaB(x) + get_OmegaCDM(x) + get_OmegaR(x) + get_OmegaLambda(x) + get_OmegaK(x) + get_OmegaNu(x) << "\n";
        //std::cout << "Hp/e^xH: " << Hp_of_x(x) / (exp(x) * H_of_x(x)) << "\n";
        //double detadx = eta_of_x_spline.deriv_x(x);
        //std::cout << "detadxHp/c: " << detadx * Hp_of_x(x) / Constants.c << "\n";
    }

    // Compute time variables at matter-radiation equality
    double a_mr_eq = OmegaRad / OmegaM;
    double x_mr_eq = log(a_mr_eq);
    double z_mr_eq = 1 / a_mr_eq - 1;
    double t_mr_eq = t_of_x(x_mr_eq) / (365 * 24 * 60 * 60);
    std::cout << "Matter-radiation equality:" << "\n";
    std::cout << "x = " << x_mr_eq << "\n";
    std::cout << "z = " << z_mr_eq << "\n";
    std::cout << "t = " << t_mr_eq << " yrs." << "\n";
    std::cout << "\n";

    // Compute time variables at matter-dark energy equality
    double a_mde_eq = pow(OmegaM / OmegaLambda, 1.0/3.0);
    double x_mde_eq = log(a_mde_eq);
    double z_mde_eq = 1 / a_mde_eq - 1;
    double t_mde_eq = t_of_x(x_mde_eq) / (1e9 * 365 * 24 * 60 * 60);
    std::cout << "Matter-dark energy equality:" << "\n";
    std::cout << "x = " << x_mde_eq << "\n";
    std::cout << "z = " << z_mde_eq << "\n";
    std::cout << "t = " << t_mde_eq << " Gyrs." << "\n";
    std::cout << "\n";

    // Compute time variables at time when the Universe started accelerating
    double a_acc = pow(OmegaM / (2 * OmegaLambda), 1.0 / 3.0);
    double z_acc = 1 / a_acc - 1;
    double x_acc = log(a_acc);
    double t_acc = t_of_x(x_acc) / (1e9 * 365 * 24 * 60 * 60);
    std::cout << "Acceleration of Universe starts: " << "\n";
    std::cout << "x = " << x_acc << "\n";
    std::cout << "z = " << z_acc << "\n";
    std::cout << "t = " << t_acc << " Gyrs." << "\n";

    // Write time variables for matter-radiation equality, matter-dark energy equality and
    // the time when the acceleration started to file
    std::ofstream write(filename_times.c_str());
    write << x_mr_eq << " " << z_mr_eq << " " << a_mr_eq << " " << t_mr_eq << "\n";
    write << x_mde_eq << " " << z_mde_eq << " " << a_mde_eq << " " << t_mde_eq << "\n";
    write << x_acc << " " << z_acc << " " << a_acc << " " << t_acc << "\n";
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = x_start;
  const double x_max =  x_end;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                          << " ";
    fp << eta_of_x(x)/Constants.Mpc  << " ";
    fp << H_of_x(x)/get_H0()         << " ";
    fp << Hp_of_x(x)                 << " ";
    fp << dHpdx_of_x(x)/Hp_of_x(x)   << " ";
    fp << ddHpddx_of_x(x)/Hp_of_x(x) << " ";
    fp << get_OmegaB(x)              << " ";
    fp << get_OmegaCDM(x)            << " ";
    fp << get_OmegaLambda(x)         << " ";
    fp << get_OmegaR(x)              << " ";
    fp << get_OmegaNu(x)             << " ";
    fp << get_OmegaK(x)              << " ";
    fp << dHpdx_analytical(x)        << " ";
    fp << ddHpddx_analytical(x)      << " ";
    fp << t_of_x(x)                  << " ";
    fp << comoving_distance(x) / Constants.Mpc           << " ";
    fp << angular_distance(x) / Constants.Mpc            << " ";
    fp << luminosity_distance(x) / Constants.Mpc         << " ";
    fp << H_of_x(x)                 << " ";
    fp << eta_of_x(x)               << " ";
    fp << x_of_t(x)                 << " ";
    fp << exp(x_of_t(x))            << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

