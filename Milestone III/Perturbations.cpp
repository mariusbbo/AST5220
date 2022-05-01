#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));

  // Arrays for storing all solutions 
  Vector2D y_solution(Constants.n_ell_tot_full, Vector(n_x * n_k));
  Vector Psi_array(n_x * n_k);
  Vector Pi_array(n_x * n_k);

  // Loop over all wavenumbers
  for (int ik = 0; ik < n_k; ik++) {
      
    // Progress bar...
    if ((10 * ik) / n_k != (10 * ik + 10) / n_k) {
        std::cout << (100 * ik + 100) / n_k << "% " << std::flush;
        if (ik == n_k - 1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find index of x-value and x-value to integrate to
    double idx_end_tight = get_tight_coupling_time(k);
    double x_end_tight = x_array[idx_end_tight];

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system    
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    Vector x_array_tc = Utils::linspace(x_start, x_end_tight, idx_end_tight);
    
    // Set initial conditions
    Vector y_ic{ y_tight_coupling_ini };

    // Solve the ODE
    ODESolver tight_coupling_ode;
    tight_coupling_ode.solve(dydx_tight_coupling, x_array_tc, y_ic);

    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Get last values from tight coupling
    auto y_tight_coupling_end = tight_coupling_ode.get_final_data();

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
     auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling_end, x_end_tight, k);

    // The full ODE system
    ODESolver full_ode;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    // Add +1 to the length to make the final arrays with solutions n_x long after 
    // removing the first elements in the full system below
    Vector x_array_full = Utils::linspace(x_end_tight, x_end, n_x-idx_end_tight+1);

    // Set initial conditions
    Vector y_ic_full{ y_full_ini };

    // Solve the ODE
    full_ode.solve(dydx_full, x_array_full, y_ic_full);

    // Get solutions
    Vector2D y_tc_solution = tight_coupling_ode.get_data_transpose();
    Vector2D y_full_solution = full_ode.get_data_transpose();

    // Constants
    int n_ell_theta_tc = Constants.n_ell_theta_tc;
    int n_ell_theta = Constants.n_ell_theta;
    int n_tc = Constants.n_scalars + Constants.n_ell_theta_tc;
    int n_tot = Constants.n_ell_tot_full;

    // Remove first values in full system to avoid error in the transition between
    // the regimes when creating the splines
    for (int i = 0; i < n_tot; i++) {
        y_full_solution[i].erase(y_full_solution[i].begin());
    }

    // Vector for storing thetas from l=2 to l=7 in tight coupling regime
    Vector2D Theta_l(n_ell_theta - n_ell_theta_tc, Vector(idx_end_tight));

    // Compute and store values for thetas from l=2 to l=7 in tight coupling regime
    for (int l = n_ell_theta_tc; l < n_ell_theta; l++) {
        for (int j = 0; j < idx_end_tight; j++) {

            int l_ind = l - n_ell_theta_tc; // index of theta in the Theta_l 2D vector
            double x = x_array_tc[j];
            double Theta1_j = y_tc_solution[Constants.ind_start_theta_tc + 1][j];

            // Theta_2
            if (l == n_ell_theta_tc) {
                Theta_l[l_ind][j] = -20.0 * Constants.c * k / (45.0 * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x)) * Theta1_j;
            }
            // Theta_l
            else {
                Theta_l[l_ind][j] = -l / (2.0 * l + 1) * Constants.c * k / (cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x)) * Theta_l[l_ind - 1][j];
            }
        }
    }

    // Add solutions from tight coupling to solutions from the full system
    for (int n_i = 0; n_i < n_tot; n_i++) {
        if (n_i < n_tc) {
            y_full_solution[n_i].insert(y_full_solution[n_i].begin(), y_tc_solution[n_i].begin(), y_tc_solution[n_i].end());
        }
        else {
            y_full_solution[n_i].insert(y_full_solution[n_i].begin(), Theta_l[n_i - n_tc].begin(), Theta_l[n_i - n_tc].end());
        }
        
        // Add values to array storing all solutions
        for (int ix = 0; ix < n_x; ix++) {
            y_solution[n_i][ix + ik * n_x] = y_full_solution[n_i][ix];
        }
    }

    // Compute Psi and add values to array storing Pi
    for (int ix = 0; ix < n_x; ix++) {
        if (ix < idx_end_tight){
            Psi_array[ix + ik * n_x] = -y_full_solution[Constants.ind_Phi][ix] - 12.0 * pow(cosmo->H_of_x(0.0) / 
                                        (Constants.c * k * exp(x_array_tc[ix])), 2) * 
                                        cosmo->get_OmegaR(0.0) * y_full_solution[Constants.ind_start_theta + 2][ix];
        }
        else {
            Psi_array[ix + ik * n_x] = -y_full_solution[Constants.ind_Phi][ix] - 12.0 * pow(cosmo->H_of_x(0.0) / 
                                       (Constants.c * k * exp(x_array_full[ix-idx_end_tight])), 2)
                                       * cosmo->get_OmegaR(0.0) * y_full_solution[Constants.ind_start_theta + 2][ix];
        }
        Pi_array[ix + ik * n_x] = y_full_solution[Constants.ind_start_theta + 2][ix];
    }

  }
  Utils::EndTiming("integrateperturbation");


  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================

  delta_cdm_spline.create(x_array, k_array, y_solution[Constants.ind_deltacdm_tc]);
  delta_b_spline.create(x_array, k_array, y_solution[Constants.ind_deltab_tc]);
  v_cdm_spline.create(x_array, k_array, y_solution[Constants.ind_vcdm_tc]);
  v_b_spline.create(x_array, k_array, y_solution[Constants.ind_vb_tc]);
  Phi_spline.create(x_array, k_array, y_solution[Constants.ind_Phi]);
  Psi_spline.create(x_array, k_array, Psi_array);
  Pi_spline.create(x_array, k_array, Pi_array);
  
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);

  for (int l = 0; l < Constants.n_ell_theta; l++) {
      Theta_spline[l].create(x_array, k_array, y_solution[Constants.ind_start_theta + l]);
  }

}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  
  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // Cosmological parameter
  double Hp = cosmo->Hp_of_x(x);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)  
  double Psi = -2.0 / 3.0;
  Phi = -Psi;
  delta_cdm = -3.0 / 2.0 * Psi;
  delta_b = -3.0 / 2.0 * Psi;
  v_cdm = -Constants.c * k / (2.0 * Hp) * Psi;
  v_b = -Constants.c * k / (2.0 * Hp) * Psi;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -Psi / 2.0;
  Theta[1] = Constants.c * k / (6.0 * Hp) * Psi;
  
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // Constants and cosmological parameters
  double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  v_cdm = v_cdm_tc;
  delta_b = delta_b_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -20.0 * c * k / (45.0 * Hp * dtaudx) * Theta[1];
  
  for (int l = n_ell_theta_tc + 1; l < n_ell_theta; l++) {
      Theta[l] = -l / (2.0 * l + 1) * c * k / (Hp * dtaudx) * Theta[l - 1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double idx_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  double npts = n_x;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // Cosmological parameters
  double Hp = 0.0;
  double dtaudx = 0.0;

  // End tight coupling no later than z=2000
  double z_rec_start = 2000;
  double x_rec_start = log(1 / (1 + z_rec_start));
  
  // Check the conditions for tight coupling from Callin
  for (int i = 0; i < npts; i++) {
      Hp = cosmo->Hp_of_x(x_array[i]);
      dtaudx = rec->dtaudx_of_x(x_array[i]);
      
      if (abs(k / Hp / dtaudx) > 0.1 && abs(dtaudx) < 10 || x_array[i] > x_rec_start) {
          idx_tight_coupling_end = i;
          break;
      };
  }
  return idx_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================

  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++) {
    const double x = x_array[ix];

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];
     
      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================

      const double Hp       = cosmo->Hp_of_x(x);
      const double dHpdx    = cosmo->dHpdx_analytical(x);
      const double ddHpddx  = cosmo->ddHpddx_analytical(x);
      const double tau      = rec->tau_of_x(x);
      const double dtaudx   = rec->dtaudx_of_x(x);
      const double ddtauddx  = rec->ddtauddx_of_x(x);
      const double g        = rec->g_tilde_of_x(x);
      const double dgdx     = rec->dgdx_tilde_of_x(x);
      const double ddgddx   = rec->ddgddx_tilde_of_x(x);

      // Perturbation parameters and derivatives
      double v_b = get_v_b(x, k);
      double dv_bdx = v_b_spline.deriv_x(x, k);
      double Theta0 = get_Theta(x, k, 0);
      double Theta1 = get_Theta(x, k, 1);
      double dTheta1dx = Theta_spline[1].deriv_x(x, k);
      double Theta3 = get_Theta(x, k, 3);
      double dTheta3dx = Theta_spline[3].deriv_x(x, k);
      double Psi = get_Psi(x, k);
      double dPsidx = Psi_spline.deriv_x(x, k);
      double dPhidx = Phi_spline.deriv_x(x, k);
      double Pi = get_Pi(x, k);
      double dPidx = Pi_spline.deriv_x(x, k);
      double ddPiddx = 2.0 * Constants.c * k / (5.0 * Hp) * (-dHpdx / Hp * Theta1 + dTheta1dx) + 3.0 / 10.0 * (ddtauddx * Pi + dtaudx * dPidx)
                     - 3.0 * Constants.c * k / (5.0 * Hp) * (-dHpdx / Hp * Theta3 + dTheta3dx);
      //double ddPiddx = Pi_spline.deriv_xx(x, k);

      
      // Temperatur source      
      ST_array[index] = g * (Theta0 + Psi + Pi / 4.0) + exp(-tau) * (dPsidx - dPhidx)
                      - 1.0 / (Constants.c * k) * (dHpdx * g * v_b + Hp * dgdx * v_b + Hp * g * dv_bdx)
                      + 3.0 / (4.0 * pow(Constants.c * k, 2)) * (pow(dHpdx, 2) * g * Pi + Hp * ddHpddx * g * Pi
                      + 3 * Hp * dHpdx * (dgdx * Pi + g * dPidx)
                      + pow(Hp, 2) * (ddgddx * Pi + 2.0 * dgdx * dPidx + g * ddPiddx));

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}





//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  
  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  
  // Constants and cosmological parameters
  double c = Constants.c;
  double a = exp(x);
  double H0 = cosmo->H_of_x(0);
  double OmegaB = cosmo->get_OmegaB(0);
  double OmegaCDM = cosmo->get_OmegaCDM(0);
  double OmegaR = cosmo->get_OmegaR(0);

  // Hp=aH and first derivative
  double Hp = cosmo->Hp_of_x(x);
  double dHpdx = cosmo->dHpdx_analytical(x);
  
  // First and second derivative of optical depth
  double dtaudx = rec->dtaudx_of_x(x);
  double ddtauddx = rec->ddtauddx_of_x(x);

  // Variables
  double ck_over_Hp = c * k / Hp;
  double R = 4.0 * OmegaR / (3.0 * OmegaB * a);
  double Theta2 = -20.0 * ck_over_Hp / (45.0 * dtaudx) * Theta[1];
  double Psi = -Phi - 12.0 * pow(H0 / (c * k * a), 2) * OmegaR * Theta2;
  
  // Rhs differential equations tight coupling
  dPhidx = Psi - pow(ck_over_Hp, 2) / 3.0 * Phi + 0.5 * pow(H0 / Hp, 2) * 
           (OmegaCDM * delta_cdm / a + OmegaB * delta_b / a + 4.0 * OmegaR * Theta[0] / pow(a, 2));
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  double q = (-((1 - R) * dtaudx + (1 + R) * ddtauddx) * (3.0 * Theta[1] + v_b) - ck_over_Hp * Psi +
             (1 - dHpdx / Hp) * ck_over_Hp * (-Theta[0] + 2.0 * Theta2) - ck_over_Hp * dThetadx[0]) / ((1 + R) * dtaudx + dHpdx / Hp - 1);
  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp * Psi;
  ddelta_bdx = ck_over_Hp * v_b - 3.0 * dPhidx;
  dv_bdx = (-v_b - ck_over_Hp * Psi + R * (q + ck_over_Hp * (-Theta[0] + 2.0 * Theta2) - ck_over_Hp * Psi)) / (1 + R);
  dThetadx[1] = (q - dv_bdx) / 3.0;

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_theta         = Constants.n_ell_theta;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

  // Cosmological parameters
  double a = exp(x);
  double H0 = cosmo->H_of_x(0);
  double Hp = cosmo->Hp_of_x(x);
  double eta = cosmo->eta_of_x(x);
  double OmegaR = cosmo->get_OmegaR(0);
  double OmegaB = cosmo->get_OmegaB(0);
  double OmegaCDM = cosmo->get_OmegaCDM(0);
  double dtaudx = rec->dtaudx_of_x(x);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  
  // Constants and variables
  double c = Constants.c;
  double ck_over_Hp = c * k / Hp;
  double R = 4.0 * OmegaR / (3.0 * OmegaB * a);
  double Psi = -Phi - 12.0 * pow(H0 / (c * k * a), 2) * OmegaR * Theta[2];

  // Rhs differential equations full system
  dPhidx = Psi - pow(ck_over_Hp, 2) / 3.0 * Phi + 0.5 * pow(H0 / Hp, 2) * 
           (OmegaCDM * delta_cdm / a + OmegaB * delta_b / a + 4.0 * OmegaR * Theta[0] / pow(a, 2));
  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp * Psi;
  ddelta_bdx = ck_over_Hp * v_b - 3.0 * dPhidx;
  dv_bdx = -v_b - ck_over_Hp * Psi + dtaudx * R * (3.0 * Theta[1] + v_b);
  
  dThetadx[0] = -ck_over_Hp * Theta[1] - dPhidx;
  dThetadx[1] = ck_over_Hp / 3.0 * Theta[0] - 2.0 * ck_over_Hp / 3.0 * Theta[2] + ck_over_Hp / 3.0 * Psi + dtaudx * (Theta[1] + v_b / 3.0);
  int l_max = n_ell_theta - 1;
  dThetadx[l_max] = ck_over_Hp * Theta[l_max - 1] - c * (l_max + 1) / (Hp * eta) * Theta[l_max] + dtaudx * Theta[l_max];

  for (int l = n_ell_theta_tc; l < n_ell_theta - 1; l++) {
      dThetadx[l] = l * ck_over_Hp / (2.0 * l + 1) * Theta[l - 1] - (l + 1) * ck_over_Hp / (2.0 * l + 1) * Theta[l + 1] + 
                    dtaudx * (Theta[l] - 0.1 * Theta[2] * (l == 2));
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";    
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
