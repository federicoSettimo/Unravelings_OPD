// JC with APO technique
// Continnum limit and omega_S = 0
#include "../roqj_pop.h"

using namespace std;
using namespace Eigen;

static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

// Parameters...
int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10, n = 10;
double tmin = 0., tmax = 10., dt = 0.001, omega_c = 1., g = .05, nx = .5*(double)n, nz = (double)n;
bool printTraj = true;

// Functions of the ME
// Map \Phi^{0,x,y}
double R_p (double t) {
  if (t == 0.)
    return 0.;
  return (nx+1.)*g*(1.-cos(omega_c*t))/t;
}
double R_m (double t) {
  if (t == 0.)
    return 0.;
  return nx*g*(1.-cos(omega_c*t))/t;
}
double I_p (double t) {
  if (t == 0.)
    return 0.;
  return (nx+1.)*g*(sin(omega_c*t)/t - omega_c);
}
double I_m (double t) {
  if (t == 0.)
    return 0.;
  return nx*g*(sin(omega_c*t)/t - omega_c);
}
// Map \Phi^z
/*double R_p (double t) {
  if (t == 0.)
    return 0.;
  return (nz+1.)*g*(1.-cos(omega_c*t))/t;
}
double R_m (double t) {
  if (t == 0.)
    return 0.;
  return nz*g*(1.-cos(omega_c*t))/t;
}
double I_p (double t) {
  if (t == 0.)
    return 0.;
  return (nz+1.)*g*(sin(omega_c*t)/t - omega_c);
}
double I_m (double t) {
  if (t == 0.)
    return 0.;
  return nz*g*(sin(omega_c*t)/t - omega_c);
}*/

double gamma_p (double t) {return R_m(t);}
double gamma_m (double t) {return R_p(t);}


MatrixXcd H (double t) {
  return I_p(t)*sigma_p*sigma_m + I_m(t)*sigma_m*sigma_p;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
    return 0.*ground_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}


int main () {
    srand(time(NULL));
    // Defining the initial operators Q_i (use only one and comment out the others)
    // For Q_0:
    /*Vector2cd phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;
    Matrix2cd X0 = .5*(id - sigma_x - sigma_y - sigma_z);
    double mu_p = .5*(1.+sqrt(3)), mu_m = .5*(1.-sqrt(3));*/

    // For Q_x:
    Vector2cd phi_p = plus_state, phi_m = minus_state;
    Matrix2cd X0 = .5*sigma_x;
    double mu_p = .5, mu_m = -.5;

    // For Q_y:
    /*Vector2cd phi_p = plus_y, phi_m = minus_y;
    Matrix2cd X0 = .5*sigma_y;
    double mu_p = .5, mu_m = -.5;*/

    // For Q_z:
    /*Vector2cd phi_p = excited_state, phi_m = ground_state;
    Matrix2cd X0 = .5*sigma_z;
    double mu_p = .5, mu_m = -.5;*/

    // Now running the simulations
    qubit_roqj_pop jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true);

    //string index = "_p.txt";
    string index = "_m.txt";
    if (index == "_p.txt")
        jump.set_initial_state_R(phi_p);
    else
        jump.set_initial_state_R(phi_m);
    jump.run();
    jump.get_observable("average"+index);
    jump.get_error_observable("error"+index);
    jump.get_trajectories ("trajectories"+index);

    ofstream out_g;
    out_g.open("gammas.txt");
    for (double t = 0.; t < tmax; t += dt)
        out_g << gamma_m(t) << " " << gamma_p(t) << endl;


  return 0;
}