#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::Matrix2cd sigma_x {{0,1},{1,0}};
static Eigen::Matrix2cd sigma_y {{0,-I},{I,0}};
static Eigen::Matrix2cd sigma_z {{1,0},{0,-1}};
static Eigen::Matrix2cd sigma_p {{0,1},{0,0}};
static Eigen::Matrix2cd sigma_m {{0,0},{1,0}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

static Eigen::Vector2cd ground_state {{0.,1.}};
static Eigen::Vector2cd excited_state {{1.,0.}};
static Eigen::Vector2cd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::Vector2cd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};
static Eigen::Vector2cd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

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


Matrix2cd H (double t) {
  return I_p(t)*sigma_p*sigma_m + I_m(t)*sigma_m*sigma_p;
}

Matrix2cd J (const Matrix2cd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p;
}

Matrix2cd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m;
}

Vector2cd Phi (const Vector2cd &psi, double t, bool jumped) {
    return 0.*ground_state;
}

double observable (const Matrix2cd &rho) {return real((rho*sigma_x).trace());}

Matrix2cd proj (const Vector2cd &psi) {return psi*psi.adjoint();}

// Function to compute the stochastic increment dW (Wiener process)
complex<double> dW(std::mt19937 &gen, std::normal_distribution<double> &dist) {
    return dist(gen) + I * dist(gen);
}

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

    Vector2cd psi0;
    string index = "_p.txt";
    //string index = "_m.txt";
    if (index == "_p.txt")
        psi0 = phi_p;
    else
        psi0 = phi_m;

    ofstream out_avg, out_traj;
    out_avg.open("QSD_average"+index);
    out_traj.open("QSD_trajectories"+index);




    // Initialize random number generator for Wiener process
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, sqrt(dt));

    vector<Vector2cd> psi_i;
    for (int i = 0; i < N_ensemble; ++i)
        psi_i.push_back(psi0);

    for (double t = 0; t < tmax; t+=dt) {
        // Average state
        Matrix2cd rho = Matrix2cd::Zero();
        // Jump operatrs
        Matrix2cd L_p = sqrt(gamma_p(t))*sigma_p, L_m = sqrt(gamma_m(t))*sigma_m;
        for (int i = 0; i < N_ensemble; ++i) {
            rho += proj(psi_i[i])/((double)N_ensemble);
            if (i < Ntraj)
                out_traj << observable(proj(psi_i[i])) << " ";

            Vector2cd dpsi_det = -I*H(t)*psi_i[i]*dt +
                        (L_p*psi_i[i]).conjugate().dot(L_p*psi_i[i])*psi_i[i]*dt +
                        (L_m*psi_i[i]).conjugate().dot(L_m*psi_i[i])*psi_i[i]*dt -
                        L_p.adjoint()*L_p*psi_i[i]*dt -
                        L_m.adjoint()*L_m*psi_i[i]*dt,
                dpsi_stoch = L_p*psi_i[i]*dW(gen, dist) +
                        L_m*psi_i[i]*dW(gen, dist);
            psi_i[i] += dpsi_det + dpsi_stoch;
            psi_i[i].normalize();
        }
        out_avg << observable(rho) << endl;
        out_traj << endl;
    }

    return 0;
}