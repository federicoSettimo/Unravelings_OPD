// JC with APO technique - no use of roqj
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

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::MatrixXcd sigma_x {{0,1},{1,0}};
static Eigen::MatrixXcd sigma_y {{0,-I},{I,0}};
static Eigen::MatrixXcd sigma_z {{1,0},{0,-1}};
static Eigen::MatrixXcd sigma_p {{0,1},{0,0}};
static Eigen::MatrixXcd sigma_m {{0,0},{1,0}};
static Eigen::MatrixXcd id {{1,0},{0,1}};

// ------------------------- SOME STATES -------------------------
static Eigen::VectorXcd ground_state {{0.,1.}};
static Eigen::VectorXcd excited_state {{1.,0.}};
static Eigen::VectorXcd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::VectorXcd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};
static Eigen::Vector2cd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}

Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B + B*A;}

// Parameters...
int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10, n0 = 1, n1 = 0;
double tmin = 0., tmax = 10., dt = 0.001, g = .5, omega_S = 1., omega_E = .1, nx = .5*(double)(n0+n1), nz = (double)n1;
bool printTraj = true;

// Functions of the ME
// Map \Phi^{0,x,y}
/*double R_p (double t) {return g*g*(nx+1.)*sin((omega_S-omega_E)*t)/(omega_S-omega_E);}
double R_m (double t) {return g*g*nx*sin((omega_S-omega_E)*t)/(omega_S-omega_E);}
double I_p (double t) {return g*g*(nx+1.)*(1.-cos((omega_S-omega_E)*t))/(omega_S-omega_E);}
double I_m (double t) {return g*g*nx*(1.-cos((omega_S-omega_E)*t))/(omega_S-omega_E);}*/
// Map \Phi^z
double R_p (double t) {return g*g*(nz+1.)*sin((omega_S-omega_E)*t)/(omega_S-omega_E);}
double R_m (double t) {return g*g*nz*sin((omega_S-omega_E)*t)/(omega_S-omega_E);}
double I_p (double t) {return g*g*(nz+1.)*(1.-cos((omega_S-omega_E)*t))/(omega_S-omega_E);}
double I_m (double t) {return g*g*nz*(1.-cos((omega_S-omega_E)*t))/(omega_S-omega_E);}

double gamma_p (double t) {return R_m(t);}
double gamma_m (double t) {return R_p(t);}
//double gamma_p (double t) {return gamma_m(t);}

// There should be an Hamiltonian, but let's just put it to zero because of phase covariance
Matrix2cd H (double t) {
  return I_p(t)*sigma_p*sigma_m + I_m(t)*sigma_m*sigma_p;
  //return 0.*sigma_z;
}

Matrix2cd J (const Matrix2cd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p;
}

Matrix2cd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m;
}

double observable (const Matrix2cd &rho) {return real((rho*sigma_z).trace());}

Matrix2cd L_t (const Matrix2cd &rho, double t) {
    return -I*comm(H(t),rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

void unravel (const string &alpha, const string &index); // Unr on Q_alpha, pos or neg part

double rand01 () {return (double)rand()/((double)RAND_MAX);}

int main () {
  srand(time(NULL));
  
  unravel("0","+");
  unravel("0","-");
  unravel("x","+");
  unravel("x","-");
  unravel("y","+");
  unravel("y","-");
  unravel("z","+");
  unravel("z","-");

  return 0;
}

void unravel (const string &alpha, const string &index) {
  // Setting psi0
  Vector2cd psi0, phi_p, phi_m;
  if (alpha == "0") {
    phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;
  }
  else if (alpha == "x") {
    phi_p = plus_state; phi_m = minus_state;
  }
  else if (alpha == "y") {
    phi_p = plus_y; phi_m = minus_y;
  }
  else if (alpha == "z") {
    phi_p = excited_state; phi_m = ground_state;
  }
  else {
    cerr << "Wrong format of alpha\n";
    exit(EXIT_FAILURE);
  }
  if (index == "+") {
    psi0 = phi_p;
    cout << "Unraveling the positive part of Q_" << alpha << endl;
  }
  else if (index == "-") {
    psi0 = phi_m;
    cout << "Unraveling the negative part of Q_" << alpha << endl;
  }
  else {
    cerr << "Wrong format of index\n";
    exit(EXIT_FAILURE);
  }
  cout << "\tUsing " << Ncopies << " copies, each with " << N_ensemble << " ensemble members\n";

  // Out files and initial occupations
  string out_index = "_Phi_z_Q_"+alpha+"_"+index+"_ent.txt";
  ofstream out_ex, out_avg, out_rates, out_err_obs, params;
  out_ex.open("exact"+out_index);
  out_avg.open("avg"+out_index);
  out_rates.open("gamma"+out_index);
  out_err_obs.open("err_obs"+out_index);
  params.open("params.txt");
  params << tmax << endl << Ncopies << endl << N_ensemble << endl;

  int num_timesteps = (int)(tmax/dt) + 1;
  MatrixXd obs_z (Ncopies, num_timesteps);

  // Repeating the cycle Ncopies times...
  for (int icopy = 0; icopy < Ncopies; ++icopy) {
    int delay_print = 50;
    if ((icopy + 1) % delay_print == 0)
      cout << "\tRunning copy " << icopy+1 << "/" << Ncopies << "...\n";
    // Storing the time evolution of all observables
    // Memory inefficient, maybe change?
    VectorXd zs (num_timesteps);

    int Npsi = N_ensemble, N0 = 0, N1 = 0, timestep = 0;
    int N_psi_0 = 0, N_psi_1 = 0, N_0_1_dir = 0, N_0_1_inv = 0, N_1_0_dir = 0, N_1_0_inv = 0, N_0_psi_inv = 0, N_1_psi_inv = 0;
    Matrix2cd rho_ex = projector(psi0);
    Vector2cd psi = psi0;
    for (double t = 0.; t < tmax; t += dt) {
      // Printing the realization-independent things
      if (icopy == 0) {
        out_rates << gamma_p(t) << " " << gamma_m(t) << endl;
        out_ex << observable(rho_ex) << endl;
        // Updating the exact solution
        rho_ex += L_t(rho_ex, t)*dt;
      }

      // Calculating and storing the observables obtained from the realizations
      double fpsi = (double)Npsi/((double)N_ensemble), f1 = (double)N1/((double)N_ensemble), f0 = (double)N0/((double)N_ensemble);
      Matrix2cd rho = fpsi*projector(psi) + f0*projector(ground_state) + f1*projector(excited_state);
      zs[timestep] = observable(rho);

      // Now doing jumps etc
      //double l_psi_0 = norm(psi(0))*gamma_m(t)*dt, l_psi_1 = norm(psi(1))*gamma_p(t)*dt;
      double l_psi_0 = norm(psi.dot(excited_state))*gamma_m(t)*dt, l_psi_1 = norm(psi.dot(ground_state))*gamma_p(t)*dt;
      double l_0_1 = gamma_p(t)*dt, l_1_0 = gamma_m(t)*dt;
      Matrix2cd K = H(t) - .5*I*Gamma(t);
      psi = (id - dt*I*K)*psi;
      psi.normalize();
      int Npsi_old = Npsi, N0_old = N0, N1_old = N1;
      /*cout << "t = " << t << ", f_psi = " << fpsi << ", f_0 = " << f0 << ", f_1 = " << f1 << endl;
      cout << "Direct jump probabilities:\n\tpsi->0: " << l_psi_0 << "\n\tpsi->1: " << l_psi_1 << "\n\t0->1: " << l_0_1 << "\n\t1->0: " << l_1_0 << endl;
      cout << "Reverse jump probabilities:\n\t0->psi: " << max(-l_psi_0*fpsi/f0,0.) << "\n\t1->psi: " << max(-l_psi_1*fpsi/f1,0.) << "\n\t0->1: " << max(-l_1_0*f1/f0,0.) << "\n\t1->0: " << max(-l_0_1*f0/f1,0.) << endl << endl << endl;*/


      for (int i = 0; i < Npsi_old; ++i) {
        double z = rand01();
        if (z <= l_psi_0) {
          Npsi--; N0++; N_psi_0++;
        }
        else if (z >= abs(l_psi_0) && z <= abs(l_psi_0) + l_psi_1) {
          Npsi--; N1++; N_psi_1++;
        }
      }

      // In |0>
      for (int i = 0; i < N0_old; ++i) {
        double z = rand01();
        if (z <= l_0_1) {
          N0--; N1++; N_0_1_dir++;
        }
        //if (l_psi_0 < 0. && z >= abs(l_0_1) && z <= abs(l_0_1) + abs(l_psi_0)*(double)Npsi_old/((double)N0_old) && f0 != 0.) {
        else if (l_psi_0 < 0. && z >= abs(l_0_1) && z <= abs(l_0_1) + abs(l_psi_0)*(double)Npsi/((double)N0) && f0 != 0.) {
          N0--; Npsi++; N_0_psi_inv++;
        }
        //if (l_1_0 < 0. && z > 1. - abs(l_1_0)*(double)N1_old/((double)N0_old) && f0 != 0.) {
        else if (l_1_0 < 0. && z > 1. - abs(l_1_0)*(double)N1/((double)N0) && f0 != 0.) {
          N0--; N1++; N_0_1_inv++;
        }
      }

      // In |1>
      for (int i = 0; i < N1_old; ++i) {
        double z = rand01();
        if (z <= l_1_0) {
          N1--; N0++; N_1_0_dir++;
        }
        //if (l_psi_1 < 0. && z >= abs(l_1_0) && z <= abs(l_1_0) + abs(l_psi_1)*(double)Npsi_old/((double)N1_old) && f1 != 0.) {
        else if (l_psi_1 < 0. && z >= abs(l_1_0) && z <= abs(l_1_0) + abs(l_psi_1)*(double)Npsi/((double)N1) && f1 != 0.) {
          N1--; Npsi++; N_1_psi_inv++;
        }
        //if (l_0_1 < 0. && z > 1. - abs(l_0_1)*(double)N0_old/((double)N1_old) && f1 != 0.) {
        else if (l_0_1 < 0. && z > 1. - abs(l_0_1)*(double)N0/((double)N1) && f1 != 0.) {
          N1--; N0++; N_1_0_inv++;
        }
      }

      timestep++;
    }
    obs_z.row(icopy) = zs;
  }

  // Now taking averages and std_dev
  VectorXd zs_avg (num_timesteps);
  VectorXd zs_err (num_timesteps);
  VectorXd f_psi_err (num_timesteps), f_0_err (num_timesteps), f_1_err (num_timesteps);
  for (int i = 0; i < num_timesteps; ++i) {
    zs_avg(i) = obs_z.col(i).mean();

    zs_err(i) = sqrt((obs_z.col(i).array() - zs_avg(i)).square().sum() / (obs_z.col(i).size() - 1));

    // printing...
    out_avg << zs_avg(i) << endl;
    out_err_obs << zs_err(i) << endl;
  }
}