#include "../roqj_pop.h"

using namespace std;
using namespace Eigen;

static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

// Parameters...
int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 5;
double tmin = 0., tmax = M_PI, dt = 0.001;
bool printTraj = true;


MatrixXcd H (double t) {
  return cos(t)*sigma_x - sin(t)*sigma_y;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return 0.*id;
}

MatrixXcd Gamma (double t) {
  return 0.*id;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
  return 0.*excited_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}


int main (int argc, char **argv) {
    srand(time(NULL));

    string str_pm = argv[1];
    // For Q_0
    Vector2cd phi_p = excited_state;
    Vector2cd phi_m = ground_state;


    // Now running the simulations
    qubit_roqj_pop jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true);

    string index = "_z_"+str_pm+".txt";
    if (str_pm == "p")
        jump.set_initial_state_R(phi_p);
    else
        jump.set_initial_state_R(phi_m);
    jump.run();
    jump.get_observable("average"+index);
    jump.get_error_observable("error"+index);

  return 0;
}