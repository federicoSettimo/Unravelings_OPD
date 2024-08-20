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

Matrix4cd id {{1., 0., 0., 0.},
              {0., 1., 0., 0.},
              {0., 0., 1., 0.},
              {0., 0., 0., 1.}};
Matrix4cd Q1 {{0., 1., 0., 0.},
              {1., 0., 0., 0.},
              {0., 0., 0., 0.},
              {0., 0., 0., 0.}};
Vector4cd ket0 {{0., 0., 0., 1.}}, ket1 {{0., 0., 1., 0.}}, ket2 {{0., 1., 0., 0.}}, ket3 {{1., 0., 0., 0.}};
Vector4cd ketPlus = (ket3 + ket2)/sqrt(2.), ketMinus = (ket3 - ket2)/sqrt(2.);

// Parameters...
int N_ensemble = 1000, d = 4, Ntraj = 10;
double tmin = 0., tmax = 3., dt = 0.001;

Matrix4cd A (double t) {
    Matrix4cd a {{0., -tan(.25*t), real((I*cos(.5*t) - 2.*sin(.5*t))/(4.*cos(.5*t) + 2.*I*sin(.5*t))), -.5*tan(.25*t)},
                {-tan(.25*t), 0., -3.*tan(.25*t)/4., -.5*tan(.25*t)},
                {-real((I*cos(.5*t) + 2.*sin(.5*t))/(4.*cos(.5*t) - 2.*I*sin(.5*t))), -3.*tan(.25*t)/4., 0., -.25*tan(.25*t)},
                {-.5*tan(.25*t), -.5*tan(.25*t), -.25*tan(.25*t), 0.}};
    return a;
}

Matrix4cd B (double t) {
    Matrix4cd b {{0., 0., imag((I*cos(.5*t) - 2.*sin(.5*t))/(4.*cos(.5*t) + 2.*I*sin(.5*t))), 0.},
                {0., 0., .25, 0},
                {-imag((I*cos(.5*t) + 2.*sin(.5*t))/(4.*cos(.5*t) - 2.*I*sin(.5*t))), -.25, 0., -.25},
                {0., 0., .25, 0}};
    return b;
}

Matrix4cd P (int i) {
    if (i < 0 || i > 3) {
        cerr << "Projecting outside the Hilbert space...\n";
        exit(EXIT_FAILURE);
    }
    Matrix4cd p = Matrix4cd::Zero();
    if (i == 0)
        p(3,3) = 1;
    else if (i == 1)
        p(2,2) = 1;
    else if (i == 2)
        p(1,1) = 1;
    else p(0,0) = 1;
    return p;
}

Matrix4cd H (double t) {
    Matrix4cd h = Matrix4cd::Zero(), Bt = B(t);
    for (int l = 0; l < 4; ++l) {
        complex<double> epsilon_l = 0.;
        for (int k = 0; k < 4; k++)
            epsilon_l += Bt(k,l);
        h += epsilon_l * P(l);
    }
    return h;
}

MatrixXcd f {{1./sqrt(2.), 1./sqrt(6.), 1./(2.*sqrt(3.))},
            {-1./sqrt(2.), 1./sqrt(6.), 1./(2.*sqrt(3.))},
            {0., -sqrt(2./3.), 1./(2.*sqrt(3.))},
            {0., 0., -sqrt(3.)/2.}};

Matrix3cd K (double t) {
    return f.transpose() * A(t) * f;
}

Matrix4cd S (int l) {
    if (l < 0 || l > 2) {
        cerr << "Wrong index for S...\n";
        exit(EXIT_FAILURE);
    }
    Matrix4cd s = Matrix4cd::Zero();
    for (int k = 0; k < l+1; k++) {
        s += P(k);
    }
    return (s - (l+1.)*P(l+1))/sqrt((double)(l+1.)*(l+2.));
}

Vector3d Lindblad_rates (double t) {
    Matrix3cd Kt = K(t);
    ComplexEigenSolver<Matrix3cd> eigs;
    eigs.compute(Kt);
    Vector3cd eigval = eigs.eigenvalues();
    return eigval.real();
}

vector<Matrix4cd> Lindblad_operators (double t) {
    Matrix3cd Kt = K(t);
    ComplexEigenSolver<Matrix3cd> eigs;
    eigs.compute(Kt);
    Matrix3cd U = eigs.eigenvectors();
    Matrix4cd L0 = Matrix4cd::Zero(), L1 = Matrix4cd::Zero(), L2 = Matrix4cd::Zero();
    for (int i = 0; i < 3; ++i) {
        L0 += U(i,0)*S(i);
        L1 += U(i,1)*S(i);
        L2 += U(i,2)*S(i);
    }
    vector<Matrix4cd> ops;
    ops.push_back(L0);
    ops.push_back(L1);
    ops.push_back(L2);
    return ops;
}

Matrix4cd comm (const Matrix4cd &A, const Matrix4cd &B) {return A*B - B*A;}
Matrix4cd anticomm (const Matrix4cd &A, const Matrix4cd &B) {return A*B + B*A;}
Matrix4cd proj (const Vector4cd &psi) {return psi * psi.adjoint();}

Matrix4cd Lindbladian (const Matrix4cd &rho, double t) {
    Vector3d gamma = Lindblad_rates(t);
    vector<Matrix4cd> L = Lindblad_operators(t);
    Matrix4cd J = Matrix4cd::Zero(), Gamma = Matrix4cd::Zero();
    for (int i = 0; i < 3; i++) {
        J += gamma(i)*L[i]*rho*L[i].adjoint();
        Gamma += gamma(i)*L[i].adjoint()*L[i];
    }
    return -I*comm(H(t), rho) + J - .5*anticomm(Gamma, rho);
}

double observable (const Matrix4cd &rho) {return real(rho(1,0));}

void get_exact_solution (const Matrix4cd &rho0, string file_index) {
    string out_file = "exact_"+file_index+".txt";
    ofstream out(out_file);

    Matrix4cd rho_t = rho0;
    for (double t = 0.; t < tmax; t += dt) {
        out << observable(rho_t) << endl;
        rho_t += Lindbladian(rho_t, t)*dt;
        //rho_t /= rho_t.trace();
    }

    return;
}

void MCWF (const Vector4cd &psi0, string file_index) {
    string out_file = "MCWF_"+file_index+".txt", out_file_traj = "MCWF_traj_"+file_index+".txt";
    ofstream out(out_file), out_traj(out_file_traj), out_rates("rates.txt");

    vector<Vector4cd> psi;
    for (int i = 0; i < N_ensemble; ++i)
        psi.push_back(psi0);
    
    for (double t = 0.; t < tmax; t += dt) {
        Matrix4cd rho = Matrix4cd::Zero();

        Matrix4cd Gamma = Matrix4cd::Zero();
        Vector3d gamma = Lindblad_rates(t);
        vector<Matrix4cd> L = Lindblad_operators(t);
        for (int i = 0; i < 3; ++i)
            Gamma += gamma(i)*L[i].adjoint()*L[i];
        Matrix4cd Kt = H(t) - .5*I*Gamma;

        for (int i = 0; i < N_ensemble; ++i) {
            rho += proj(psi[i])/((double)N_ensemble);
            if (i < Ntraj)
                out_traj << observable(proj(psi[i])) << " ";

            double z = rand()/((double)RAND_MAX),
                   pj_0 = gamma(0)*dt*(L[0]*psi[i]).squaredNorm(),
                   pj_1 = gamma(1)*dt*(L[1]*psi[i]).squaredNorm(),
                   pj_2 = gamma(2)*dt*(L[2]*psi[i]).squaredNorm();
            if (z < pj_0)
                psi[i] = L[0]*psi[i];
            else if (z < pj_0 + pj_1)
                psi[i] = L[1]*psi[i];
            else if (z < pj_0 + pj_1 + pj_2)
                psi[i] = L[2]*psi[i];
            else
                psi[i] -= .5*I*dt*Kt*psi[i];
            psi[i].normalize();
        }
        out << observable(rho) << endl;
        out_traj << endl;
        out_rates << gamma(0) << " " << gamma(1) << " " << gamma(2) << endl;
    }
}

complex<double> dW(std::mt19937 &gen) {
  double adjustedSigma = sqrt(dt) / std::sqrt(2.0);
  std::normal_distribution<double> distribution(0.0, adjustedSigma);
  double realPart = distribution(gen);
  double imagPart = distribution(gen);

  // Combine them into a complex number
  return std::complex<double>(realPart, imagPart);
}

void QSD (const Vector4cd &psi0, string file_index) {
    string out_file = "QSD_"+file_index+".txt", out_file_traj = "QSD_traj_"+file_index+".txt";
    ofstream out(out_file), out_traj(out_file_traj);

    std::random_device rd;
    std::mt19937 gen(rd());

    vector<Vector4cd> psi;
    for (int i = 0; i < N_ensemble; ++i)
        psi.push_back(psi0);
    
    for (double t = 0.; t < tmax; t += dt) {
        Matrix4cd rho = Matrix4cd::Zero();

        Vector3d gamma = Lindblad_rates(t);
        vector<Matrix4cd> L = Lindblad_operators(t);
        Matrix4cd L0 = sqrt(gamma(0))*L[0], L1 = sqrt(gamma(1))*L[1], L2 = sqrt(gamma(2))*L[2];

        for (int i = 0; i < N_ensemble; ++i) {
            rho += proj(psi[i])/((double)N_ensemble);
            if (i < Ntraj)
                out_traj << observable(proj(psi[i])) << " ";

            complex<double> avg_L_0 = psi[i].dot(L0*psi[i]), avg_L_1 = psi[i].dot(L1*psi[i]), avg_L_2 = psi[i].dot(L2*psi[i]),
                            avg_L_0_adj = psi[i].dot(L0.adjoint()*psi[i]), avg_L_1_adj = psi[i].dot(L1.adjoint()*psi[i]), avg_L_2_adj = psi[i].dot(L2.adjoint()*psi[i]);
            Vector4cd dpsi = -I*H(t)*psi[i]*dt +
                            (avg_L_0_adj*L0 + L0.adjoint()*L0 + avg_L_0*avg_L_0_adj*id)*psi[i]*dt +
                            (avg_L_1_adj*L1 + L1.adjoint()*L1 + avg_L_1*avg_L_1_adj*id)*psi[i]*dt +
                            (avg_L_2_adj*L2 + L2.adjoint()*L2 + avg_L_2*avg_L_2_adj*id)*psi[i]*dt +
                            (L0 - avg_L_0*id)*psi[i]*dW(gen) +
                            (L1 - avg_L_1*id)*psi[i]*dW(gen) +
                            (L2 - avg_L_2*id)*psi[i]*dW(gen);
            psi[i] += dpsi;
            psi[i].normalize();
        }
        out << observable(rho) << endl;
        out_traj << endl;
    }
}


int main () {
    get_exact_solution(proj(ketPlus), "p");
    get_exact_solution(proj(ketMinus), "m");
    get_exact_solution(Q1, "0");

    MCWF(ketPlus, "p");
    MCWF(ketMinus, "m");

    QSD(ketPlus, "p");
    QSD(ketMinus, "m");

    ofstream out("tmax.txt");
    out << tmax;

    return 0;
}