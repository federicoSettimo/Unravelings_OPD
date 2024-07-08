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


static Eigen::VectorXcd ground_state {{0.,1.}};
static Eigen::VectorXcd excited_state {{1.,0.}};
static Eigen::VectorXcd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::VectorXcd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};
static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

// Parameters...
int N_ensemble = 10000;
double tmin = 0., tmax = 2.*M_PI, dt = 0.001;

// Dynamics parameters
double g = .5, omegaS = 0., omegaE = 1., xi = 1.7;


// 0
double gamma_0 (double t) {
  if (t == 0) return 0;
  return .5*(-((pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*
       ((2*cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
             (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
             sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE) - 
            2*(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))*
             (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
             sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          (2.*pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*
            sqrt(pow(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2) + 
              pow(sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2))) - 
         (2*pow(g,2)*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*
            sqrt(pow(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2) + 
              pow(sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2)))/
          (pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*omegaE)))/
     sqrt(pow(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2) + 
       pow(sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2))));
}

// x
double gamma_x (double t) {
  if (t == 0) return 0;
  return .5*(-0.5*(2*((pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE) + 
         (1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)))*
       ((-2*pow(g,2)*cos((omegaE*t)/2.)*(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))*
            sin((omegaE*t)/2.))/(pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*omegaE) + 
         cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
          ((pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2. + 
            (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2.) - 
         ((g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
            sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)) - 
         (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)) + 
      2*(-(sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)/
            pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))) - 
         (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))*
       (-((cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
              (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2)))/
            pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))) - 
         (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
          (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2)) + 
         (2*pow(g,2)*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          (pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*omegaE) - 
         ((pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2. + 
            (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2.)*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)))/
    (pow((pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
           pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
         cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE) + 
        (1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
         pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)),2) + 
      pow(-(sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)/
           pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))) - 
        (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
           pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
         sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2)));
}

// y
double gamma_y (double t) {
  if (t == 0) return 0;
  return .5*(-0.5*(2*((-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                 pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE) - 
         sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)/
          pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)))*
       (-((cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
              (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2)))/
            pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))) + 
         cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
          (-0.5*(pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE))) + 
            (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2.) + 
         (2*pow(g,2)*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          (pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*omegaE) - 
         (-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                 pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)) + 
      2*((1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)) + 
         (-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                 pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))*
       ((-2*pow(g,2)*cos((omegaE*t)/2.)*(1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))*
            sin((omegaE*t)/2.))/(pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2))*omegaE) + 
         (-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                 pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
            pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
          cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)*
          (g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2)) - 
         ((g*xi*pow(cos((omegaE*t)/2.),2) - g*xi*pow(sin((omegaE*t)/2.),2))*
            sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
          pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)) + 
         (-0.5*(pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE))) + 
            (pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                   pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.)*
               ((-4*pow(g,2)*pow(cos((omegaE*t)/2.),3)*sin((omegaE*t)/2.))/omegaE + 
                 (4*pow(g,2)*cos((omegaE*t)/2.)*pow(sin((omegaE*t)/2.),3))/omegaE - 
                 4*g*cos((omegaE*t)/2.)*sin((omegaE*t)/2.)*(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE)))/2.)*
          sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)))/
    (pow((-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
           pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
         cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE) - 
        sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE)/
         pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)),2) + 
      pow((1 + cos((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE))/
         pow(exp(1),(2*pow(g,2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2)) + 
        (-pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
                pow(-xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.) + 
           pow(exp(1),((-4*pow(g,2)*pow(cos((omegaE*t)/2.),2)*pow(sin((omegaE*t)/2.),2))/pow(omegaE,2) - 
               pow(xi + (2*g*pow(sin((omegaE*t)/2.),2))/omegaE,2))/2.))*
         sin((2*g*xi*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE),2)));
}

// z
double gamma_z (double t) {
  if (t == 0) return 0;
  return .5*(2*pow(g,2)*cos((omegaE*t)/2.)*sin((omegaE*t)/2.))/omegaE;
}

void unravel (string str_init_state, string str_pm);

Matrix2cd proj (const Vector2cd &psi) {return psi*psi.adjoint();}

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}

Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}

double observable (const Matrix2cd &rho) {return real((rho*sigma_x).trace());}

int main () {
    srand(time(NULL));

    unravel("0","p");
    unravel("0","m");
    unravel("x","p");
    unravel("x","m");
    unravel("y","p");
    unravel("y","m");
    unravel("z","p");
    unravel("z","m");

    ofstream out;
    out.open("params.txt");
    out << g << endl << omegaS << endl << omegaE << endl << tmax << endl << dt;
}

void unravel (string str_init_state, string str_pm) {
    string index = "_"+str_init_state+"_"+str_pm+"_1.txt";
    cout << "Unraveling Phi^" << str_init_state << "_t[Q_" << str_init_state << "^";
    if (str_pm == "p")
        cout << "+";
    else cout << "-";
    cout << "]...\n";

    // Initializing the initial state...
    // For Q_0
    Vector2cd phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;

    // For Q_x:
    if (str_init_state == "x") {
      phi_p = plus_state; phi_m = minus_state;
    }

    // For Q_y:
    else if (str_init_state == "y") {
      phi_p = plus_y; phi_m = minus_y;
    }

    // For Q_z:
    else if (str_init_state == "z") {
      phi_p = excited_state; phi_m = ground_state;
    }

    Vector2cd psi, psibar;
    if (str_pm == "p")
        psi = phi_p;
    else psi = phi_m;
    psibar = (sigma_z*psi).normalized();
    Matrix2cd rho_ex = proj(psi);

    // Out files
    ofstream out_ex, out_avg, out_err, out_gamma;
    out_ex.open("analytic"+index);
    out_avg.open("average"+index);
    out_err.open("error"+index);
    out_gamma.open("gamma_"+str_init_state+"_1.txt");

    int Npsi = N_ensemble, Npsibar = 0;
    for (double t = 0.; t < tmax; t += dt) {
        double gamma_deph;
        if (str_init_state == "0")
            gamma_deph = gamma_0(t);
        else if (str_init_state == "x")
            gamma_deph = gamma_x(t);
        else if (str_init_state == "y")
            gamma_deph = gamma_y(t);
        else gamma_deph = gamma_z(t);
        int Npsi_old = Npsi, Npsibar_old = Npsibar;
        double fpsi = (double)Npsi/((double)N_ensemble);

        // Printing...
        out_ex << observable(rho_ex) << endl;
        out_avg << observable(fpsi*proj(psi) + (1.-fpsi)*proj(psibar)) << endl;
        out_err << 0 << endl;
        out_gamma << gamma_deph << endl;

        // Updating rho_ex
        Matrix2cd L = gamma_deph*(sigma_z*rho_ex*sigma_z - .5*anticomm(sigma_z*sigma_z, rho_ex));
        rho_ex = rho_ex + L*dt;
        rho_ex = rho_ex / rho_ex.trace();

        // Now doing the jumps

        // States in psi
        double p_jump;
        if (gamma_deph >= 0.)
            p_jump = gamma_deph*dt;
        else p_jump = -gamma_deph*dt*(1.-fpsi)/fpsi;
        for (int i = 0; i < Npsi_old; ++i) {
            double z = rand()/((double)RAND_MAX);
            if (z < p_jump) {
                Npsi--; Npsibar++;
            }
        }

        // States in psibar
        if (gamma_deph >= 0.)
            p_jump = gamma_deph*dt;
        else p_jump = -gamma_deph*dt*fpsi/(1.-fpsi);
        for (int i = 0; i < Npsibar_old; ++i) {
            double z = rand()/((double)RAND_MAX);
            if (z < p_jump) {
                Npsibar--; Npsi++;
            }
        }
    }
}