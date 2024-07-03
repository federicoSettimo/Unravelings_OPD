#include "roqj_pop.h"
// ------------------------------------- qubit_roqj_pop -------------------------------------
qubit_roqj_pop::qubit_roqj_pop (int N_states, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print, verbose, threshold);
}

void qubit_roqj_pop::run () {
  if (_verbose)
    print_info();

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << 2;
  params.close();

  MatrixXd matrix_observables = MatrixXd::Zero(_num_timesteps, _N_copies);
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    VectorXd this_obs = run_single_iterations(i==0);
    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = this_obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = matrix_observables.row(i).mean();
    _sigma_observable[i] = sqrt((matrix_observables.row(i).array() - _observable[i]).square().sum() / (matrix_observables.row(i).size() - 1));
  }
}

VectorXd qubit_roqj_pop::run_single_iterations (bool verbose) const {
  VectorXd observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_states copies of the initial state
  int N_init = _N_states, N_1 = 0, N_2 = 0;

  // Exact solution
  MatrixXcd rho_ex = projector(_initial_state);
  ofstream out_ex, traj;
  if (verbose) {
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
  }
  
  // _eig_1 and _eig_2 must be constant, otherwise the ensemble size is not finite
  Vector2cd initial_state_t = _initial_state;
  
  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Prints and evolves the exact solution
    if (verbose) {
      out_ex << observable(rho_ex) << endl;
      rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
      rho_ex /= rho_ex.trace();
    }

    // Average state
    Matrix2cd rho = (double)N_init/((double)_N_states)*projector(initial_state_t) + (double)N_1/((double)_N_states)*projector(_eig_1) + (double)N_2/((double)_N_states)*projector(_eig_2);
    rho = rho/rho.trace();
    observables[n_observable] = observable(rho);
    n_observable++;

    int N_1_old = N_1, N_2_old = N_2, N_init_old = N_init;

    // Let's compute all the eigenvalues
    MatrixXcd R = J(projector(initial_state_t),t) + 0.5*(Phi(initial_state_t, t, false)*(initial_state_t.adjoint()) + initial_state_t*(Phi(initial_state_t, t, false).adjoint()));
    ComplexEigenSolver<MatrixXcd> eigs;
    eigs.compute(R);
    Vector2cd eigval = eigs.eigenvalues();
    double lambda_1_init = real(eigval[0])*_dt, lambda_2_init = real(eigval[1])*_dt;
    MatrixXcd eigvec = eigs.eigenvectors();
    // Lets see if the first eigenvalue is actually _eig_1 or it is eig_2. If eig_2, swap
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_init;
      lambda_1_init = lambda_2_init;
      lambda_2_init = tmp;
    }

    R = J(projector(_eig_1),t) + 0.5*(Phi(_eig_1, t, true)*(_eig_1.adjoint()) + _eig_1*(Phi(_eig_1, t, true).adjoint()));
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    double lambda_1_eig_1 = real(eigval[0])*_dt, lambda_2_eig_1 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_eig_1;
      lambda_1_eig_1 = lambda_2_eig_1;
      lambda_2_eig_1 = tmp;
    }

    R = J(projector(_eig_2),t) + 0.5*(Phi(_eig_2, t, true)*(_eig_2.adjoint()) + _eig_2*(Phi(_eig_2, t, true).adjoint()));
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    double lambda_1_eig_2 = real(eigval[0])*_dt, lambda_2_eig_2 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_eig_2;
      lambda_1_eig_2 = lambda_2_eig_2;
      lambda_2_eig_2 = tmp;
    }

    // Cycle in the members in initial_state
    for (int i = 0; i < N_init_old; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (lambda_1_init > 0 && z <= lambda_1_init) {
        N_init--; N_1++;
      }
      else if (lambda_2_init > 0 && z <= lambda_1_init + lambda_2_init) {
        N_init--; N_2++;
      }
    }

    // Cycle in the members in eig_1
    for (int i = 0; i < N_1_old; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (lambda_2_eig_1 > 0 && z <= lambda_2_eig_1) {
        N_1--; N_2++;
      }
      else if (lambda_1_init < 0 && z <= -lambda_1_init*(double)N_init/((double)N_1) + abs(lambda_2_eig_1)) {
        N_init++; N_1--;
        if (_verbose)
          cout << "\tReverse jump from eig_1 to initial_state at time " << t << endl;
      }
      else if (lambda_1_eig_2 < 0 && z >= 1. + lambda_1_eig_2*(double)N_2/((double)N_1)) {
        N_2++; N_1--;
        if (_verbose)
          cout << "\tReverse jump from eig_1 to eig_2 at time " << t << endl;
      }
    }

    // Cycle in the members in eig_2
    for (int i = 0; i < N_2_old; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (lambda_1_eig_2 > 0 && z <= lambda_1_eig_2) {
        N_2--; N_1++;
      }
      else if (lambda_2_init < 0 && z <= -lambda_2_init*(double)N_init/((double)N_2) + abs(lambda_1_eig_2)) {
        N_init++; N_2--;
        if (_verbose)
          cout << "\tReverse jump from eig_2 to initial_state at time " << t << endl;
      }
      else if (lambda_2_eig_1 < 0 && z >= 1. + lambda_2_eig_1*(double)N_2/((double)N_1)) {
        N_1++; N_2--;
        if (_verbose)
            cout << "\tReverse jump from eig_2 to eig_1 at time " << t << endl;
      }
    }

    MatrixXcd K = H(t)  - .5*complex<double>(0.,1.)*Gamma(t);
    initial_state_t = initial_state_t - I*_dt*K*initial_state_t - .5*_dt*Phi(initial_state_t,t,false);
    initial_state_t *= exp(-I*arg(initial_state_t(0)));
    initial_state_t.normalize();
  }
  if (_verbose)
    cout << endl;
  return observables;
}


// Prints the trajectories
void qubit_roqj_pop::get_trajectories (string file_out) {
  if (!_print_trajectory || _N_traj_print <= 0 || file_out == "")
    return;
  if (_verbose)
    cout << "Printing " << _N_traj_print << " trajectories\n";
  
  ofstream out;
  out.open(file_out);

  std::vector<VectorXcd> psi;
  for (int i = 0; i < _N_traj_print; ++i)
    psi.push_back(_initial_state);

  int N_init = _N_traj_print, N_1 = 0, N_2 = 0;

  Vector2cd initial_state_t = _initial_state, initial_state_t_dt;

  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Let's compute all the eigenvalues
    MatrixXcd R = J(projector(initial_state_t),t) + 0.5*(Phi(initial_state_t, t, false)*(initial_state_t.adjoint()) + initial_state_t*(Phi(initial_state_t, t, false).adjoint()));
    ComplexEigenSolver<MatrixXcd> eigs;
    eigs.compute(R);
    Vector2cd eigval = eigs.eigenvalues();
    double lambda_1_init = real(eigval[0])*_dt, lambda_2_init = real(eigval[1])*_dt;
    MatrixXcd eigvec = eigs.eigenvectors();
    // Lets see if the first eigenvalue is actually _eig_1 or it is eig_2. If eig_2, swap
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_init;
      lambda_1_init = lambda_2_init;
      lambda_2_init = tmp;
    }

    R = J(projector(_eig_1),t) + 0.5*(Phi(_eig_1, t, true)*(_eig_1.adjoint()) + _eig_1*(Phi(_eig_1, t, true).adjoint()));
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    double lambda_1_eig_1 = real(eigval[0])*_dt, lambda_2_eig_1 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_eig_1;
      lambda_1_eig_1 = lambda_2_eig_1;
      lambda_2_eig_1 = tmp;
    }

    R = J(projector(_eig_2),t) + 0.5*(Phi(_eig_2, t, true)*(_eig_2.adjoint()) + _eig_2*(Phi(_eig_2, t, true).adjoint()));
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    double lambda_1_eig_2 = real(eigval[0])*_dt, lambda_2_eig_2 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1_eig_2;
      lambda_1_eig_2 = lambda_2_eig_2;
      lambda_2_eig_2 = tmp;
    }

    MatrixXcd K = H(t)  - .5*complex<double>(0.,1.)*Gamma(t);
    initial_state_t_dt = initial_state_t - I*_dt*K*initial_state_t - .5*_dt*Phi(initial_state_t,t,false);
    initial_state_t_dt *= exp(-I*arg(initial_state_t_dt(0)));
    initial_state_t_dt.normalize();

    for (int i = 0; i < _N_traj_print; ++i) {
      out << observable(projector(psi[i])) << " ";
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (psi[i] == initial_state_t) {
        if (lambda_1_init > 0 && z <= lambda_1_init) {
          N_init--; N_1++; psi[i] = _eig_1;
        }
        else if (lambda_2_init > 0 && z <= lambda_1_init + lambda_2_init) {
          N_init--; N_2++; psi[i] = _eig_2;
        }
        else {
          psi[i] = initial_state_t_dt;
        }
      }

      else if (psi[i] == _eig_1) {
        if (lambda_2_eig_1 > 0 && z <= lambda_2_eig_1) {
          N_1--; N_2++; psi[i] = _eig_2;
        }
        else if (lambda_1_init < 0 && z <= -lambda_1_init*(double)N_init/((double)N_1) + lambda_2_eig_1) {
          N_init++; N_1--; psi[i] = initial_state_t_dt;
        }
        else if (lambda_1_eig_2 < 0 && z >= 1. + lambda_1_eig_2*(double)N_2/((double)N_1)) {
          N_2++; N_1--; psi[i] = _eig_2;
        }
      }

      else if (psi[i] == _eig_2) {
        if (lambda_1_eig_2 > 0 && z <= lambda_1_eig_2) {
          N_2--; N_1++; psi[i] = _eig_1;
        }
        else if (lambda_2_init < 0 && z <= -lambda_2_init*(double)N_init/((double)N_2) + lambda_1_eig_2) {
          N_init++; N_2--; psi[i] = initial_state_t_dt;
        }
        else if (lambda_2_eig_1 < 0 && z >= 1. + lambda_2_eig_1*(double)N_2/((double)N_1)) {
          N_1++; N_2--; psi[i] = _eig_1;
        }
      }

    }
    out << endl;

    initial_state_t = initial_state_t_dt;
  }
}











// ------------------------------------- qubit_roqj_pop_mixed -------------------------------------
qubit_roqj_pop_mixed::qubit_roqj_pop_mixed (int N_states, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print, verbose, threshold);
}

// Setting the ensemble
void qubit_roqj_pop_mixed::set_ensemble (const vector<pair<double, Vector2cd>> &ensemble) {
  _ensemble = ensemble;
  double sum_p = 0.;
  for (auto & elem : _ensemble) {
    if (elem.first < 0.)
      elem.first = -elem.first;
    sum_p += elem.first;
  }
  if (sum_p == 0.) {
    set_ensemble();
    return;
  }
  if (sum_p != 1.) {
    for (auto & elem : _ensemble)
      elem.first /= sum_p;
  }
}

void qubit_roqj_pop_mixed::set_ensemble (const vector<double> &probabilities, const vector<Vector2cd> &states) {
  if (probabilities.size() != states.size()) {
    cerr << "Different number of probabilities and states -- EXIT\n";
    exit(EXIT_FAILURE);
  }
  _ensemble.clear();
  vector<pair<double, Vector2cd>> ens;
  for (int i = 0; i < probabilities.size(); ++i)
    ens.push_back(pair<double, Vector2cd>(probabilities[i], states[i]));
  set_ensemble(ens);
}

void qubit_roqj_pop_mixed::set_ensemble () {
  Vector2cd init_state = Vector2cd::Ones(_dim_Hilbert_space).normalized();
  _ensemble.clear();
  _ensemble.push_back(pair<double, Vector2cd>(1., init_state));
}

void qubit_roqj_pop_mixed::add_ensemble (const pair<double, Vector2cd> &state) {
  pair<double, Vector2cd> ens(state);
  if (ens.first == 0.) return;
  if (ens.first < 0.) ens.first = -ens.first;
  _ensemble.push_back(ens);

  double normalization = 0.;
  for (auto & i : _ensemble)
    normalization += i.first;
  if (normalization != 1.) {
    for (auto & i : _ensemble)
      i.first /= normalization;
  }
}

void qubit_roqj_pop_mixed::add_ensemble (double prob, const Vector2cd &state) {
  add_ensemble(pair<double, Vector2cd>(prob, state));
}

// Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
void qubit_roqj_pop_mixed::run () {
  _roqj.initialize(_N_states, _t_i, _t_f, _dt, 1, 2, false, 0, false, _threshold);
  if (_verbose) {
    _roqj.print_info();
    print_ens();
  }

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_copies*_N_traj_print << endl << _dim_Hilbert_space;
  params.close();

  MatrixXd matrix_observables = MatrixXd::Zero(_num_timesteps, _N_copies);
  // Cycle on the copies
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    
    // Cycle on the ensemble members'
    VectorXd obs = VectorXd::Zero(_num_timesteps);
    int n = 1;
    for (auto & elem : _ensemble) {
      _roqj.set_initial_state(elem.second);
      if (_verbose)
        cout << "\tEnsemble member " << n++ << "/" << _ensemble.size() << "...\n";
      obs += elem.first * _roqj.run_single_iterations(false);
    }
    if (_verbose) cout << endl;

    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = matrix_observables.row(i).mean();
    _sigma_observable[i] = sqrt((matrix_observables.row(i).array() - _observable[i]).square().sum() / (matrix_observables.row(i).size() - 1));
  }
}


VectorXd qubit_roqj_pop_mixed::get_exact_sol (string file_out) {
  ofstream out;
  bool verbose = file_out != "";
  if (verbose)
    out.open(file_out);

  MatrixXcd rho = MatrixXcd::Zero(_dim_Hilbert_space, _dim_Hilbert_space);
  for (auto & elem : _ensemble)
    rho += elem.first*projector(elem.second);

  VectorXd obs = VectorXd(_num_timesteps);
  int n = 0;

  for (double t = _t_i; t <= _t_f; t += _dt) {
    double this_obs = observable(rho);
    if (verbose)
      out << this_obs << endl;
    obs[n] = this_obs;
    n++;
    rho = rho + (-I*comm(H(t),rho) + J(rho,t) - 0.5*anticomm(Gamma(t),rho))*_dt;
  }

  return obs;
}

void qubit_roqj_pop_mixed::print_ens () {
  cout << "Ensemble:\n";
  cout << "\tState\t| Probability\n";
  for (auto & elem : _ensemble) {
    cout << "----------------|----------------\n";
    cout << elem.second << "\t| " << elem.first << endl;
  }
  cout << endl;
  return;
}

void qubit_roqj_pop_mixed::get_trajectories (string file_out) {
  if (!_print_trajectory || _N_traj_print <= 0 || file_out == "")
    return;
  if (_verbose)
    cout << "Printing " << _N_traj_print << " trajectories\n";
  
  ofstream out;
  out.open(file_out);

  std::vector<VectorXcd> psi;
  for (int j = 0; j < _ensemble.size(); ++j) {
    for (int i = 0; i < _N_traj_print; ++i)
      psi.push_back(_ensemble[j].second);
  }

  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    for (int j = 0; j < _ensemble.size(); ++j) {
      for (int i = 0; i < _N_traj_print; ++i) {
        int k = i + j*_N_traj_print;
        out << observable(projector(psi[k])) << " ";

        MatrixXcd R = J(projector(psi[k]),t) + 0.5*(Phi(psi[k], t, false)*(psi[k].adjoint()) + psi[k]*(Phi(psi[k], t, false).adjoint()));
        
        // Draws a random number and calculates whether the evolution is deterministic or via a jump
        double z = (double)rand()/((double)RAND_MAX);

        if (z < real(R.trace())*_dt) // Jump
          psi[k] = this->jump(R,z, psi[k], t);
        else {// Free evolution
          MatrixXcd K = H(t)  - .5*complex<double>(0.,1.)*Gamma(t);
          psi[k] = psi[k] - I*_dt*K*psi[k] - .5*_dt*Phi(psi[k],t,false);
        }
        psi[k] *= exp(-I*arg(psi[k](0)));
        psi[k].normalize();
      }
    }
    out << endl;
  }
}


void qubit_roqj_pop::set_initial_state_R (const VectorXcd &psi) {
  set_initial_state(psi);

  ComplexEigenSolver<MatrixXcd> eigs;
  VectorXcd Phi0 = Phi(_initial_state, _t_i, false);
  MatrixXcd R = J(projector(_initial_state),_t_i) + 0.5*(Phi0*(_initial_state.adjoint()) + _initial_state*(Phi0.adjoint()));
  eigs.compute(R);
  MatrixXcd eigvec = eigs.eigenvectors();
  _eig_1 = eigvec.col(0);
  _eig_2 = eigvec.col(1);
  
  // check that both are not nan
  if (isnan(real(_eig_2(0))) || isnan(real(-_eig_2(0)))) {
    if (_eig_1(0) != zero)
      _eig_2 << -conj(_eig_1(1))/conj(_eig_1(0)), 1;
    else _eig_2 << 1, -conj(_eig_1(0))/conj(_eig_1(1));
    _eig_2 = _eig_2.normalized();
  }
  else if (isnan(real(_eig_1(0))) || isnan(real(-_eig_1(0)))) {
    if (_eig_2(0) != zero)
      _eig_1 << -conj(_eig_2(1))/conj(_eig_2(0)), 1;
    else _eig_1 << 1, -conj(_eig_2(0))/conj(_eig_2(1));
    _eig_1 = _eig_1.normalized();
  }
  if (_verbose)
    cout << "Post-jump states:\n" << _eig_1 << endl << _eig_2 << endl;
  return;
}