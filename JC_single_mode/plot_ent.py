import matplotlib.pyplot as plt
import numpy as np

# Function to read data from a file
def read_file(filename):
    with open(filename, 'r') as file:
        data = np.loadtxt(file)
    return data

# Reading files
# First: evol with Phi^0 = Phi^x = Phi^y
avg_Phi_x_Q_0_p = read_file("avg_Phi_x_Q_0_+_ent.txt")
err_Phi_x_Q_0_p = read_file("err_obs_Phi_x_Q_0_+_ent.txt")
exact_Phi_x_Q_0_p = read_file("exact_Phi_x_Q_0_+_ent.txt")
avg_Phi_x_Q_x_p = read_file("avg_Phi_x_Q_x_+_ent.txt")
err_Phi_x_Q_x_p = read_file("err_obs_Phi_x_Q_x_+_ent.txt")
exact_Phi_x_Q_x_p = read_file("exact_Phi_x_Q_x_+_ent.txt")
avg_Phi_x_Q_y_p = read_file("avg_Phi_x_Q_y_+_ent.txt")
err_Phi_x_Q_y_p = read_file("err_obs_Phi_x_Q_y_+_ent.txt")
exact_Phi_x_Q_y_p = read_file("exact_Phi_x_Q_y_+_ent.txt")
avg_Phi_x_Q_z_p = read_file("avg_Phi_x_Q_z_+_ent.txt")
err_Phi_x_Q_z_p = read_file("err_obs_Phi_x_Q_z_+_ent.txt")
exact_Phi_x_Q_z_p = read_file("exact_Phi_x_Q_z_+_ent.txt")

avg_Phi_x_Q_0_m = read_file("avg_Phi_x_Q_0_-_ent.txt")
err_Phi_x_Q_0_m = read_file("err_obs_Phi_x_Q_0_-_ent.txt")
exact_Phi_x_Q_0_m = read_file("exact_Phi_x_Q_0_-_ent.txt")
avg_Phi_x_Q_x_m = read_file("avg_Phi_x_Q_x_-_ent.txt")
err_Phi_x_Q_x_m = read_file("err_obs_Phi_x_Q_x_-_ent.txt")
exact_Phi_x_Q_x_m = read_file("exact_Phi_x_Q_x_-_ent.txt")
avg_Phi_x_Q_y_m = read_file("avg_Phi_x_Q_y_-_ent.txt")
err_Phi_x_Q_y_m = read_file("err_obs_Phi_x_Q_y_-_ent.txt")
exact_Phi_x_Q_y_m = read_file("exact_Phi_x_Q_y_-_ent.txt")
avg_Phi_x_Q_z_m = read_file("avg_Phi_x_Q_z_-_ent.txt")
err_Phi_x_Q_z_m = read_file("err_obs_Phi_x_Q_z_-_ent.txt")
exact_Phi_x_Q_z_m = read_file("exact_Phi_x_Q_z_-_ent.txt")

mu_0_p = (np.sqrt(3)+1)/2
mu_0_m = (np.sqrt(3)-1)/2
avg_Phi_x_Q_0 = mu_0_p*avg_Phi_x_Q_0_p - mu_0_m*avg_Phi_x_Q_0_m
err_Phi_x_Q_0 = .5*np.sqrt(err_Phi_x_Q_0_p**2 + err_Phi_x_Q_0_m**2)
exact_Phi_x_Q_0 = mu_0_p*exact_Phi_x_Q_0_p - mu_0_m*exact_Phi_x_Q_0_m
avg_Phi_x_Q_x = .5*(avg_Phi_x_Q_x_p - avg_Phi_x_Q_x_m)
err_Phi_x_Q_x = .5*np.sqrt(err_Phi_x_Q_x_p**2 + err_Phi_x_Q_x_m**2)
exact_Phi_x_Q_x = .5*(exact_Phi_x_Q_x_p - exact_Phi_x_Q_x_m)
avg_Phi_x_Q_y = .5*(avg_Phi_x_Q_y_p - avg_Phi_x_Q_y_m)
err_Phi_x_Q_y = .5*np.sqrt(err_Phi_x_Q_y_p**2 + err_Phi_x_Q_y_m**2)
exact_Phi_x_Q_y = .5*(exact_Phi_x_Q_y_p - exact_Phi_x_Q_y_m)
avg_Phi_x_Q_z = .5*(avg_Phi_x_Q_z_p - avg_Phi_x_Q_z_m)
err_Phi_x_Q_z = .5*np.sqrt(err_Phi_x_Q_z_p**2 + err_Phi_x_Q_z_m**2)
exact_Phi_x_Q_z = .5*(exact_Phi_x_Q_z_p - exact_Phi_x_Q_z_m)

# Then: evol with Phi^z
avg_Phi_z_Q_0_p = read_file("avg_Phi_z_Q_0_+_ent.txt")
err_Phi_z_Q_0_p = read_file("err_obs_Phi_z_Q_0_+_ent.txt")
exact_Phi_z_Q_0_p = read_file("exact_Phi_z_Q_0_+_ent.txt")
avg_Phi_z_Q_x_p = read_file("avg_Phi_z_Q_x_+_ent.txt")
err_Phi_z_Q_x_p = read_file("err_obs_Phi_z_Q_x_+_ent.txt")
exact_Phi_z_Q_x_p = read_file("exact_Phi_z_Q_x_+_ent.txt")
avg_Phi_z_Q_y_p = read_file("avg_Phi_z_Q_y_+_ent.txt")
err_Phi_z_Q_y_p = read_file("err_obs_Phi_z_Q_y_+_ent.txt")
exact_Phi_z_Q_y_p = read_file("exact_Phi_z_Q_y_+_ent.txt")
avg_Phi_z_Q_z_p = read_file("avg_Phi_z_Q_z_+_ent.txt")
err_Phi_z_Q_z_p = read_file("err_obs_Phi_z_Q_z_+_ent.txt")
exact_Phi_z_Q_z_p = read_file("exact_Phi_z_Q_z_+_ent.txt")

avg_Phi_z_Q_0_m = read_file("avg_Phi_z_Q_0_-_ent.txt")
err_Phi_z_Q_0_m = read_file("err_obs_Phi_z_Q_0_-_ent.txt")
exact_Phi_z_Q_0_m = read_file("exact_Phi_z_Q_0_-_ent.txt")
avg_Phi_z_Q_x_m = read_file("avg_Phi_z_Q_x_-_ent.txt")
err_Phi_z_Q_x_m = read_file("err_obs_Phi_z_Q_x_-_ent.txt")
exact_Phi_z_Q_x_m = read_file("exact_Phi_z_Q_x_-_ent.txt")
avg_Phi_z_Q_y_m = read_file("avg_Phi_z_Q_y_-_ent.txt")
err_Phi_z_Q_y_m = read_file("err_obs_Phi_z_Q_y_-_ent.txt")
exact_Phi_z_Q_y_m = read_file("exact_Phi_z_Q_y_-_ent.txt")
avg_Phi_z_Q_z_m = read_file("avg_Phi_z_Q_z_-_ent.txt")
err_Phi_z_Q_z_m = read_file("err_obs_Phi_z_Q_z_-_ent.txt")
exact_Phi_z_Q_z_m = read_file("exact_Phi_z_Q_z_-_ent.txt")

avg_Phi_z_Q_0 = mu_0_p*avg_Phi_z_Q_0_p - mu_0_m*avg_Phi_z_Q_0_m
err_Phi_z_Q_0 = .5*np.sqrt(err_Phi_z_Q_0_p**2 + err_Phi_z_Q_0_m**2)
exact_Phi_z_Q_0 = mu_0_p*exact_Phi_z_Q_0_p - mu_0_m*exact_Phi_z_Q_0_m
avg_Phi_z_Q_x = .5*(avg_Phi_z_Q_x_p - avg_Phi_z_Q_x_m)
err_Phi_z_Q_x = .5*np.sqrt(err_Phi_z_Q_x_p**2 + err_Phi_z_Q_x_m**2)
exact_Phi_z_Q_x = .5*(exact_Phi_z_Q_x_p - exact_Phi_z_Q_x_m)
avg_Phi_z_Q_y = .5*(avg_Phi_z_Q_y_p - avg_Phi_z_Q_y_m)
err_Phi_z_Q_y = .5*np.sqrt(err_Phi_z_Q_y_p**2 + err_Phi_z_Q_y_m**2)
exact_Phi_z_Q_y = .5*(exact_Phi_z_Q_y_p - exact_Phi_z_Q_y_m)
avg_Phi_z_Q_z = .5*(avg_Phi_z_Q_z_p - avg_Phi_z_Q_z_m)
err_Phi_z_Q_z = .5*np.sqrt(err_Phi_z_Q_z_p**2 + err_Phi_z_Q_z_m**2)
exact_Phi_z_Q_z = .5*(exact_Phi_z_Q_z_p - exact_Phi_z_Q_z_m)

markevery = int(len(avg_Phi_x_Q_z_p)/50)
markersize=4
params = open("params.txt")
tmax = float(params.readline())
t = np.linspace(0, tmax, avg_Phi_x_Q_z_p.shape[0])


fig, ax = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey = True)
fig.tight_layout()


ax[0].plot(t, exact_Phi_x_Q_0, color='black', label=r'$\Phi^0_t[Q_0]$', marker='o', markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_0, err_Phi_x_Q_0, color='black', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_x_Q_x, color='red', label=r'$\Phi^x_t[Q_x]$', marker='x', markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_x, err_Phi_x_Q_x, color='red', marker='x', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_x_Q_y, color='green', label=r'$\Phi^y_t[Q_y]$', marker='^', markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_y, err_Phi_x_Q_y, color='green', marker='^', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_z_Q_z, color='blue', label=r'$\Phi^z_t[Q_z]$', marker='s', markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_z_Q_z, err_Phi_z_Q_z, color='blue', marker='s', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

if False:
    ax[0].plot(t, exact_Phi_x_Q_0_p, color='black', label=r'$\Phi^0_t[Q_0^+]$')
    ax[0].errorbar(t, avg_Phi_x_Q_0_p, err_Phi_x_Q_0_p, color='black', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_Phi_x_Q_x_p, color='red', label=r'$\Phi^x_t[Q_x^+]$')
    ax[0].errorbar(t, avg_Phi_x_Q_x_p, err_Phi_x_Q_x_p, color='red', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_Phi_x_Q_y_p, color='green', label=r'$\Phi^y_t[Q_y^+]$')
    ax[0].errorbar(t, avg_Phi_x_Q_y_p, err_Phi_x_Q_y_p, color='green', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_Phi_z_Q_z_p, color='blue', label=r'$\Phi^z_t[Q_z^+]$')
    ax[0].errorbar(t, avg_Phi_z_Q_z_p, err_Phi_z_Q_z_p, color='blue', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

    ax[1].plot(t, exact_Phi_x_Q_0_m, color='black', label=r'$\Phi^0_t[Q_0^-]$')
    ax[1].errorbar(t, avg_Phi_x_Q_0_m, err_Phi_x_Q_0_m, color='black', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[1].plot(t, exact_Phi_x_Q_x_m, color='red', label=r'$\Phi^x_t[Q_x^-]$')
    ax[1].errorbar(t, avg_Phi_x_Q_x_m, err_Phi_x_Q_x_m, color='red', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[1].plot(t, exact_Phi_x_Q_y_m, color='green', label=r'$\Phi^y_t[Q_y^-]$')
    ax[1].errorbar(t, avg_Phi_x_Q_y_m, err_Phi_x_Q_y_m, color='green', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[1].plot(t, exact_Phi_z_Q_z_m, color='blue', label=r'$\Phi^x_t[Q_z^-]$')
    ax[1].errorbar(t, avg_Phi_z_Q_z_m, err_Phi_z_Q_z_m, color='blue', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)



# Now rho_S and some repreparations
w_0 = 1
w_x = 1
w_y = 1
w_z = 1
exact_rho_S = w_0*exact_Phi_x_Q_0 + w_x*exact_Phi_x_Q_x + w_y*exact_Phi_x_Q_y + w_z*exact_Phi_z_Q_z
avg_rho_S = w_0*avg_Phi_x_Q_0 + w_x*avg_Phi_x_Q_x + w_y*avg_Phi_x_Q_y + w_z*avg_Phi_z_Q_z
err_rho_S = .25*np.sqrt(err_Phi_x_Q_0**2 + err_Phi_x_Q_x**2 + err_Phi_x_Q_y**2 + err_Phi_z_Q_z**2)
ax[1].plot(t, exact_rho_S, color='black', label=r'$\mathcal{R}=\operatorname{id}$', marker='o', markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_S, err_rho_S, color='black', marker='o', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R[X] = \sigma_x X \sigma_x
r00 = 1
r0y = 2
r0z = 2
rxx = 1
ryy = -1
rzz = -1
exact_R1 = w_0*(r00*exact_Phi_x_Q_0 + r0y*exact_Phi_x_Q_y + r0z*exact_Phi_x_Q_z) + w_x*rxx*exact_Phi_x_Q_x + w_y*ryy*exact_Phi_x_Q_y + w_z*rzz*exact_Phi_z_Q_z;
avg_R1 = w_0*(r00*avg_Phi_x_Q_0 + r0y*avg_Phi_x_Q_y + r0z*avg_Phi_x_Q_z) + w_x*rxx*avg_Phi_x_Q_x + w_y*ryy*avg_Phi_x_Q_y + w_z*rzz*avg_Phi_z_Q_z;
ax[1].plot(t, exact_R1, color='green', label=r'$\mathcal{R}=\sigma_x\cdot\sigma_x$', marker='x', markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_R1, err_rho_S, color='green', marker='x', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R[X] = <0|X|0> |-><-|
if False:
    rz0 = -.5
    rzy = -.5
    rzz = -.5
    r0z = 1
    r0y = 1
    r0z = 1
    exact_R2 = (w_0*(r00*exact_Phi_x_Q_0 + r0y*exact_Phi_x_Q_y + r0z*exact_Phi_x_Q_z) + w_z*(rz0*exact_Phi_z_Q_0 + rzy*exact_Phi_z_Q_y + rzz*exact_Phi_z_Q_z))/(.5)
    avg_R2 = (w_0*(r00*avg_Phi_x_Q_0 + r0y*avg_Phi_x_Q_y + r0z*avg_Phi_x_Q_z) + w_z*(rz0*avg_Phi_z_Q_0 + rzy*avg_Phi_z_Q_y + rzz*avg_Phi_z_Q_z))/(.5)
    ax[1].plot(t, exact_R2, color='blue', label=r'$\mathcal{R}=\langle0\vert\cdot\vert0\rangle\,\vert-\rangle\langle-\vert$', marker='s', markersize=markersize, markevery=len(t))
    ax[1].errorbar(t, avg_R2, err_rho_S, color='blue', marker='s', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R[X] = <0|X|0> |0><0|
rz0 = -.5
rzx = -.5
rzy = -.5
r00 = 1
r0x = 1
r0y = 1
exact_R2 = (w_0*(r00*exact_Phi_x_Q_0 + r0x*exact_Phi_x_Q_x + r0y*exact_Phi_x_Q_y) + w_z*(rz0*exact_Phi_z_Q_0 + rzx*exact_Phi_z_Q_x + rzy*exact_Phi_z_Q_y))/(.5)
avg_R2 = (w_0*(r00*avg_Phi_x_Q_0 + r0x*avg_Phi_x_Q_x + r0y*avg_Phi_x_Q_y) + w_z*(rz0*avg_Phi_z_Q_0 + rzx*avg_Phi_z_Q_x + rzy*avg_Phi_z_Q_y))/(.5)
ax[1].plot(t, exact_R2, color='blue', label=r'$\mathcal{R}=\langle0\vert\cdot\vert0\rangle\,\vert0\rangle\langle0\vert$', marker='s', markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_R2, err_rho_S, color='blue', marker='s', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R[X] = <1|X|1> |1><1|
rz0 = .5
rzx = .5
rzy = .5
rzz = 1
exact_R3 = w_z*(rz0*exact_Phi_z_Q_0 + rzx*exact_Phi_z_Q_x + rzy*exact_Phi_z_Q_y + rzz*exact_Phi_z_Q_z)/(.5)
avg_R3 = w_z*(rz0*avg_Phi_z_Q_0 + rzx*avg_Phi_z_Q_x + rzy*avg_Phi_z_Q_y + rzz*avg_Phi_z_Q_z)/(.5)
ax[1].plot(t, exact_R3, color='red', label=r'$\mathcal{R}=\langle1\vert\cdot\vert1\rangle\,\vert1\rangle\langle1\vert$', marker='^', markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_R3, err_rho_S, color='red', marker='^', markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)



fontSize = 25
ax[0].legend(loc='lower right', fontsize=fontSize-1)
ax[1].legend(loc='lower right', fontsize=fontSize-1)

ax[0].set_title(r'$\Phi^\alpha_t[Q_\alpha]$', fontsize=fontSize)
ax[1].set_title(r'tr$_E[U(t)(\mathcal{R}\otimes\operatorname{id})\rho_{SE} U^\dagger(t)]$', fontsize=fontSize)


ax[0].set_ylabel(r'tr$[X \sigma_z]$', fontsize=fontSize)
for i in range(2):
    ax[i].set_xlabel(r'$t$', fontsize=fontSize)
    ax[i].tick_params(axis='both', which='major', labelsize=fontSize-1)

plt.show()