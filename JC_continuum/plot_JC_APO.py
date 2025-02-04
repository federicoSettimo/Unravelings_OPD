import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import sys

filein = open("params.txt")
Ncopies = int(filein.readline())
Nensemble = int(filein.readline())
ti = float(filein.readline())
tf = float (filein.readline())
dt = float(filein.readline())
print_traj = bool(filein.readline())
Ntraj = int(filein.readline())
dimH = int(filein.readline())
Npoints = int((tf-ti)/dt)

t = np.arange(ti,tf,dt)
markevery = int(Npoints/20)

fig, ax = plt.subplots(1,3, figsize=(16,4), sharex=True, sharey=True)
axx = ax[2].inset_axes([.1,.1,.4,.4])
fig.tight_layout()


trajectories_p = np.zeros((Ntraj, Npoints))
exact_p = np.zeros(Npoints)
avg_obs_p = np.zeros(Npoints)
err_obs_p = np.zeros(Npoints)
gp = np.zeros(Npoints)
gm = np.zeros(Npoints)
non_P_div = np.zeros(Npoints)
if print_traj == True:
    filein = open("trajectories_p.txt")
f_exact = open("analytic_p.txt")
f_avg = open("average_p.txt")
f_err = open("error_p.txt")
f_g = open("gammas.txt")
for i in range(Npoints):
    exact_p[i] = f_exact.readline()
    avg_obs_p[i] = f_avg.readline()
    err_obs_p[i] = f_err.readline()
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            trajectories_p[j,i] = x
            j+=1
    g = f_g.readline()
    gm[i] = g.split()[0]
    gp[i] = g.split()[1]
if print_traj == True:
    for i in range(Ntraj):
        ax[0].plot(t, trajectories_p[i,:], alpha=.1, color='red')
ax[0].plot(t,exact_p,color='black')
ax[0].plot(t,avg_obs_p, marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')
axx.plot(t, gm, label=r'$\gamma_-^x$')
axx.plot(t, gp, '--', label=r'$\gamma_+^x$')
axx.axhline(0, color="black", linewidth=.5)


# Same with QSD
QSD_trajectories_p = np.zeros((Ntraj, Npoints))
QSD_avg_obs_p = np.zeros(Npoints)
if print_traj == True:
    filein = open("QSD_trajectories_p.txt")
f_avg = open("QSD_average_p.txt")
for i in range(Npoints):
    QSD_avg_obs_p[i] = f_avg.readline()
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            QSD_trajectories_p[j,i] = x
            j+=1
if print_traj == True:
    for i in range(Ntraj):
        ax[0].plot(t, QSD_trajectories_p[i,:], alpha=.1, color='blue')
ax[0].plot(t,QSD_avg_obs_p, marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0)


trajectories_m = np.zeros((Ntraj, Npoints))
exact_m = np.zeros(Npoints)
avg_obs_m = np.zeros(Npoints)
err_obs_m = np.zeros(Npoints)
exact = np.zeros(Npoints)
avg_obs = np.zeros(Npoints)
err_obs = np.zeros(Npoints)
if print_traj == True:
    filein = open("trajectories_m.txt")
f_exact = open("analytic_m.txt")
f_avg = open("average_m.txt")
f_err = open("error_m.txt")
for i in range(Npoints):
    exact_m[i] = f_exact.readline()
    avg_obs_m[i] = f_avg.readline()
    err_obs_m[i] = f_err.readline()
    exact[i] = .5*(exact_p[i]-exact_m[i])
    avg_obs[i] = .5*(avg_obs_p[i]-avg_obs_m[i])
    err_obs[i] = np.sqrt(err_obs_p[i]**2 + err_obs_m[i]**2)
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            trajectories_m[j,i] = x
            j+=1
if print_traj == True:
    for i in range(Ntraj):
        ax[1].plot(t, trajectories_m[i,:], alpha=.1, color='red')
ax[1].plot(t,exact_m,color='black')
ax[1].plot(t,avg_obs_m, marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')
ax[2].plot(t,exact,color='black')
ax[2].plot(t,avg_obs, marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')

# Same with QSD
QSD_trajectories_m = np.zeros((Ntraj, Npoints))
QSD_avg_obs_m = np.zeros(Npoints)
if print_traj == True:
    filein = open("QSD_trajectories_m.txt")
f_avg = open("QSD_average_m.txt")
for i in range(Npoints):
    QSD_avg_obs_m[i] = f_avg.readline()
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            QSD_trajectories_m[j,i] = x
            j+=1
if print_traj == True:
    for i in range(Ntraj):
        ax[1].plot(t, QSD_trajectories_m[i,:], alpha=.1, color='blue')
ax[1].plot(t,QSD_avg_obs_m, marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0)
ax[2].plot(t,.5*(QSD_avg_obs_p-QSD_avg_obs_m), marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0)

fontSize = 14
for i in range(3):
    ax[i].set_xlabel(r'$t$', fontsize=fontSize)
ax[0].legend(loc="lower right", fontsize=fontSize)
ax[1].legend(loc="upper right", fontsize=fontSize)
ax[2].legend(loc="lower right", fontsize=fontSize)
axx.legend(loc = "upper right", fontsize=fontSize-1)
ax[0].set_title(r"$Q_\alpha^+$", fontsize=fontSize)
ax[1].set_title(r"$Q_\alpha^-$", fontsize=fontSize)
ax[2].set_title(r"$Q_\alpha = \frac{1}{2}(Q_\alpha^+ - Q_\alpha^-)$", fontsize=fontSize)
axx.set_title("Rates", fontsize=fontSize-1)
ax[0].set_ylabel(r"tr $[X \sigma_x]$", fontsize=fontSize)

for i in range(3):
    ax[i].tick_params(axis='both', which='major', labelsize=fontSize-1)
axx.tick_params(axis='both', which='major', labelsize=fontSize-2)

plt.show()