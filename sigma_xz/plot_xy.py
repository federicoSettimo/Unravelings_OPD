import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Function to read data from a file
def read_file(filename):
    with open(filename, 'r') as file:
        data = np.loadtxt(file)
    return data

params = read_file("params.txt")
gammas = read_file("gammas.txt")

# Q_0
exact_0_p = read_file("analytic_0_p.txt")
exact_0_m = read_file("analytic_0_m.txt")
avg_0_p = read_file("average_0_p.txt")
avg_0_m = read_file("average_0_m.txt")
err_0_p = read_file("error_0_p.txt")
err_0_m = read_file("error_0_m.txt")
traj_0_p = read_file("trajectories_0_p.txt")
traj_0_m = read_file("trajectories_0_m.txt")
mu_0_p = (np.sqrt(3)+1)/2
mu_0_m = (np.sqrt(3)-1)/2
exact_0 = mu_0_p*exact_0_p - mu_0_m*exact_0_m
avg_0 = mu_0_p*avg_0_p - mu_0_m*avg_0_m
err_0 = .5*np.sqrt(err_0_p**2 + err_0_m**2)

# Q_x
exact_x_p = read_file("analytic_x_p.txt")
exact_x_m = read_file("analytic_x_m.txt")
avg_x_p = read_file("average_x_p.txt")
avg_x_m = read_file("average_x_m.txt")
err_x_p = read_file("error_x_p.txt")
err_x_m = read_file("error_x_m.txt")
traj_x_p = read_file("trajectories_x_p.txt")
traj_x_m = read_file("trajectories_x_m.txt")
exact_x = .5*(exact_x_p - exact_x_m)
avg_x = .5*(avg_x_p - avg_x_m)
err_x = .5*np.sqrt(err_x_p**2 + err_x_m**2)

# Q_y
exact_y_p = read_file("analytic_y_p.txt")
exact_y_m = read_file("analytic_y_m.txt")
avg_y_p = read_file("average_y_p.txt")
avg_y_m = read_file("average_y_m.txt")
err_y_p = read_file("error_y_p.txt")
err_y_m = read_file("error_y_m.txt")
traj_y_p = read_file("trajectories_y_p.txt")
traj_y_m = read_file("trajectories_y_m.txt")
exact_y = .5*(exact_y_p - exact_y_m)
avg_y = .5*(avg_y_p - avg_y_m)
err_y = .5*np.sqrt(err_y_p**2 + err_y_m**2)

# Q_z
exact_z_p = read_file("analytic_z_p.txt")
exact_z_m = read_file("analytic_z_m.txt")
exact_z = .5*(exact_z_p - exact_z_m)

rho_t = exact_0 + exact_x + exact_y + exact_z
avg_rho_t = avg_0 + avg_x + avg_y + exact_z
err_rho_t = .5*np.sqrt(err_0**2 + err_x**2 + err_y**2)

t = np.linspace(0, params[3], len(exact_0_p))
Ntraj = int(params[6])


fig, ax = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey = True)
fig.tight_layout()

colors = sns.color_palette("husl", 5)
colors[4] = "black"
markers = ['o', 'x', '^', 's', 'd']
markevery = int(len(t)/50)
markersize=5.5


# Plotting Q^+
ax[0].plot(t, exact_0_p, color=colors[0], label=r'$0$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_0_p, err_0_p, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
for i in range(Ntraj):
    ax[0].plot(t, traj_0_p[:,i], color=colors[0], alpha=.1)
ax[0].plot(t, exact_x_p, color=colors[1], label=r'$x$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_x_p, err_x_p, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
for i in range(Ntraj):
    ax[0].plot(t, traj_x_p[:,i], color=colors[1], alpha=.1)
ax[0].plot(t, exact_y_p, color=colors[2], label=r'$y$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_y_p, err_y_p, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
for i in range(Ntraj):
    ax[0].plot(t, traj_y_p[:,i], color=colors[2], alpha=.1)
ax[0].plot(t, exact_z_p, color=colors[3], label=r'$z$', marker=markers[3], markersize=markersize, markevery=markevery)

# Plotting Q
ax[1].plot(t, exact_0, color=colors[0], label=r'$0$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_0, err_0, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_x, color=colors[1], label=r'$x$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_x, err_x, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_y, color=colors[2], label=r'$y$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_y, err_y, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_z, color=colors[3], label=r'$z$', marker=markers[3], markersize=markersize, markevery=markevery)
ax[1].plot(t, rho_t, color=colors[4], label=r'$\rho_S(t)$', marker=markers[4], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_t, err_rho_t, color=colors[4], marker=markers[4], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)


axx = ax[1].inset_axes([.1,.6,.35,.33])
axx.plot(t, gammas[:,0], label = r'$\gamma_-$')
axx.plot(t, gammas[:,1], label = r'$\gamma_+$')
axx.axhline(0, color="black", linewidth=.5)

fontSize = 25
ax[0].legend(loc='lower right', fontsize=fontSize-1)
ax[1].legend(loc='lower right', fontsize=fontSize-1)
axx.legend(loc = "lower right", fontsize=fontSize-2)

ax[0].set_title(r'$Q_\alpha^+(t)$', fontsize=fontSize)
ax[1].set_title(r'$Q_\alpha(t)$ and $\rho_S(t)$', fontsize=fontSize)
axx.set_title("Rates", fontsize=fontSize-1)

ax[0].set_ylabel(r'tr$[X \sigma_z]$', fontsize=fontSize)
for i in range(2):
    ax[i].set_xlabel(r'$t$', fontsize=fontSize)
    ax[i].tick_params(axis='both', which='major', labelsize=fontSize-1)
axx.tick_params(axis='both', which='major', labelsize=fontSize-2)

plt.show()