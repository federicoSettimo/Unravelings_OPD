import matplotlib.pyplot as plt
import numpy as np

# Function to read data from a file
def read_file(filename):
    with open(filename, 'r') as file:
        data = np.loadtxt(file)
    return data

params = read_file("params.txt")

# Q_0
exact_0_p = read_file("analytic_0_p_1.txt")
exact_0_m = read_file("analytic_0_m_1.txt")
avg_0_p = read_file("average_0_p_1.txt")
avg_0_m = read_file("average_0_m_1.txt")
err_0_p = read_file("error_0_p_1.txt")
err_0_m = read_file("error_0_m_1.txt")
gamma_0 = read_file("gamma_0_1.txt")
mu_0_p = (np.sqrt(3)+1)/2
mu_0_m = (np.sqrt(3)-1)/2
exact_0 = mu_0_p*exact_0_p - mu_0_m*exact_0_m
avg_0 = mu_0_p*avg_0_p - mu_0_m*avg_0_m
err_0 = .5*np.sqrt(err_0_p**2 + err_0_m**2)

# Q_x
exact_x_p = read_file("analytic_x_p_1.txt")
exact_x_m = read_file("analytic_x_m_1.txt")
avg_x_p = read_file("average_x_p_1.txt")
avg_x_m = read_file("average_x_m_1.txt")
err_x_p = read_file("error_x_p_1.txt")
err_x_m = read_file("error_x_m_1.txt")
gamma_x = read_file("gamma_x_1.txt")
exact_x = .5*(exact_x_p - exact_x_m)
avg_x = .5*(avg_x_p - avg_x_m)
err_x = .5*np.sqrt(err_x_p**2 + err_x_m**2)

# Q_y
exact_y_p = read_file("analytic_y_p_1.txt")
exact_y_m = read_file("analytic_y_m_1.txt")
avg_y_p = read_file("average_y_p_1.txt")
avg_y_m = read_file("average_y_m_1.txt")
err_y_p = read_file("error_y_p_1.txt")
err_y_m = read_file("error_y_m_1.txt")
gamma_y = read_file("gamma_y_1.txt")
exact_y = .5*(exact_y_p - exact_y_m)
avg_y = .5*(avg_y_p - avg_y_m)
err_y = .5*np.sqrt(err_y_p**2 + err_y_m**2)

# Q_z
exact_z_p = read_file("analytic_z_p_1.txt")
exact_z_m = read_file("analytic_z_m_1.txt")
avg_z_p = read_file("average_z_p_1.txt")
avg_z_m = read_file("average_z_m_1.txt")
err_z_p = read_file("error_z_p_1.txt")
err_z_m = read_file("error_z_m_1.txt")
gamma_z = read_file("gamma_z_1.txt")
exact_z = .5*(exact_z_p - exact_z_m)
avg_z = .5*(avg_z_p - avg_z_m)
err_z = .5*np.sqrt(err_z_p**2 + err_z_m**2)

xi = 2
w0 = 1
wx = 1 + np.exp(-.5*xi*xi)
wy = 1
wz = 1
rho_t = w0*exact_0 + wx*exact_x + wy*exact_y + wz*exact_z
avg_rho_t = w0*avg_0 + wx*avg_x + wy*avg_y + wz*exact_z
err_rho_t = .5*np.sqrt(err_0**2 + err_x**2 + err_y**2)

t = np.linspace(0, params[3], len(exact_0_p))

fig, ax = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey = True)
fig.tight_layout()

colors = ["blue", "red", "green", "yellow", "black"]
markers = ['o', 'x', '^', 's', 'd']

markevery = int(len(exact_0_p)/50)
markersize=4
fontSize = 25

# Plotting Q for xi 1
if False:
    ax[0].plot(t, exact_0, color=colors[0], label=r'$0$', marker=markers[0], markersize=markersize, markevery=len(t))
    ax[0].errorbar(t, avg_0, err_0, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_x, color=colors[1], label=r'$x$', marker=markers[1], markersize=markersize, markevery=len(t))
    ax[0].errorbar(t, avg_x, err_x, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_y, color=colors[2], label=r'$y$', marker=markers[2], markersize=markersize, markevery=len(t))
    ax[0].errorbar(t, avg_y, err_y, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, exact_z, color=colors[3], label=r'$z$', marker=markers[3], markersize=markersize, markevery=len(t))
    ax[0].errorbar(t, avg_z, err_z, color=colors[3], marker=markers[3], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
    ax[0].plot(t, rho_t, color=colors[4], label=r'$\rho_S$', marker=markers[4], markersize=markersize, markevery=len(t))
    ax[0].errorbar(t, avg_rho_t, err_rho_t, color=colors[4], marker=markers[4], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)


ax[0].plot(t, exact_0_p, color=colors[0], label=r'$0$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_0_p, err_0_p, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_x_p, color=colors[1], label=r'$x$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_x_p, err_x_p, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_y_p, color=colors[2], label=r'$y$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_y_p, err_y_p, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_z_p, color=colors[3], label=r'$z$', marker=markers[3], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_z_p, err_z_p, color=colors[3], marker=markers[3], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, rho_t, color=colors[4], label=r'$\rho_S$', marker=markers[4], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_rho_t, err_rho_t, color=colors[4], marker=markers[4], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

ax[1].plot(t, exact_0_m, color=colors[0], label=r'$0$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_0_m, err_0_m, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_x_m, color=colors[1], label=r'$x$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_x_m, err_x_p, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_y_m, color=colors[2], label=r'$y$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_y_m, err_y_m, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, exact_z_m, color=colors[3], label=r'$z$', marker=markers[3], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_z_m, err_z_m, color=colors[3], marker=markers[3], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[1].plot(t, rho_t, color=colors[4], label=r'$\rho_S$', marker=markers[4], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_t, err_rho_t, color=colors[4], marker=markers[4], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

axx = ax[0].inset_axes([.1,.6,.35,.33])
axx.plot(t, gamma_0, color=colors[0])
axx.plot(t, gamma_x, color=colors[1])
axx.plot(t, gamma_y, color=colors[2])
axx.plot(t, gamma_z, color=colors[3])
axx.axhline(0, color="black", linewidth=.5)
axx.set_title("Rates", fontsize=fontSize-1)
axx.tick_params(axis='both', which='major', labelsize=fontSize-2)

ax[0].legend(loc='lower right', fontsize=fontSize-1)
ax[1].legend(loc='lower right', fontsize=fontSize-1)

ax[0].set_title(r'$Q_\alpha(t)$ and $\rho_S(t)$, $\xi = 2$', fontsize=fontSize)
ax[1].set_title(r'$Q_\alpha(t)$ and $\rho_S(t)$, $\xi = $', fontsize=fontSize)

ax[0].set_ylabel(r'tr$[X \sigma_x]$', fontsize=fontSize)
for i in range(2):
    ax[i].set_xlabel(r'$t$', fontsize=fontSize)
    ax[i].tick_params(axis='both', which='major', labelsize=fontSize-1)

plt.show()