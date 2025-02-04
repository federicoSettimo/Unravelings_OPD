import matplotlib.pyplot as plt
import seaborn as sns
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

# Now reading zero discord
avg_Phi_0_Q_0_p_0d = read_file("avg_Phi_0_Q_0_+_zero_discord.txt")
err_Phi_0_Q_0_p_0d = read_file("err_obs_Phi_0_Q_0_+_zero_discord.txt")
exact_Phi_0_Q_0_p_0d = read_file("exact_Phi_0_Q_0_+_zero_discord.txt")
avg_Phi_0_Q_x_p_0d = read_file("avg_Phi_0_Q_x_+_zero_discord.txt")
err_Phi_0_Q_x_p_0d = read_file("err_obs_Phi_0_Q_x_+_zero_discord.txt")
exact_Phi_0_Q_x_p_0d = read_file("exact_Phi_0_Q_x_+_zero_discord.txt")
avg_Phi_0_Q_y_p_0d = read_file("avg_Phi_0_Q_y_+_zero_discord.txt")
err_Phi_0_Q_y_p_0d = read_file("err_obs_Phi_0_Q_y_+_zero_discord.txt")
exact_Phi_0_Q_y_p_0d = read_file("exact_Phi_0_Q_y_+_zero_discord.txt")
avg_Phi_0_Q_z_p_0d = read_file("avg_Phi_0_Q_z_+_zero_discord.txt")
err_Phi_0_Q_z_p_0d = read_file("err_obs_Phi_0_Q_z_+_zero_discord.txt")
exact_Phi_0_Q_z_p_0d = read_file("exact_Phi_0_Q_z_+_zero_discord.txt")

avg_Phi_0_Q_0_m_0d = read_file("avg_Phi_0_Q_0_-_zero_discord.txt")
err_Phi_0_Q_0_m_0d = read_file("err_obs_Phi_0_Q_0_-_zero_discord.txt")
exact_Phi_0_Q_0_m_0d = read_file("exact_Phi_0_Q_0_-_zero_discord.txt")
avg_Phi_0_Q_x_m_0d = read_file("avg_Phi_0_Q_x_-_zero_discord.txt")
err_Phi_0_Q_x_m_0d = read_file("err_obs_Phi_0_Q_x_-_zero_discord.txt")
exact_Phi_0_Q_x_m_0d = read_file("exact_Phi_0_Q_x_-_zero_discord.txt")
avg_Phi_0_Q_y_m_0d = read_file("avg_Phi_0_Q_y_-_zero_discord.txt")
err_Phi_0_Q_y_m_0d = read_file("err_obs_Phi_0_Q_y_-_zero_discord.txt")
exact_Phi_0_Q_y_m_0d = read_file("exact_Phi_0_Q_y_-_zero_discord.txt")
avg_Phi_0_Q_z_m_0d = read_file("avg_Phi_0_Q_z_-_zero_discord.txt")
err_Phi_0_Q_z_m_0d = read_file("err_obs_Phi_0_Q_z_-_zero_discord.txt")
exact_Phi_0_Q_z_m_0d = read_file("exact_Phi_0_Q_z_-_zero_discord.txt")

mu_0_p = (np.sqrt(3)+1)/2
mu_0_m = (np.sqrt(3)-1)/2
avg_Phi_0_tilde_Q_0_0d = mu_0_p*avg_Phi_0_Q_0_p_0d - mu_0_m*avg_Phi_0_Q_0_m_0d
err_Phi_0_tilde_Q_0_0d = .5*np.sqrt(err_Phi_0_Q_0_p_0d**2 + err_Phi_0_Q_0_m_0d**2)
exact_Phi_0_tilde_Q_0_0d = mu_0_p*exact_Phi_0_Q_0_p_0d - mu_0_m*exact_Phi_0_Q_0_m_0d
avg_Phi_0_tilde_Q_x_0d = .5*(avg_Phi_0_Q_x_p_0d - avg_Phi_0_Q_x_m_0d)
err_Phi_0_tilde_Q_x_0d = .5*np.sqrt(err_Phi_0_Q_x_p_0d**2 + err_Phi_0_Q_x_m_0d**2)
exact_Phi_0_tilde_Q_x_0d = .5*(exact_Phi_0_Q_x_p_0d - exact_Phi_0_Q_x_m_0d)
avg_Phi_0_tilde_Q_y_0d = .5*(avg_Phi_0_Q_y_p_0d - avg_Phi_0_Q_y_m_0d)
err_Phi_0_tilde_Q_y_0d = .5*np.sqrt(err_Phi_0_Q_y_p_0d**2 + err_Phi_0_Q_y_m_0d**2)
exact_Phi_0_tilde_Q_y_0d = .5*(exact_Phi_0_Q_y_p_0d - exact_Phi_0_Q_y_m_0d)
avg_Phi_0_tilde_Q_z_0d = .5*(avg_Phi_0_Q_z_p_0d - avg_Phi_0_Q_z_m_0d)
err_Phi_0_tilde_Q_z_0d = .5*np.sqrt(err_Phi_0_Q_z_p_0d**2 + err_Phi_0_Q_z_m_0d**2)
exact_Phi_0_tilde_Q_z_0d = .5*(exact_Phi_0_Q_z_p_0d - exact_Phi_0_Q_z_m_0d)

# Then: evol with Phi^z
avg_Phi_1_Q_0_p_0d = read_file("avg_Phi_1_Q_0_+_zero_discord.txt")
err_Phi_1_Q_0_p_0d = read_file("err_obs_Phi_1_Q_0_+_zero_discord.txt")
exact_Phi_1_Q_0_p_0d = read_file("exact_Phi_1_Q_0_+_zero_discord.txt")
avg_Phi_1_Q_x_p_0d = read_file("avg_Phi_1_Q_x_+_zero_discord.txt")
err_Phi_1_Q_x_p_0d = read_file("err_obs_Phi_1_Q_x_+_zero_discord.txt")
exact_Phi_1_Q_x_p_0d = read_file("exact_Phi_1_Q_x_+_zero_discord.txt")
avg_Phi_1_Q_y_p_0d = read_file("avg_Phi_1_Q_y_+_zero_discord.txt")
err_Phi_1_Q_y_p_0d = read_file("err_obs_Phi_1_Q_y_+_zero_discord.txt")
exact_Phi_1_Q_y_p_0d = read_file("exact_Phi_1_Q_y_+_zero_discord.txt")
avg_Phi_1_Q_z_p_0d = read_file("avg_Phi_1_Q_z_+_zero_discord.txt")
err_Phi_1_Q_z_p_0d = read_file("err_obs_Phi_1_Q_z_+_zero_discord.txt")
exact_Phi_1_Q_z_p_0d = read_file("exact_Phi_1_Q_z_+_zero_discord.txt")

avg_Phi_1_Q_0_m_0d = read_file("avg_Phi_1_Q_0_-_zero_discord.txt")
err_Phi_1_Q_0_m_0d = read_file("err_obs_Phi_1_Q_0_-_zero_discord.txt")
exact_Phi_1_Q_0_m_0d = read_file("exact_Phi_1_Q_0_-_zero_discord.txt")
avg_Phi_1_Q_x_m_0d = read_file("avg_Phi_1_Q_x_-_zero_discord.txt")
err_Phi_1_Q_x_m_0d = read_file("err_obs_Phi_1_Q_x_-_zero_discord.txt")
exact_Phi_1_Q_x_m_0d = read_file("exact_Phi_1_Q_x_-_zero_discord.txt")
avg_Phi_1_Q_y_m_0d = read_file("avg_Phi_1_Q_y_-_zero_discord.txt")
err_Phi_1_Q_y_m_0d = read_file("err_obs_Phi_1_Q_y_-_zero_discord.txt")
exact_Phi_1_Q_y_m_0d = read_file("exact_Phi_1_Q_y_-_zero_discord.txt")
avg_Phi_1_Q_z_m_0d = read_file("avg_Phi_1_Q_z_-_zero_discord.txt")
err_Phi_1_Q_z_m_0d = read_file("err_obs_Phi_1_Q_z_-_zero_discord.txt")
exact_Phi_1_Q_z_m_0d = read_file("exact_Phi_1_Q_z_-_zero_discord.txt")

avg_Phi_1_Q_0_0d = mu_0_p*avg_Phi_1_Q_0_p_0d - mu_0_m*avg_Phi_1_Q_0_m_0d
err_Phi_1_Q_0_0d = .5*np.sqrt(err_Phi_1_Q_0_p_0d**2 + err_Phi_1_Q_0_m_0d**2)
exact_Phi_1_Q_0_0d = mu_0_p*exact_Phi_1_Q_0_p_0d - mu_0_m*exact_Phi_1_Q_0_m_0d
avg_Phi_1_Q_x_0d = .5*(avg_Phi_1_Q_x_p_0d - avg_Phi_1_Q_x_m_0d)
err_Phi_1_Q_x_0d = .5*np.sqrt(err_Phi_1_Q_x_p_0d**2 + err_Phi_1_Q_x_m_0d**2)
exact_Phi_1_Q_x_0d = .5*(exact_Phi_1_Q_x_p_0d - exact_Phi_1_Q_x_m_0d)
avg_Phi_1_Q_y_0d = .5*(avg_Phi_1_Q_y_p_0d - avg_Phi_1_Q_y_m_0d)
err_Phi_1_Q_y_0d = .5*np.sqrt(err_Phi_1_Q_y_p_0d**2 + err_Phi_1_Q_y_m_0d**2)
exact_Phi_1_Q_y_0d = .5*(exact_Phi_1_Q_y_p_0d - exact_Phi_1_Q_y_m_0d)
avg_Phi_1_Q_z_0d = .5*(avg_Phi_1_Q_z_p_0d - avg_Phi_1_Q_z_m_0d)
err_Phi_1_Q_z_0d = .5*np.sqrt(err_Phi_1_Q_z_p_0d**2 + err_Phi_1_Q_z_m_0d**2)
exact_Phi_1_Q_z_0d = .5*(exact_Phi_1_Q_z_p_0d - exact_Phi_1_Q_z_m_0d)



params = open("params.txt")
tmax = float(params.readline())
t = np.linspace(0, tmax, avg_Phi_x_Q_z_p.shape[0])


fig, ax = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey = True)
fig.tight_layout()

#colors = sns.color_palette("muted", 4)
colors = ["goldenrod", "navy", "red", "darkgreen"]
markers = ['o', 'x', '^', 's', 'd']
markevery = int(len(avg_Phi_x_Q_z_p)/30)
markersize=5.5


ax[0].plot(t, exact_Phi_x_Q_0, color=colors[0], label=r'$\Phi^0_t[Q_0]$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_0, err_Phi_x_Q_0, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_x_Q_x, color=colors[1], label=r'$\Phi^x_t[Q_x]$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_x, err_Phi_x_Q_x, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_x_Q_y, color=colors[2], label=r'$\Phi^y_t[Q_y]$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_x_Q_y, err_Phi_x_Q_y, color=colors[2], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)
ax[0].plot(t, exact_Phi_z_Q_z, color=colors[3], label=r'$\Phi^z_t[Q_z]$', marker=markers[3], markersize=markersize, markevery=len(t))
ax[0].errorbar(t, avg_Phi_z_Q_z, err_Phi_z_Q_z, color=colors[3], marker=markers[3], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)





#colors = sns.color_palette("husl")
colors = ["navy", "red", "darkgreen"]

# Now rho_S and some repreparations
w_0 = 1
w_x = 1
w_y = 1
w_z = 1
exact_rho_S = w_0*exact_Phi_x_Q_0 + w_x*exact_Phi_x_Q_x + w_y*exact_Phi_x_Q_y + w_z*exact_Phi_z_Q_z
avg_rho_S = w_0*avg_Phi_x_Q_0 + w_x*avg_Phi_x_Q_x + w_y*avg_Phi_x_Q_y + w_z*avg_Phi_z_Q_z
err_rho_S = .25*np.sqrt(err_Phi_x_Q_0**2 + err_Phi_x_Q_x**2 + err_Phi_x_Q_y**2 + err_Phi_z_Q_z**2)
ax[1].plot(t, exact_rho_S, color=colors[0], label=r'$\mathcal{R}=\operatorname{id}$', marker=markers[0], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_S, err_rho_S, color=colors[0], marker=markers[0], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R 0 discord - p = 0.5
p = .5
avg_Phi_x_Q_0_0d = p*avg_Phi_0_tilde_Q_0_0d + (1-p)*avg_Phi_1_Q_0_0d
err_Phi_x_Q_0_0d = p*err_Phi_0_tilde_Q_0_0d + (1-p)*err_Phi_1_Q_0_0d
exact_Phi_x_Q_0_0d = p*exact_Phi_0_tilde_Q_0_0d + (1-p)*exact_Phi_1_Q_0_0d
avg_Phi_x_Q_x_0d = p*avg_Phi_0_tilde_Q_x_0d + (1-p)*avg_Phi_1_Q_x_0d
err_Phi_x_Q_x_0d = p*err_Phi_0_tilde_Q_x_0d + (1-p)*err_Phi_1_Q_x_0d
exact_Phi_x_Q_x_0d = p*exact_Phi_0_tilde_Q_x_0d + (1-p)*exact_Phi_1_Q_x_0d
avg_Phi_x_Q_y_0d = p*avg_Phi_0_tilde_Q_y_0d + (1-p)*avg_Phi_1_Q_y_0d
err_Phi_x_Q_y_0d = p*err_Phi_0_tilde_Q_y_0d + (1-p)*err_Phi_1_Q_y_0d
exact_Phi_x_Q_y_0d = p*exact_Phi_0_tilde_Q_y_0d + (1-p)*exact_Phi_1_Q_y_0d
avg_Phi_z_Q_z_0d = avg_Phi_1_Q_z_0d
err_Phi_z_Q_z_0d = err_Phi_1_Q_z_0d
exact_Phi_z_Q_z_0d = exact_Phi_1_Q_z_0d
exact_rho_S = w_0*exact_Phi_x_Q_0_0d + w_x*exact_Phi_x_Q_x_0d + w_y*exact_Phi_x_Q_y_0d + w_z*exact_Phi_z_Q_z_0d
avg_rho_S = w_0*avg_Phi_x_Q_0_0d + w_x*avg_Phi_x_Q_x_0d + w_y*avg_Phi_x_Q_y_0d + w_z*avg_Phi_z_Q_z_0d
err_rho_S = .25*np.sqrt(err_Phi_x_Q_0_0d**2 + err_Phi_x_Q_x_0d**2 + err_Phi_x_Q_y_0d**2 + err_Phi_z_Q_z_0d**2)
ax[1].plot(t, exact_rho_S, color=colors[1], label=r'$\mathcal{R}_0^{0.5}$', marker=markers[1], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_S, err_rho_S, color=colors[1], marker=markers[1], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R 0 discord - p = 0.9
p = .9
avg_Phi_x_Q_0_0d = p*avg_Phi_0_tilde_Q_0_0d + (1-p)*avg_Phi_1_Q_0_0d
err_Phi_x_Q_0_0d = p*err_Phi_0_tilde_Q_0_0d + (1-p)*err_Phi_1_Q_0_0d
exact_Phi_x_Q_0_0d = p*exact_Phi_0_tilde_Q_0_0d + (1-p)*exact_Phi_1_Q_0_0d
avg_Phi_x_Q_x_0d = p*avg_Phi_0_tilde_Q_x_0d + (1-p)*avg_Phi_1_Q_x_0d
err_Phi_x_Q_x_0d = p*err_Phi_0_tilde_Q_x_0d + (1-p)*err_Phi_1_Q_x_0d
exact_Phi_x_Q_x_0d = p*exact_Phi_0_tilde_Q_x_0d + (1-p)*exact_Phi_1_Q_x_0d
avg_Phi_x_Q_y_0d = p*avg_Phi_0_tilde_Q_y_0d + (1-p)*avg_Phi_1_Q_y_0d
err_Phi_x_Q_y_0d = p*err_Phi_0_tilde_Q_y_0d + (1-p)*err_Phi_1_Q_y_0d
exact_Phi_x_Q_y_0d = p*exact_Phi_0_tilde_Q_y_0d + (1-p)*exact_Phi_1_Q_y_0d
avg_Phi_z_Q_z_0d = avg_Phi_1_Q_z_0d
err_Phi_z_Q_z_0d = err_Phi_1_Q_z_0d
exact_Phi_z_Q_z_0d = exact_Phi_1_Q_z_0d
exact_rho_S = w_0*exact_Phi_x_Q_0_0d + w_x*exact_Phi_x_Q_x_0d + w_y*exact_Phi_x_Q_y_0d + w_z*exact_Phi_z_Q_z_0d
avg_rho_S = w_0*avg_Phi_x_Q_0_0d + w_x*avg_Phi_x_Q_x_0d + w_y*avg_Phi_x_Q_y_0d + w_z*avg_Phi_z_Q_z_0d
err_rho_S = .25*np.sqrt(err_Phi_x_Q_0_0d**2 + err_Phi_x_Q_x_0d**2 + err_Phi_x_Q_y_0d**2 + err_Phi_z_Q_z_0d**2)
ax[1].plot(t, exact_rho_S, '--', color=colors[1], label=r'$\mathcal{R}_0^{0.9}$', marker=markers[2], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_rho_S, err_rho_S, color=colors[1], marker=markers[2], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)

# R factorized
r00 = 1
r0x = 1
r0y = 1
r0z = 1
exact_R1 = w_0*(r00*exact_Phi_x_Q_0 + r0x*exact_Phi_x_Q_y + r0y*exact_Phi_x_Q_y + r0z*exact_Phi_x_Q_z)
avg_R1 = w_0*(r00*avg_Phi_x_Q_0 + r0x*avg_Phi_x_Q_x + r0y*avg_Phi_x_Q_y + r0z*avg_Phi_x_Q_z)
ax[1].plot(t, exact_R1, color=colors[2], label=r'$\mathcal{R}_{\text{f}}$', marker=markers[3], markersize=markersize, markevery=len(t))
ax[1].errorbar(t, avg_R1, err_rho_S, color=colors[2], marker=markers[3], markevery=markevery, errorevery=markevery, markersize=markersize, linewidth=0, elinewidth=1)



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