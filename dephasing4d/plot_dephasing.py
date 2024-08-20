import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

filein = open("tmax.txt")
tmax = float(filein.readline())

# Function to read data from a file
def read_file(filename):
    with open(filename, 'r') as file:
        data = np.loadtxt(file)
    return data

def read_file_columns(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            for i, number in enumerate(numbers):
                if len(data) <= i:
                    data.append([])
                data[i].append(float(number))  # Convert to float or int as needed
    return data

exact_p = read_file("exact_p.txt")
exact_m = read_file("exact_m.txt")
exact_0 = read_file("exact_0.txt")

MCWF_p = read_file("MCWF_p.txt")
MCWF_m = read_file("MCWF_m.txt")
MCWF_traj_p = read_file_columns("MCWF_traj_p.txt")
MCWF_traj_m = read_file_columns("MCWF_traj_m.txt")

QSD_p = read_file("QSD_p.txt")
QSD_m = read_file("QSD_m.txt")
QSD_traj_p = read_file_columns("QSD_traj_p.txt")
QSD_traj_m = read_file_columns("QSD_traj_m.txt")

rates = read_file_columns("rates.txt")

Npoints = len(exact_0)
t = np.linspace(0,tmax,Npoints)
markevery = int(Npoints/20)
Ntraj = len(MCWF_traj_p)

fig, ax = plt.subplots(1,3, figsize=(16,4), sharex=True, sharey=True)
axx = ax[2].inset_axes([.1,.1,.4,.4])
fig.tight_layout()

ax[0].plot(t,exact_p,color='black')
ax[0].plot(t,MCWF_p, marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')
ax[0].plot(t,QSD_p, marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0, fillstyle='none')
for i in range(Ntraj):
    ax[0].plot(t, MCWF_traj_p[i], alpha=.1, color='red')
    ax[0].plot(t, QSD_traj_p[i], alpha=.1, color='blue')

ax[1].plot(t,exact_m,color='black')
ax[1].plot(t,MCWF_m, marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')
ax[1].plot(t,QSD_m, marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0, fillstyle='none')
for i in range(Ntraj):
    ax[1].plot(t, MCWF_traj_m[i], alpha=.1, color='red')
    ax[1].plot(t, QSD_traj_m[i], alpha=.1, color='blue')

ax[2].plot(t, .5*exact_0, color='black')
ax[2].plot(t, .5*(MCWF_p - MCWF_m), marker='o', markersize=6, color='red', label="MCWF", markevery=markevery, linewidth=0, fillstyle='none')
ax[2].plot(t, .5*(QSD_p - QSD_m), marker='x', markersize=6, color='blue', label="QSD", markevery=markevery, linewidth=0, fillstyle='none')

for i in range(3):
    axx.plot(t, rates[i])
axx.axhline(0, color='black', linewidth=.5)

fontSize = 14
for i in range(3):
    ax[i].set_xlabel(r'$t$', fontsize=fontSize)
ax[0].legend(loc="lower right", fontsize=fontSize)
ax[1].legend(loc="upper right", fontsize=fontSize)
ax[2].legend(loc="lower right", fontsize=fontSize)
#axx.legend(loc = "upper right", fontsize=fontSize-1)
ax[0].set_title(r"$Q_\alpha^+$", fontsize=fontSize)
ax[1].set_title(r"$Q_\alpha^-$", fontsize=fontSize)
ax[2].set_title(r"$Q_\alpha = \frac{1}{2}(Q_\alpha^+ - Q_\alpha^-)$", fontsize=fontSize)
axx.set_title("Rates", fontsize=fontSize-1)
ax[0].set_ylabel(r"Re$\langle 0\vert X(t)\vert 1\rangle$", fontsize=fontSize)

for i in range(3):
    ax[i].tick_params(axis='both', which='major', labelsize=fontSize-1)
axx.tick_params(axis='both', which='major', labelsize=fontSize-2)

plt.show()