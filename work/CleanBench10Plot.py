import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

path_to_folder = "Bench10"
path1 = path_to_folder + "/exact/MatrixSolution.out"
path2 = path_to_folder + "/perturbed/MatrixSolution.out"
path3 = path_to_folder + "/unperturbed/MatrixSolution.out"
path4 = path_to_folder + "/perturbed2/MatrixSolution.out"
dexact = np.loadtxt(path1, delimiter=";")
dp = np.loadtxt(path2, delimiter=";")
dup = np.loadtxt(path3, delimiter=";")
dp2 = np.loadtxt(path4, delimiter=";")

fig, ax = plt.subplots(2,2)
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

n = 5
m1 = 200
m2 = 10
c1 = n*m1
c11 = n*m1 + 3
c2 = n*m2
c21 = n*m2+3

ax[0,0].plot(dexact[:,0],dexact[:,3], "k", linewidth=3)
ax[0,0].plot(dp[:, 0], dp[:, 3], "b", linewidth=3)
ax[0,0].plot(dup[:, 0], dup[:, 3], "r--", linewidth=3)
ax[0,0].plot(dp2[:, 0], dp2[:, 3], "g-.", linewidth=3)
ax[0,0].legend([ 'Exact','Classical perturbation theory', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
ax[0,0].tick_params(axis='both', which='major', labelsize=15)
t = ax[0,0].yaxis.get_offset_text()
t.set_size(15)
ax[0,0].set_ylabel("Potential (m${}^2/$s${}^2$)", fontsize = 20, labelpad = 15)

ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp[:,3])/abs(max(dexact[:,3])) * 100.0, "b", linewidth=3)
ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dup[:,3])/abs(max(dexact[:,3])) * 100.0, "r", linewidth=3)
ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp2[:,3])/abs(max(dexact[:,3])) * 100.0, "g-.", linewidth=3)
ax[1,0].legend(['Classical perturbation theory', 'Unperturbed solution', 'Perturb a'],fontsize = MEDIUM_SIZE, loc='upper left')
ax[1,0].tick_params(axis='both', which='major', labelsize=15)
t = ax[1,0].yaxis.get_offset_text()
t.set_size(15)
ax[1,0].set_xlabel("Radius (m)", fontsize = 20, labelpad = 15)
ax[1,0].set_ylabel("Relative difference (%)", fontsize = 20, labelpad = 15)

ax[0,1].plot(dexact[:,c1],dexact[:,c11], "k", linewidth=3)
ax[0,1].plot(dp[:, c1], dp[:, c11], "b", linewidth=3)
ax[0,1].plot(dup[:, c1], dup[:, c11], "r--", linewidth=3)
ax[0,1].plot(dp2[:, c1], dp2[:, c11], "g-.", linewidth=3)
ax[0,1].legend([ 'Exact','Classical perturbation theory', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
ax[0,1].tick_params(axis='both', which='major', labelsize=15)
t = ax[0,1].yaxis.get_offset_text()
t.set_size(15)
ax[0,1].set_ylabel("Potential (m${}^2/$s${}^2$)", fontsize = 20, labelpad = 15)

ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp[:,c11])/abs(max(dexact[:,c11])) * 100.0, "b", linewidth=3)
ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dup[:,c11])/abs(max(dexact[:,c11])) * 100.0, "r", linewidth=3)
ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp2[:,c11])/abs(max(dexact[:,c11])) * 100.0, "g-.", linewidth=3)
ax[1,1].legend(['Classical perturbation theory', 'Unperturbed solution','Perturb a'],fontsize = MEDIUM_SIZE, loc='upper left')
ax[1,1].tick_params(axis='both', which='major', labelsize=15)
t = ax[1,1].yaxis.get_offset_text()
t.set_size(15)
ax[1,1].set_xlabel("Radius (m)", fontsize = 20, labelpad = 15)
ax[1,1].set_ylabel("Relative difference (%)", fontsize = 20, labelpad = 15)

fig.align_ylabels(ax)

plt.show()

