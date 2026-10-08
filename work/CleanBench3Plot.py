import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

path_to_folder = "Bench3"
path1 = path_to_folder + "/ExactSolution.out"
path3 = path_to_folder + "/MatrixSolution.out"
d2 = np.loadtxt(path3, delimiter=";")
d0 = np.loadtxt(path1, delimiter=";")

# initialise plot
fig, ax = plt.subplots(2,1)
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

# potential 
ax[0].plot(d0[:,0],d0[:,1], "k", linewidth=3)
ax[0].plot(d2[:, 0], d2[:, 3], "r--", linewidth=3)
ax[0].plot(d0[:,0],d0[:,1], "ko")
ax[0].legend([ 'Exact','Matrix-free pseudospectral','Nodes'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=15)
t = ax[0].yaxis.get_offset_text()
t.set_size(15)
ax[0].set_ylabel("Potential (m${}^2/$s${}^2$)", fontsize = 20, labelpad = 15)

# error in potential
ax[1].plot(d0[:,0],abs(d0[:,1] - d2[:,3])/abs(max(d0[:,1])) * 100.0, "b", linewidth=3)
ax[1].tick_params(axis='both', which='major', labelsize=15)
t = ax[1].yaxis.get_offset_text()
t.set_size(15)
ax[1].set_xlabel("Radius (m)", fontsize = 20, labelpad = 15)
ax[1].set_ylabel("Relative difference (%)", fontsize = 20, labelpad = 15)

fig.align_ylabels(ax)

plt.show()

