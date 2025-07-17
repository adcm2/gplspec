#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

path_to_folder = "Bench7"
path1 = path_to_folder + "/zeroradius/MatrixSolution.out"
d = np.loadtxt(path1, delimiter=";")

# initialise plot
fig, ax = plt.subplots(1,2)
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

errd = abs(d[:,2] - d[:,1])/max(abs(d[:,1]))
nidx = 3
b = errd[nidx]/(d[nidx,0]*d[nidx,0])
prederr = b * d[:,0] * d[:,0]

# potential 
ax[0].plot(d[:,0],d[:,2], "b", linewidth=3)
ax[0].plot(d[:,0],d[:,1], "r--", linewidth=3)
ax[0].legend([ 'Perturbation','Exact'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0].set_ylabel("Potential at centre of planet", fontsize = BIGGER_SIZE)
ax[0].set_xlabel("Perturbation h", fontsize = BIGGER_SIZE)


ax[1].plot(d[:,0], errd * 100.0, "b", linewidth=3)
ax[1].plot(d[:,0], prederr * 100.0, "r--", linewidth=3)
ax[1].legend([ 'Error','Quadratic error prediction'], fontsize = MEDIUM_SIZE)
ax[1].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1].set_xlabel("Perturbation h", fontsize = BIGGER_SIZE)
ax[1].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)
ax[1].set_xscale('log')
ax[1].set_yscale('log')


# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()

