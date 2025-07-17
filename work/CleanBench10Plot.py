#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("Lagrange.out", delimiter=";")
path_to_folder = "Bench10"
path1 = path_to_folder + "/exact/MatrixSolution.out"
# path2 = path_to_folder + "/IntegralSolution.out"
path2 = path_to_folder + "/perturbed/MatrixSolution.out"
path3 = path_to_folder + "/unperturbed/MatrixSolution.out"
path4 = path_to_folder + "/perturbed2/MatrixSolution.out"
# d = np.loadtxt("cleanbench1.out", delimiter=";")
# d0 = np.loadtxt(path1, delimiter=";")
# d1 = np.loadtxt(path2, delimiter=";")
dexact = np.loadtxt(path1, delimiter=";")
dp = np.loadtxt(path2, delimiter=";")
dup = np.loadtxt(path3, delimiter=";")
dp2 = np.loadtxt(path4, delimiter=";")

# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(2,2)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# row, col = d1.shape
# l = int(math.sqrt((col - 1)/2)-1)
# print("Number of cols: ", col)
# print("Check: ", math.sqrt((col - 1)/2))
# print("l: ", l)
# print("First few values d0", d0[0:5,0])
# print("First few values d2", d2[0:5,0])

n = 5
m1 = 200
m2 = 10
c1 = n*m1
c11 = n*m1 + 3
c2 = n*m2
c21 = n*m2+3

# potential 
ax[0,0].plot(dexact[:,0],dexact[:,3], "k", linewidth=3)
ax[0,0].plot(dp[:, 0], dp[:, 3], "b", linewidth=3)
ax[0,0].plot(dup[:, 0], dup[:, 3], "r--", linewidth=3)
ax[0,0].plot(dp2[:, 0], dp2[:, 3], "g-.", linewidth=3)
ax[0,0].legend([ 'Exact','Classical perturbation theory', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
ax[0,0].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[0,0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0,0].set_ylabel("Potential", fontsize = BIGGER_SIZE)

ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp[:,3])/abs(max(dexact[:,3])) * 100.0, "b", linewidth=3)
ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dup[:,3])/abs(max(dexact[:,3])) * 100.0, "r", linewidth=3)
ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp2[:,3])/abs(max(dexact[:,3])) * 100.0, "g-.", linewidth=3)
ax[1,0].legend(['Classical perturbation theory', 'Unperturbed solution', 'Perturb a'],fontsize = MEDIUM_SIZE)
ax[1,0].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1,0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1,0].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
ax[1,0].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)


ax[0,1].plot(dexact[:,c1],dexact[:,c11], "k", linewidth=3)
ax[0,1].plot(dp[:, c1], dp[:, c11], "b", linewidth=3)
ax[0,1].plot(dup[:, c1], dup[:, c11], "r--", linewidth=3)
ax[0,1].plot(dp2[:, c1], dp2[:, c11], "g-.", linewidth=3)
ax[0,1].legend([ 'Exact','Classical perturbation theory', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
ax[0,1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[0,1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0,1].set_ylabel("Potential", fontsize = BIGGER_SIZE)

ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp[:,c11])/abs(max(dexact[:,c11])) * 100.0, "b", linewidth=3)
ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dup[:,c11])/abs(max(dexact[:,c11])) * 100.0, "r", linewidth=3)
ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp2[:,c11])/abs(max(dexact[:,c11])) * 100.0, "g-.", linewidth=3)
ax[1,1].legend(['Classical perturbation theory', 'Unperturbed solution','Perturb a'],fontsize = MEDIUM_SIZE)
ax[1,1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1,1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1,1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
ax[1,1].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)
# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()

