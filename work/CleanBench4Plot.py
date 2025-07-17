#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("Lagrange.out", delimiter=";")
path_to_folder = "Bench4"
path1 = path_to_folder + "/ExactSolution.out"
# path2 = path_to_folder + "/IntegralSolution.out"
path3 = path_to_folder + "/MatrixSolution.out"
# d = np.loadtxt("cleanbench1.out", delimiter=";")
# d0 = np.loadtxt(path1, delimiter=";")
# d1 = np.loadtxt(path2, delimiter=";")
d2 = np.loadtxt(path3, delimiter=";")
d0 = np.loadtxt(path1, delimiter=";")
# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(2,1)
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


# potential 
ax[0].plot(d0[:,0],d0[:,1], "k", linewidth=3)
# ax[0].plot(d1[:, 0], d1[:, 3], "b", linewidth=3)
ax[0].plot(d2[:, 0], d2[:, 3], "b", linewidth=3)
ax[0].plot(d2[:, 5], d2[:, 8], "r--", linewidth=3)
# ax[0].plot(d1[:, 0], d1[:, 3], "ko")
ax[0].plot(d0[:,0],d0[:,1], "ko")
ax[0].legend([ 'Exact','Matrix-free pseudospectral','Nodes'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0].set_ylabel("Potential", fontsize = BIGGER_SIZE)

# error in potential
# interpolation
# f = sc.interpolate.CubicSpline(d1[:,0],d1[:,3])
# f2 = sc.interpolate.CubicSpline(d2[:,0],d2[:,3])

# x_new = np.linspace(0,d1[-1,0],120)
# y1 = f(x_new)
# y2 = f2(x_new)

# ax[1].plot(x_new, abs(y1 - y2)/abs(max(y1)) * 100.0, "b", linewidth=3)
ax[1].plot(d0[:,0],abs(d0[:,1] - d2[:,3])/abs(max(d0[:,1])) * 100.0, "b", linewidth=3)
# ax[1].plot(x_new,y1,"b",linewidth=3)
# ax[1].plot(x_new,y2,"r--",linewidth=3)
# ax[1].plot(d1[:,0], abs(d2[:,1] - d1[:,1])/abs(max(d1[:,1])) * 100.0, "r", linewidth=3)
ax[1].legend(['Relative difference (%)'],fontsize = MEDIUM_SIZE)
ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)


# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()

