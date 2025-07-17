#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("Lagrange.out", delimiter=";")
path_to_folder = "Bench6/lm"
path1 = path_to_folder + "/IntegralSolution.out"
# path2 = path_to_folder + "/IntegralSolution.out"
path2 = path_to_folder + "/MatrixSolution.out"
path3 = path_to_folder + "/SensitivityKernel.out"
# d = np.loadtxt("cleanbench1.out", delimiter=";")
# d0 = np.loadtxt(path1, delimiter=";")
# d1 = np.loadtxt(path2, delimiter=";")
d0 = np.loadtxt(path1, delimiter=";")
d2 = np.loadtxt(path2, delimiter=";")
d3 = np.loadtxt(path3, delimiter=";")

# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(2,3)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

row, col = d0.shape
l = int(math.sqrt((col - 1)/2)-1)


# potential 
ax[0,0].plot(d0[:,0],d0[:,1], "k", linewidth=3)
ax[0,0].plot(d2[:, 0], d2[:, 1], "b", linewidth=3)
ax[0,0].plot(d2[:, 0], d2[:, 2], "r--", linewidth=3)
# ax[0,0].plot(d1[:, 0], d1[:, 3], "ko")
ax[0,0].plot(d0[:,0],d0[:,1], "ko")
ax[0,0].legend([ 'Exact','MFP (0,0) REAL', 'MFP (0,0) IMAG','Nodes'], fontsize = MEDIUM_SIZE)
ax[0,0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0,0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0,0].set_ylabel("Potential", fontsize = BIGGER_SIZE)

# error in potential
n = 5
m1 = 1
m2 = 500
c1 = n*m1
c11 = n*m1 + 3
c2 = n*m2
c21 = n*m2+3

ax[1,0].plot(d0[:,0],abs(d0[:,1] - d2[:,1])/abs(max(d0[:,1])) * 100.0, "b", linewidth=3)
ax[1,0].legend(['(0,0)'],fontsize = MEDIUM_SIZE)
ax[1,0].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1,0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1,0].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
ax[1,0].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)

# potential 
c1 = 9
c2 = 10
ax[0,1].plot(d0[:,0],d0[:,c1], "k", linewidth=3)
ax[0,1].plot(d0[:,0],d0[:,c2], "g", linewidth=3)
ax[0,1].plot(d2[:, 0], d2[:, c1], "b--", linewidth=3)
ax[0,1].plot(d2[:, 0], d2[:, c2], "r--", linewidth=3)
ax[0,1].plot(d0[:,0],d0[:,c1], "ko")
ax[0,1].legend([ 'INT (1,-1) REAL','INT (1,-1) IMAG','MFP (1,-1) REAL', 'MFP (1,-1) IMAG','Nodes'], fontsize = MEDIUM_SIZE)
ax[0,1].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0,1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0,1].set_ylabel("Potential", fontsize = BIGGER_SIZE)

ax[1,1].plot(d0[:,0],abs(d0[:,3] - d2[:,3])/max(abs(d0[:,3])) * 100.0, "b", linewidth=3)
ax[1,1].legend(['(1,-1)'],fontsize = MEDIUM_SIZE)
ax[1,1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1,1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1,1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
ax[1,1].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)


colnum = 1
maxval = 0
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            maxval = max(maxval,max(abs(d0[:,colnum])))
            colnum = colnum + 1

# potential 
colnum = 1
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            ax[0,2].plot(d0[:,0],abs(d0[:,colnum] - d2[:,colnum])/maxval * 100.0, linewidth=3)
            ax[1,2].plot(d0[:,0],d3[:,colnum],linewidth=3)
            colnum = colnum + 1

ax[0,2].set_yscale('log')
ax[0,2].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0,2].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0,2].set_ylabel("Rel. Diff (%)", fontsize = BIGGER_SIZE)

# ax[1,2].plot(d0[:,0],abs(d0[:,5] - d2[:,5])/max(abs(d0[:,5])) * 100.0, "b", linewidth=3)

# # ax[1,1].plot(d0[:,c1],abs(d0[:,c11] - d2[:,c11])/abs(max(d0[:,c11])) * 100.0, "r", linewidth=3)
# # ax[1,1].plot(d0[:,c2],abs(d0[:,c21] - d2[:,c21])/abs(max(d0[:,c21])) * 100.0, "k", linewidth=3)
# # # ax[1,1].plot(x_new,y1,"b",linewidth=3)
# # ax[1,1].plot(x_new,y2,"r--",linewidth=3)
# # ax[1,1].plot(d1[:,0], abs(d2[:,1] - d1[:,1])/abs(max(d1[:,1])) * 100.0, "r", linewidth=3)
# # ax[1,2].legend(['Angle 1', 'Angle 2', 'Angle 3'],fontsize = MEDIUM_SIZE)
# ax[1,2].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
# t = ax[1,2].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[1,2].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
# ax[1,2].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)

# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()

