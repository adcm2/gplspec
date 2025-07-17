#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("Lagrange.out", delimiter=";")
path_to_folder = "Bench6/lm"
# path1 = path_to_folder + "/IntegralSolution.out"
# path2 = path_to_folder + "/IntegralSolution.out"
# path2 = path_to_folder + "/MatrixSolution.out"
path3 = path_to_folder + "/SensitivityKernel.out"

# d0 = np.loadtxt(path1, delimiter=";")
# d2 = np.loadtxt(path2, delimiter=";")
d3 = np.loadtxt(path3, delimiter=";")

# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(1,1)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

row, col = d3.shape
l = int(math.sqrt((col - 1)/2)-1)

mymaxval = 0
colnum = 1
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        maxtmp = max(abs(d3[:,colnum]))
        mymaxval = max(maxtmp,mymaxval)
        colnum = colnum + 1
        maxtmp = max(abs(d3[:,colnum]))
        mymaxval = max(maxtmp,mymaxval)
        colnum = colnum + 1


# potential 
colnum = 1
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        tmpmax= max(abs(d3[:,colnum]))
        if tmpmax/mymaxval > 0.01:
            ax.plot(d3[:,0],d3[:,colnum],linewidth=3, label="Re" + str(idxl) + str(idxm))
        else:
            ax.plot(d3[:,0],d3[:,colnum],linewidth=3)
        colnum = colnum + 1
        tmpmax= max(abs(d3[:,colnum]))
        if tmpmax/mymaxval > 0.01:
            ax.plot(d3[:,0],d3[:,colnum],linewidth=3, label="Im" + str(idxl) + str(idxm))
        else:
            ax.plot(d3[:,0],d3[:,colnum],linewidth=3)
        colnum = colnum + 1

# ax.set_yscale('log')
ax.tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax.yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax.set_ylabel("Sensitivity Kernel", fontsize = BIGGER_SIZE)
ax.set_xlabel("Radius", fontsize = BIGGER_SIZE)
ax.legend(fontsize=BIGGER_SIZE)

plt.show()

