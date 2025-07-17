#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
import math

# import data
path0 = "Phobos/Cartesian/MatrixSolution.out"
path1 = "Phobos/Spherical/MatrixSolution.out"
path2 = "Phobos/Rotated/MatrixSolution.out"
path3 = "Phobos/Sensitivity/SensitivityKernel.out"

dphobos = np.loadtxt(path0, delimiter=";")
dphobos2 = np.loadtxt(path1,delimiter=";")
dphobos3 = np.loadtxt(path2,delimiter=";")
dphobos4 = np.loadtxt(path3,delimiter=";")

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

# find l, etc
row, col = dphobos.shape
numl = int((col-1)/5)
l = int((math.sqrt(1 + 2*numl)-1)/2)
print("Number of rows: ", row)
print("Number of cols: ", col)
print("l: ", l)

# find index of outer radius
# this assumes that the length norm is the referential radius of the planet
idxnextlayer = 0
for idx in range(0,row-1):
    if (dphobos[idx+1,0] - dphobos[idx,0] == 1):
        idxnextlayer = idx +1
lenln = l+1  

fig, ax = plt.subplots(1,1)
colnum = 1
maxval = 0
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            maxval = max(maxval,max(abs(dphobos4[:,colnum])))
            colnum = colnum + 1


# potential 
colnum = 1
lmax = min(11,l+1)
for idxl in range(0,lmax):
    for idxm in range (-idxl,idxl+1):
        tmpmax= max(abs(dphobos4[:,colnum]))
        if tmpmax/maxval > 0.01:
            ax.plot(dphobos4[:,0],dphobos4[:,colnum],linewidth=3, label="Re" + str(idxl) + str(idxm))
        else:
            ax.plot(dphobos4[:,0],dphobos4[:,colnum],linewidth=3)
        colnum = colnum + 1
        tmpmax= max(abs(dphobos4[:,colnum]))
        if tmpmax/maxval > 0.01:
            ax.plot(dphobos4[:,0],dphobos4[:,colnum],linewidth=3, label="Im" + str(idxl) + str(idxm))
        else:
            ax.plot(dphobos4[:,0],dphobos4[:,colnum],linewidth=3)
        colnum = colnum + 1

ax.set_ylabel("Sensitivity Kernel", fontsize = BIGGER_SIZE)
ax.set_xlabel("Referential radius", fontsize = BIGGER_SIZE)
ax.legend(fontsize=BIGGER_SIZE)
plt.show()