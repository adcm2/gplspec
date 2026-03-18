#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
import math

import matplotlib.tri as tri

# import data
path0 = "Phobos/Cartesian/MatrixSolution.out"
path3 = "Phobos/Sensitivity/SensitivityKernel.out"
path4 = "Phobos/Sensitivity/SensitivityKernelSlice.out"

dphobos = np.loadtxt(path0, delimiter=";")
dphobos2 = np.loadtxt(path3,delimiter=";")
dphobos3 = np.loadtxt(path4,delimiter=";")

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
xlen = idxnextlayer

colnum = 1
maxval = 0
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            maxval = max(maxval,max(abs(dphobos2[:,colnum])))
            colnum = colnum + 1

row2,col2 = dphobos3.shape
numpoints = int(col2/4)
x = np.zeros(xlen * numpoints)
y = np.zeros(xlen * numpoints)
sensreal = np.zeros(xlen * numpoints)
sensimag = np.zeros(xlen * numpoints)
realmax = 0
imagmax = 0
idx = 0
for idxr in range (0,xlen):
    for idxp in range (0,numpoints):
        x[idx] = dphobos3[idxr,idxp * 4] * np.cos(dphobos3[idxr,idxp * 4 + 1] )
        y[idx] = dphobos3[idxr,idxp * 4] * np.sin(dphobos3[idxr,idxp * 4 + 1] )
        sensreal[idx] = dphobos3[idxr,idxp * 4 + 2]
        sensimag[idx] = dphobos3[idxr,idxp * 4 + 3]
        realmax = max(realmax,abs(sensreal[idx]))
        imagmax = max(imagmax,abs(sensimag[idx]))
        idx = idx + 1
print("Max value: ", realmax)
sensreal = sensreal/realmax
sensimag = sensimag/realmax
print("Number of rows: ", row2)
print("Number of cols: ", col2)
triang = tri.Triangulation(x,y)
minrad = 0.001
triang.set_mask(np.hypot(x[triang.triangles].mean(axis=1),
                         y[triang.triangles].mean(axis=1))
                < minrad)

fig2, ax2 = plt.subplots(1,1)
ax2.set_aspect('equal')
cntr2 = ax2.tricontourf(x, y, sensreal, levels=400, cmap="jet")
cbar = fig2.colorbar(cntr2,format='%.1f')
cbar.ax.set_ylabel('Sensitivity kernel (N.D.)', labelpad=30, rotation=270, fontsize=20)
cbar.ax.tick_params(axis='both', which='major', labelsize=15)
ax2.set_xlabel('x (m)', fontsize=20)
ax2.set_ylabel('y (m)', fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.set_title(r'$Y_{2,0}$ real component', fontsize=20)

plt.show()