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
path1 = "Phobos/Spherical/MatrixSolution.out"
path2 = "Phobos/Rotated/MatrixSolution.out"
path3 = "Phobos/Sensitivity/SensitivityKernel.out"
path4 = "Phobos/Sensitivity/SensitivityKernelSlice.out"


dphobos = np.loadtxt(path0, delimiter=";")
# dphobos2 = np.loadtxt(path1,delimiter=";")
# dphobos3 = np.loadtxt(path2,delimiter=";")
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


colnum = 1
maxval = 0
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            maxval = max(maxval,max(abs(dphobos2[:,colnum])))
            colnum = colnum + 1

row2,col2 = dphobos3.shape
numpoints = int(col2/4)
x = np.zeros(row2 * numpoints)
y = np.zeros(row2 * numpoints)
sensreal = np.zeros(row2 * numpoints)
sensimag = np.zeros(row2 * numpoints)
realmax = 0
imagmax = 0
idx = 0
for idxr in range (0,row2):
    for idxp in range (0,numpoints):
        x[idx] = dphobos3[idxr,idxp * 4] * np.cos(dphobos3[idxr,idxp * 4 + 1] )
        y[idx] = dphobos3[idxr,idxp * 4] * np.sin(dphobos3[idxr,idxp * 4 + 1] )
        sensreal[idx] = dphobos3[idxr,idxp * 4 + 2]
        sensimag[idx] = dphobos3[idxr,idxp * 4 + 3]
        # if(max(realmax,sensreal[idx]) > realmax):
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

# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# tpc = ax.tripcolor(triang, sensreal, shading='flat')
# fig.colorbar(tpc)
# ax.set_title('tripcolor of Delaunay triangulation, flat shading')

# plt.tick_params(left = False, right = False , labelleft = False , 
#                 labelbottom = False, bottom = False) 


# plt.rcParams['text.usetex'] = True

fig2, ax2 = plt.subplots(1,1)
ax2.set_aspect('equal')
cntr2 = ax2.tricontourf(x, y, sensreal, levels=400, cmap="jet")
# tpc = ax2.tripcolor(triang, sensreal,cmap='jet', shading='gouraud')
cbar = fig2.colorbar(cntr2)
cbar.ax.tick_params(labelsize=15)
ax2.set_title(r'$Y_{2,0}$ real component', fontsize=20)
ax2.set_yticklabels([])
ax2.set_xticklabels([])

# ax2[1].set_aspect('equal')
# tpc1 = ax2[1].tripcolor(triang, sensimag, shading='gouraud')
# fig2.colorbar(tpc1)
# ax2[1].set_title(r'$Y_{2,0}$ imag')
# ax2[1].set_yticklabels([])
# ax2[1].set_xticklabels([])


# fig2, ax2 = plt.subplots(1,2)
# ax2[0].set_aspect('equal')
# tpc = ax2[0].tripcolor(triang, sensreal, shading='gouraud')
# fig2.colorbar(tpc)
# ax2[0].set_title(r'$Y_{2,0}$ real')
# ax2[0].set_yticklabels([])
# ax2[0].set_xticklabels([])

# ax2[1].set_aspect('equal')
# tpc1 = ax2[1].tripcolor(triang, sensimag, shading='gouraud')
# fig2.colorbar(tpc1)
# ax2[1].set_title(r'$Y_{2,0}$ imag')
# ax2[1].set_yticklabels([])
# ax2[1].set_xticklabels([])

plt.show()
# potential 
# colnum = 1
# for idxl in range(0,l+1):
#     for idxm in range (-idxl,idxl+1):
#         tmpmax= max(abs(dphobos2[:,colnum]))
#         if tmpmax/maxval > 0.01:
#             ax.plot(dphobos2[:,0],dphobos2[:,colnum],linewidth=3, label="Re" + str(idxl) + str(idxm))
#         else:
#             ax.plot(dphobos2[:,0],dphobos2[:,colnum],linewidth=3)
#         colnum = colnum + 1
#         tmpmax= max(abs(dphobos2[:,colnum]))
#         if tmpmax/maxval > 0.01:
#             ax.plot(dphobos2[:,0],dphobos2[:,colnum],linewidth=3, label="Im" + str(idxl) + str(idxm))
#         else:
#             ax.plot(dphobos2[:,0],dphobos2[:,colnum],linewidth=3)
#         colnum = colnum + 1

# ax.set_ylabel("Sensitivity Kernel", fontsize = BIGGER_SIZE)
# ax.set_xlabel("Referential radius", fontsize = BIGGER_SIZE)
# ax.legend(fontsize=BIGGER_SIZE)
# plt.show()