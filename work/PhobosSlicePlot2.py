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
path5 = "Phobos/Spherical/MatrixSolutionSlice.out"
path6 = "Phobos/Spherical/MatrixSolutionRotated.out"
path6 = "Phobos/Spherical/MatrixSolutionReferentialRotated.out"
path6 = "Phobos/Spherical/MatrixSolutionReferential.out"
path6 = "Phobos/Spherical/MatrixSolutionPhysical.out"



dphobos = np.loadtxt(path0, delimiter=";")
# dphobos2 = np.loadtxt(path1,delimiter=";")
# dphobos3 = np.loadtxt(path2,delimiter=";")
dphobos2 = np.loadtxt(path3,delimiter=";")
# dphobos3 = np.loadtxt(path4,delimiter=";")
dphobos3 = np.loadtxt(path5,delimiter=";")
dphobos4 = np.loadtxt(path6,delimiter=";")

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

print("next layer: ",idxnextlayer)
colnum = 1
maxval = 0
for idxl in range(0,l+1):
    for idxm in range (-idxl,idxl+1):
        for idxn in range (0,2):
            maxval = max(maxval,max(abs(dphobos2[:,colnum])))
            colnum = colnum + 1


# number of longitude points
nlong = 2 * l
idxlow = 0
# bool evenn = false

ncolat = l + 1
if (ncolat%2 == 0):
    idxlow = int(ncolat/2 - 1)
else:
    idxlow = int((ncolat - 1)/2)

row2,col2 = dphobos3.shape
numpoints = int(col2/4)
print("numpoints: ",numpoints)
print("nlong: ",nlong)
x = np.zeros(row2 * numpoints)
y = np.zeros(row2 * numpoints)
sensreal = np.zeros(row2 * numpoints)
sensimag = np.zeros(row2 * numpoints)
xouter = np.zeros(numpoints+1)
youter = np.zeros(numpoints+1)
realmax = 0
imagmax = 0
idx = 0
idxfull = 0
if(ncolat%2==0):
    for idxr in range (0,row2):
        idxlong = 0
        for idxt in range (idxlow,idxlow+2):
            for idxp in range (0,nlong):
                idxfull = idxt * nlong + idxp
                idxlong = idxr * nlong + idxp
                x[idxlong] = x[idxlong] + 0.5 * dphobos4[idxr,idxfull * 5] * np.sin(dphobos4[idxr,idxfull * 5 + 1] ) * np.cos(dphobos4[idxr,idxfull * 5 + 2] )
                y[idxlong] = y[idxlong] + 0.5 * dphobos4[idxr,idxfull * 5] * np.sin(dphobos4[idxr,idxfull * 5 + 1] ) * np.sin(dphobos4[idxr,idxfull * 5 + 2] )
                sensreal[idxlong] = sensreal[idxlong] + 0.5 * dphobos4[idxr,idxfull * 5 + 3]
                sensimag[idxlong] = sensimag[idxlong] + 0.5 * dphobos4[idxr,idxfull * 5 + 4]
                realmax = max(realmax,abs(sensreal[idxlong]))
                imagmax = max(imagmax,abs(sensimag[idxlong]))
                if (idxr == idxnextlayer):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]
else:
    for idxr in range (0,row2):
        idxlong = 0
        for idxt in range (idxlow,idxlow+1):
            for idxp in range (0,nlong):
                idxfull = idxt * nlong + idxp
                idxlong = idxr * nlong + idxp
                x[idxlong] =  dphobos4[idxr,idxfull * 5] * np.sin(dphobos4[idxr,idxfull * 5 + 1] ) * np.cos(dphobos4[idxr,idxfull * 5 + 2] )
                y[idxlong] =  dphobos4[idxr,idxfull * 5] * np.sin(dphobos4[idxr,idxfull * 5 + 1] ) * np.sin(dphobos4[idxr,idxfull * 5 + 2] )
                sensreal[idxlong] =  dphobos4[idxr,idxfull * 5 + 3]
                sensimag[idxlong] =  dphobos4[idxr,idxfull * 5 + 4]
                realmax = max(realmax,abs(sensreal[idxlong]))
                imagmax = max(imagmax,abs(sensimag[idxlong]))
                if (idxr == idxnextlayer):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]

xouter[numpoints] = xouter[0]
youter[numpoints] = youter[0]
# sensreal = sensreal/realmax
# sensimag = sensimag/imagmax
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

plt.rcParams['text.usetex'] = True
# plt.tick_params(left = False, right = False , labelleft = False , 
#                 labelbottom = False, bottom = False) 
fig2, ax2 = plt.subplots(1,1)
ax2.set_aspect('equal')
tpc = ax2.tripcolor(triang, sensreal, shading='gouraud')
fig2.colorbar(tpc)
# ax2.set_title(r'Re$(\phi)$')
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.plot(xouter,youter,'k',linewidth=1)

# ax2[1].set_aspect('equal')
# tpc1 = ax2[1].tripcolor(triang, sensimag, shading='gouraud')
# fig2.colorbar(tpc1)
# ax2[1].set_title(r'Im$(\phi)$')
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