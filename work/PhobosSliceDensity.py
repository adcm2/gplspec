#%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
import math

import matplotlib.tri as tri
from scipy.interpolate import griddata

# import data
path0 = "Phobos/Cartesian/MatrixSolution.out"
path1 = "Phobos/Spherical/MatrixSolution.out"
path2 = "Phobos/Rotated/MatrixSolution.out"
path3 = "Phobos/Sensitivity/SensitivityKernel.out"
path4 = "Phobos/Sensitivity/SensitivityKernelSlice.out"
path5 = "Phobos/Spherical/MatrixSolutionSlice.out"
path6 = "Phobos/Spherical/MatrixSolutionRotated.out"


path6 = "Phobos/Spherical/MatrixSolutionReferentialRotated.out"
# path6 = "Phobos/Spherical/MatrixSolutionPhysicalRotated.out"
# path6 = "Phobos/Spherical/ReferentialDensitySolutionRotated.out"
# path6 = "Phobos/Spherical/PhysicalDensitySolutionRotated.out"


# "/PhysicalDensitySolutionRotated.out"
# path6 = "Phobos/Spherical/MatrixSolutionReferential.out"


dphobos = np.loadtxt(path0, delimiter=";")
# dphobos3 = np.loadtxt(path2,delimiter=";")
# dphobos2 = np.loadtxt(path3,delimiter=";")
# dphobos3 = np.loadtxt(path4,delimiter=";")
# dphobos3 = np.loadtxt(path5,delimiter=";")
dphobos4 = np.loadtxt(path6,delimiter=";")
# dphobos5 = np.loadtxt(path7,delimiter=";")

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
# colnum = 1
# maxval = 0
# for idxl in range(0,l+1):
#     for idxm in range (-idxl,idxl+1):
#         for idxn in range (0,2):
#             maxval = max(maxval,max(abs(dphobos2[:,colnum])))
#             colnum = colnum + 1


# number of longitude points
nlong = 2 * l
idxlow = 0
# bool evenn = false

ncolat = l + 1
if (ncolat%2 == 0):
    idxlow = int(ncolat/2 - 1)
else:
    idxlow = int((ncolat - 1)/2)

# print(dphobos4[:,3])
nm = idxlow * nlong + 80
print('Theta: ', dphobos4[0,1 + nm * 5])
print('Phi: ', dphobos4[0,2 + nm * 5])
print(dphobos4[:,3 + nm * 5])

# coltrial = 

# row2,col2 = dphobos3.shape
row2,col2 = dphobos4.shape
numpoints = nlong
print("numpoints: ",numpoints)
print("nlong: ",nlong)
xlen = idxnextlayer
# xlen = row2
x = np.zeros(xlen * numpoints)
y = np.zeros(xlen* numpoints)
sensreal = np.zeros(xlen* numpoints)
sensimag = np.zeros(xlen* numpoints)
xouter = np.zeros(numpoints+1)
youter = np.zeros(numpoints+1)
xouter2 = np.zeros([xlen,numpoints+1])
youter2 = np.zeros([xlen,numpoints+1])
realmax = 0
imagmax = 0
idx = 0
idxfull = 0
if(ncolat%2==0):
    for idxr in range (0,xlen):
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
                xouter2[idxr,idxp] = x[idxlong]
                youter2[idxr,idxp] = y[idxlong]
                if (idxr == xlen-1):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]
else:
    for idxr in range (0,xlen):
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
                xouter2[idxr,idxp] = x[idxlong]
                youter2[idxr,idxp] = y[idxlong]
                if (idxr == xlen-1):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]

xouter[numpoints] = xouter[0]
youter[numpoints] = youter[0]
xouter2[:,numpoints] = xouter2[:,0]
youter2[:,numpoints] = youter2[:,0]
# sensreal = sensreal/realmax
# sensimag = sensimag/imagmax
print("Number of rows: ", row2)
print("Number of cols: ", col2)
triang1 = tri.Triangulation(x,y)
triang2 = tri.Triangulation(x,y)
triang3 = tri.Triangulation(x,y)
minrad = 10.0
triang1.set_mask(np.hypot(x[triang1.triangles].mean(axis=1),
                         y[triang1.triangles].mean(axis=1))
                < minrad)
triang2.set_mask(np.hypot(x[triang2.triangles].mean(axis=1),
                         y[triang2.triangles].mean(axis=1))
                < minrad)


def apply_mask(triang, alpha=0.4):
    # Mask triangles with sidelength bigger some alpha
    triangles = triang.triangles
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    # apply masking
    triang.set_mask(maxi > alpha)



apply_mask(triang1,alpha=1450.0)
apply_mask(triang2,alpha=3000.0)





# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# tpc = ax.tripcolor(triang, sensreal, shading='flat')
# fig.colorbar(tpc)
# ax.set_title('tripcolor of Delaunay triangulation, flat shading')

# plt.rcParams['text.usetex'] = True
# plt.tick_params(left = False, right = False , labelleft = False , 
#                 labelbottom = False, bottom = False) 


fig2, ax2 = plt.subplots(1,1)
ax2.set_aspect('equal','box')
# ax2[1].set_aspect('equal')

# tpc = ax2.tripcolor(triang1, sensreal, cmap = 'jet',shading = 'gouraud' )
# ax2.tricontour(x, y, sensreal, levels=14, linewidths=0.5, colors='k')
# cntr2 = ax2.tricontourf(x, y, sensreal, levels=200, cmap="RdBu_r")
cntr2 = ax2.tricontourf(x, y, sensreal, levels=400, cmap="jet")
fig2.colorbar(cntr2, format='%.0f')
# fig2.colorbar(tpc)

############################################################
############################################################
############################################################
ax2.plot(xouter,youter,'k',linewidth=1)
for idxr in range (0,xlen - 1):
    ax2.plot(xouter2[idxr,:],youter2[idxr,:],'k--',linewidth=1)
ax2.plot(xouter2[xlen-1,:],youter2[xlen-1,:],'k',linewidth=1)
ax2.fill_between(xouter,youter,1.2*youter,color='white')
ax2.set(xlim=(-12000,12000),ylim=(-12000,12000))
############################################################
############################################################
############################################################2

# X, Y= np.meshgrid(x,y)

# print('Minimum x', x.min())
# print('Maximum x', x.max())
# print('Minimum y', y.min())
# print('Maximum y', y.max())

# xpts = np.arange(x.min(),x.max(),10.0)
# ypts = np.arange(y.min(),y.max(),10.0)
# xplot,yplot = np.meshgrid(xpts,ypts)

# znew = griddata((x,y),sensreal,(xplot,yplot),method='linear')

# ax2.pcolormesh(xplot,yplot,znew)

# # ax2.set_title(r'Re$(\phi)$')
# ax2.set_yticklabels([])
# ax2.set_xticklabels([])


# ax2[1].plot(xouter,youter,'k',linewidth=1)
# for idxr in range (0,xlen):
#     ax2[1].plot(xouter2[idxr,:],youter2[idxr,:],'k',linewidth=1)

# ax2[1].set_aspect('equal')
# tpc1 = ax2[1].tripcolor(triang1, sensimag, shading='gouraud')
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