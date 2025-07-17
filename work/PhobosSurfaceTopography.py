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

dphobos = np.loadtxt(path0, delimiter=";")
dphobos2 = np.loadtxt(path1,delimiter=";")
dphobos3 = np.loadtxt(path2,delimiter=";")

# find l, etc
row, col = dphobos.shape
numl = int((col-1)/5)
l = int((math.sqrt(1 + 2*numl)-1)/2)
print("Number of rows: ", row)
print("Number of cols: ", col)
print("l: ", l)

idxnextlayer = 0
for idx in range(0,row-1):
    if (dphobos[idx+1,0] - dphobos[idx,0] == 1):
        idxnextlayer = idx +1
lenln = l+1  
print(idxnextlayer)

# get theta, phi, r
thetamat = np.empty((lenln,2*l))
phimat = np.empty((lenln,2*l))
rmat = np.empty((lenln,2*l))
potrealmat = np.empty((lenln,2*l))

thetamat2 = np.empty((lenln,2*l))
phimat2 = np.empty((lenln,2*l))
rmat2 = np.empty((lenln,2*l))
potrealmat2 = np.empty((lenln,2*l))
potimagmat2 = np.empty((lenln,2*l))
ovidx = 0
for idxt in range(0,lenln):
    for idxp in range(0,2*l):
        # thetamat[idxt,idxp]=90 - 180/3.1415926535*dphobos2[idxnextlayer,5*ovidx+1]
        # phimat[idxt,idxp]=180 + 180/3.1415926535* dphobos2[idxnextlayer,5*ovidx+2]
        thetamat[idxt,idxp]=(3.1415926535-dphobos2[idxnextlayer,5*ovidx+1])*180.0/3.1415926535-90
        phimat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+2]*180.0/3.1415926535
        rmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx]-11100
        potrealmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+3]

        rmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx]-11100
        thetamat2[idxt,idxp]=90 - 180/3.1415926535*dphobos3[idxnextlayer,5*ovidx+1]
        # thetamat2[idxt,idxp]=3.1415926535-dphobos3[idxnextlayer,5*ovidx+1]
        phimat2[idxt,idxp]=180 + 180/3.1415926535* dphobos3[idxnextlayer,5*ovidx+2]
        # phimat2[idxt,idxp]=d/phobos3[idxnextlayer,5*ovidx+2]
        # potrealmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+3]
        # potimagmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+4]

        ovidx = ovidx+1


# plot contour map
# fig, ax = plt.subplots(1,2)
# pcm = ax[0].pcolormesh(phimat,thetamat,potrealmat,cmap='hsv_r',shading='gouraud')
# CS=ax[0].contour(phimat,thetamat,potrealmat,colors='k')
# ax[0].clabel(CS,inline=1,fontsize=10)
# # pcm.set_clim([-5000,3000])
# fig.colorbar(pcm,ax=ax[0])

# pcm1 = ax[1].pcolormesh(phimat2,thetamat2,rmat2,cmap='hsv_r',shading='gouraud')
# CS1=ax[1].contour(phimat2,thetamat2,rmat2,[-2000,-1500,-1000,-500,0,500,1000,1500,2000],colors='k')
# ax[1].clabel(CS1,inline=1,fontsize=10)
# # pcm1.set_clim([-5000,3000])
# fig.colorbar(pcm1,ax=ax[1])
# plt.show()

# plot contour map
fig, ax = plt.subplots(1,1)
pcm = ax.pcolormesh(phimat,thetamat,rmat,cmap='jet',shading='gouraud')
CS=ax.contour(phimat,thetamat,rmat,[-2000,-1000,0,1000,2000],colors='k')
ax.clabel(CS,inline=1,fontsize=15)
# pcm.set_clim([-5000,3000])
# fig.colorbar(pcm,ax)
# cbar = fig.colorbar(pcm)
# cbar.ax.set_ylabel('Potential ($m^2/s^2$)',labelpad = 30, rotation=270,fontsize=20)
ax.set_ylabel('Latitude',fontsize=20)
ax.set_xlabel('Longitude',fontsize=20)
cbar = fig.colorbar(pcm)
cbar.ax.set_ylabel('Elevation ($m$)',labelpad = 30, rotation=270,fontsize=20)

# pcm1 = ax[1].pcolormesh(phimat2,thetamat2,rmat2,cmap='hsv_r',shading='gouraud')
# CS1=ax[1].contour(phimat2,thetamat2,rmat2,[-2000,-1500,-1000,-500,0,500,1000,1500,2000],colors='k')
# ax[1].clabel(CS1,inline=1,fontsize=10)
# # pcm1.set_clim([-5000,3000])
# fig.colorbar(pcm1,ax=ax[1])
plt.show()
