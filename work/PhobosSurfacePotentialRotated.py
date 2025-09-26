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
lenln = l+1  

# find index of outer radius
# this assumes that the length norm is the referential radius of the planet
idxnextlayer = 0
for idx in range(0,row-1):
    if (dphobos[idx+1,0] - dphobos[idx,0] == 1):
        idxnextlayer = idx +1


# get theta, phi, r
thetamat = np.empty((lenln,2*l))
phimat = np.empty((lenln,2*l))
rmat = np.empty((lenln,2*l))
potrealmat = np.empty((lenln,2*l))
potimagmat = np.empty((lenln,2*l))

thetamat2 = np.empty((lenln,2*l))
phimat2 = np.empty((lenln,2*l))
rmat2 = np.empty((lenln,2*l))
potrealmat2 = np.empty((lenln,2*l))
potimagmat2 = np.empty((lenln,2*l))
ovidx = 0
for idxt in range(0,lenln):
    for idxp in range(0,2*l):
        rmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx]-11100
        thetamat[idxt,idxp]=(3.1415926535-dphobos2[idxnextlayer,5*ovidx+1])*180.0/3.1415926535-90
        phimat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+2]*180.0/3.1415926535
        potrealmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+3]
        potimagmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+4]

        rmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx]-11100
        thetamat2[idxt,idxp]=(3.1415926535-dphobos3[idxnextlayer,5*ovidx+1])*180.0/3.1415926535-90
        phimat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+2]*180.0/3.1415926535
        potrealmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+3]
        potimagmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+4]
        ovidx = ovidx+1


fig, ax = plt.subplots(1,2)
pcm = ax[0].pcolormesh(phimat,thetamat,potrealmat,cmap='jet',shading='gouraud')
ax[0].set_ylabel('Latitude',fontsize=20,)
ax[0].set_xlabel('Longitude',fontsize=20)
ax[0].tick_params(axis='both', which='major', labelsize=15)

pcm2 = ax[1].pcolormesh(phimat2,thetamat2,potrealmat2,cmap='jet',shading='gouraud')
ax[1].set_xlabel('Longitude',fontsize=20)
ax[1].tick_params(axis='both', which='major', labelsize=15)
cbar2 = fig.colorbar(pcm2)
cbar2.ax.set_ylabel('Potential (m${}^2/$s${}^2$)',labelpad = 30, rotation=270,fontsize=20)
cbar2.ax.tick_params(axis='both', which='major', labelsize=15)


plt.show()


