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
path3 = "Phobos/Perturbed/MatrixSolution.out"

dphobos = np.loadtxt(path0, delimiter=";")
dphobos2 = np.loadtxt(path1,delimiter=";")
dphobos3 = np.loadtxt(path2,delimiter=";")
dphobos4 = np.loadtxt(path3,delimiter=";")

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

thetamat_p = np.empty((lenln,2*l))
phimat_p = np.empty((lenln,2*l))
rmat_p = np.empty((lenln,2*l))
potrealmat_p = np.empty((lenln,2*l))
potimagmat_p = np.empty((lenln,2*l))

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

        rmat_p[idxt,idxp]=dphobos4[idxnextlayer,5*ovidx]-11100
        thetamat_p[idxt,idxp]=(3.1415926535-dphobos4[idxnextlayer,5*ovidx+1])*180.0/3.1415926535-90
        phimat_p[idxt,idxp]=dphobos4[idxnextlayer,5*ovidx+2]*180.0/3.1415926535
        potrealmat_p[idxt,idxp]=dphobos4[idxnextlayer,5*ovidx+3]
        potimagmat_p[idxt,idxp]=dphobos4[idxnextlayer,5*ovidx+4]
        ovidx = ovidx+1


# matrix to hold information on slice at equator
topval = 5
phimat3 = np.empty((topval,2*l))
rmat3 = np.empty((topval,2*l))
potrealmat3 = np.empty((topval,2*l))
potimagmat3 = np.empty((topval,2*l))
phimat4 = np.empty((topval,2*l))
rmat4 = np.empty((topval,2*l))
potrealmat4 = np.empty((topval,2*l))
potimagmat4 = np.empty((topval,2*l))
ovidx = 0
midl = l//2

for idxt in range(0,topval-1):
    idxstart = 5 * 2 * l * midl
    idx3 = idxstart
    idx4 = idxstart
    for idxp in range(0,2*l):
        rmat3[idxt,idxp]=dphobos2[idxt,idx3]
        idx3 = idx3 + 2
        phimat3[idxt,idxp]=dphobos2[idxt,idx3]
        idx3 = idx3 + 1
        potrealmat3[idxt,idxp]=dphobos2[idxt,idx3]
        idx3 = idx3 + 1
        potimagmat3[idxt,idxp]=dphobos2[idxt,idx3]
        idx3 = idx3 + 1

        rmat4[idxt,idxp]=dphobos3[idxt,idx4]
        idx4 = idx4 + 2
        phimat4[idxt,idxp]=dphobos3[idxt,idx4]
        idx4 = idx4 + 1
        potrealmat4[idxt,idxp]=dphobos3[idxt,idx4]
        idx4 = idx4 + 1
        potimagmat4[idxt,idxp]=dphobos3[idxt,idx4]
        idx4 = idx4 + 1


# difference
potdiff = abs(potrealmat-potrealmat_p)/np.amax(abs(potrealmat))*100
# plot contour map
fig, ax = plt.subplots(1,3)
pcm = ax[0].pcolormesh(phimat,thetamat,potrealmat,cmap='jet',shading='gouraud')
ax[0].set_ylabel('Latitude',fontsize=20)
ax[0].set_xlabel('Longitude',fontsize=20)
ax[0].set_title('Full solution',fontsize=20)

pcm2 = ax[1].pcolormesh(phimat_p,thetamat_p,potrealmat_p,cmap='jet',shading='gouraud')

ax[1].set_title('Perturbed solution',fontsize=20)
ax[1].set_xlabel('Longitude',fontsize=20)

pcm3 = ax[2].pcolormesh(phimat2,thetamat2,potdiff,cmap='jet',shading='gouraud')
ax[2].set_title('Relative error (%)',fontsize=20)
ax[2].set_xlabel('Longitude',fontsize=20)

# print average and max potdiff
print("Average relative error: ", np.mean(potdiff))
print("Max relative error: ", np.amax(potdiff))

p1max = np.amax(potrealmat)
p2max = np.amax(potrealmat_p)
p1min = np.amin(potrealmat)
p2min = np.amin(potrealmat_p)
pmax = max(p1max,p2max)
pmin = min(p1min,p2min)
pcm.set_clim([pmin,pmax])
pcm2.set_clim([pmin,pmax])

cbar = fig.colorbar(pcm,ax=ax[0])
cbar.ax.set_ylabel('Potential (m${}^2/$s${}^2$)',labelpad = 30, rotation=270,fontsize=20)

cbar2 = fig.colorbar(pcm2,ax=ax[1])
cbar2.ax.set_ylabel('Potential (m${}^2/$s${}^2$)',labelpad = 30, rotation=270,fontsize=20)

cbar3 = fig.colorbar(pcm3,ax=ax[2])
cbar3.ax.set_ylabel('Relative error (%)',labelpad = 30, rotation=270,fontsize=20)

plt.subplots_adjust(wspace=0.6)
plt.show()


