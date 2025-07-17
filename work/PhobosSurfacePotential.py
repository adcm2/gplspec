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

# find index of outer radius
# this assumes that the length norm is the referential radius of the planet
idxnextlayer = 0
for idx in range(0,row-1):
    # print("idx: ")
    # print(idx)
    # print(f"{dphobos[idx,0]:.6}")
    if (dphobos[idx+1,0] - dphobos[idx,0] == 1):
        idxnextlayer = idx +1
lenln = l+1  
# idxnextlayer = row-1
# print(idxnextlayer)
# print(row)


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
        # thetamat[idxt,idxp]=90 - 180/3.1415926535*dphobos2[idxnextlayer,5*ovidx+1]
        thetamat[idxt,idxp]=(3.1415926535-dphobos2[idxnextlayer,5*ovidx+1])*180.0/3.1415926535-90
        # phimat[idxt,idxp]=180 + 180/3.1415926535* dphobos2[idxnextlayer,5*ovidx+2]
        phimat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+2]*180.0/3.1415926535
        potrealmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+3]
        potimagmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx+4]

        rmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx]-11100
        # thetamat2[idxt,idxp]=90 - 180/3.1415926535*dphobos3[idxnextlayer,5*ovidx+1]
        thetamat2[idxt,idxp]=3.1415926535-dphobos3[idxnextlayer,5*ovidx+1]
        # phimat2[idxt,idxp]=180 + 180/3.1415926535* dphobos3[idxnextlayer,5*ovidx+2]
        phimat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+2]
        potrealmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+3]
        potimagmat2[idxt,idxp]=dphobos3[idxnextlayer,5*ovidx+4]
        ovidx = ovidx+1


# matrix to hold information on slice at equator
# thetamat = np.empty((row,2*l))
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
        # print(idxt, " ", idx3)
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

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Helvetica"
# })

fig, ax = plt.subplots(1,1)
pcm = ax.pcolormesh(phimat,thetamat,potrealmat,cmap='jet',shading='gouraud')
ax.set_ylabel('Latitude',fontsize=20)
ax.set_xlabel('Longitude',fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = fig.colorbar(pcm)
cbar.ax.set_ylabel('Potential (m${}^2/$s${}^2$)',labelpad = 30, rotation=270,fontsize=20)
cbar.ax.tick_params(axis='both', which='major', labelsize=15)
# pcm2 = ax[1].pcolormesh(phimat2,thetamat2,potrealmat2,cmap='hsv_r',shading='gouraud')
# # fig = plt.figure(),weight='bold'
# # fig = plt.figure()

#     # rho, phi = cart2pol(x,y)
#     # plt.subplot(projection="polar")
# # ax[1,0].projection("polar")
# ax[1,0].remove()
# ax3 = fig.add_subplot(2,2,3,projection="polar")
# # ax3.scatter(phimat4,rmat4)
# ax3.pcolormesh(phimat4,rmat4,potrealmat4)
#     # plt.pcolormesh(phi, rho, z)
#     # #plt.pcolormesh(th, z, r)

#     # plt.plot(phi, rho, color='k', ls='none')
#     # plt.grid()

#     # plt.show()
# # CS=ax.contour(phimat,thetamat,potrealmat,[-2000,-1500,-1000,-500,0,500,1000,1500,2000],colors='k')
# # ax.clabel(CS,inline=1,fontsize=10)
# # pcm.set_clim([-5000,3000])
# fig.colorbar(pcm,ax=ax[0,0])
# fig.colorbar(pcm2,ax=ax[0,1])

# ax = plt.subplot(projection="polar")
# ax.pcolormesh(phimat4,rmat4,potrealmat4)
# ax.scatter(phimat4,rmat4)
# ax.plot(phimat4,rmat4)
# ax.contourf(phimat4,rmat4,potrealmat4)


plt.show()


