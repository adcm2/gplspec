 #%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
import math
# import plotly.graph_objects as go

# d = np.loadtxt("Lagrange.out", delimiter=";")
path_to_folder = "Bench10"
path1 = path_to_folder + "/exact/MatrixSolution.out"
# path2 = path_to_folder + "/IntegralSolution.out"
path2 = path_to_folder + "/perturbed/MatrixSolution.out"
path3 = path_to_folder + "/unperturbed/MatrixSolution.out"
path4 = path_to_folder + "/perturbed2/MatrixSolution.out"
path5 = "Phobos/Cartesian/MatrixSolution.out"
path6 = "Phobos/Spherical/MatrixSolution.out"
# d = np.loadtxt("cleanbench1.out", delimiter=";")
# d0 = np.loadtxt(path1, delimiter=";")
# d1 = np.loadtxt(path2, delimiter=";")
dexact = np.loadtxt(path1, delimiter=";")
dp = np.loadtxt(path2, delimiter=";")
dup = np.loadtxt(path3, delimiter=";")
dp2 = np.loadtxt(path4, delimiter=";")
dphobos = np.loadtxt(path5, delimiter=";")
dphobos2 = np.loadtxt(path6,delimiter=";")
# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
# fig, ax = plt.subplots(2,2)
# plt.rcParams.update({'font.size': 18})
# fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

n = 5
m1 = 200
m2 = 10
c1 = n*m1
c11 = n*m1 + 3
c2 = n*m2
c21 = n*m2+3

row, col = dphobos.shape
numl = int((col-1)/5)
l = int((math.sqrt(1 + 2*numl)-1)/2)
print("Number of rows: ", row)
print("Number of cols: ", col)
print("l: ", l)

# get the values for one of the referential radii
x = np.empty(numl)
y = np.empty(numl)
z = np.empty(numl)

for idx in range(0,numl):
    x[idx] = dphobos[12,5 * idx]
    y[idx] = dphobos[12,5 * idx+1]
    z[idx] = dphobos[12,5 * idx+2]



idxnextlayer = 0
for idx in range(0,row-1):
    if (dphobos[idx+1,0] - dphobos[idx,0] == 1):
        idxnextlayer = idx +1

# idxnextlayer=9
# lenln = int(l/2)
lenln = l+1  

xmat = np.empty((lenln,2*l))
ymat = np.empty((lenln,2*l))
zmat = np.empty((lenln,2*l))
ovidx = 0
for idxt in range(0,lenln):
    for idxp in range(0,2*l):
        xmat[idxt,idxp]=dphobos[idxnextlayer,5*ovidx+1]
        ymat[idxt,idxp]=dphobos[idxnextlayer,5*ovidx+2]
        zmat[idxt,idxp]=dphobos[idxnextlayer,5*ovidx+3]
        # if(xmat[idxt,idxp]<0):
        #     zmat[idxt,idxp]=dphobos[15,5*ovidx+2]
        # zmat[idxt,idxp] = 1
        ovidx = ovidx+1

ovidx = 0
thetamat = np.empty((lenln,2*l))
phimat = np.empty((lenln,2*l))
rmat = np.empty((lenln,2*l))
ovidx = 0
for idxt in range(0,lenln):
    for idxp in range(0,2*l):
        thetamat[idxt,idxp]=90 - 180/3.1415926535*dphobos2[idxnextlayer,5*ovidx+1]
        phimat[idxt,idxp]=180 + 180/3.1415926535* dphobos2[idxnextlayer,5*ovidx+2]
        # if (phimat[idxt,idxp] > 360.0):
        #     phimat[idxt,idxp] = phimat[idxt,idxp] - 360
        rmat[idxt,idxp]=dphobos2[idxnextlayer,5*ovidx]-11100
        # if(xmat[idxt,idxp]<0):
        #     zmat[idxt,idxp]=dphobos[15,5*ovidx+2]
        # zmat[idxt,idxp] = 1
        ovidx = ovidx+1

# print(thetamat[0:3,0:3])
# print(phimat[0:3,0:3])
# print(rmat[0:3,0:3])
# potential 
# ax[0,0].plot(dexact[:,0],dexact[:,3], "k", linewidth=3)
# ax[0,0].plot(dp[:, 0], dp[:, 3], "b", linewidth=3)
# ax[0,0].plot(dup[:, 0], dup[:, 3], "r--", linewidth=3)
# ax[0,0].plot(dp2[:, 0], dp2[:, 3], "g-.", linewidth=3)
# ax[0,0].legend([ 'Exact','Perturbed', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
# ax[0,0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
# t = ax[0,0].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[0,0].set_ylabel("Potential", fontsize = BIGGER_SIZE)


# bounding box
blue = np.array([0., 0., 2.])
rgb = np.tile(blue, (zmat.shape[0], zmat.shape[1], 1))
ls = LightSource()
illuminated_surface=ls.shade_rgb(rgb,zmat)

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111,projection='3d')
# # ax2.plot_surface(xmat,ymat,zmat,linewidth=0,antialiased=False,facecolors=illuminated_surface)
# ax2.plot_surface(xmat,ymat,zmat,linewidth=0,antialiased=False,cmap='viridis')
# # ax2.plot_trisurf(x,y,z,cmap='viridis')
# ax2.set_box_aspect([1,1,1])
# # set_axes_equal(ax2)
# ax2.set_xlim3d([-12000,12000])
# ax2.set_ylim3d([-12000,12000])
# ax2.set_zlim3d([-12000,12000])
# # ax2.antia
# ax2.set_xlabel('X axis')
# ax2.set_ylabel('Y axis')
# ax2.set_zlabel('Z axis')

# ax3 = fig2.add_subplot(111,projection='3d')

# ax3.imshow(rmat,interpolation='bilinear')

print(rmat.max())
print(rmat.min())

fig, ax = plt.subplots()
pcm = ax.pcolormesh(phimat,thetamat,rmat,cmap='hsv_r')
CS=ax.contour(phimat,thetamat,rmat,[-2000,-1500,-1000,-500,0,500,1000,1500,2000],colors='k')
ax.clabel(CS,inline=1,fontsize=10)
# ax.invert_xaxis()
pcm.set_clim([-5000,3000])
fig.colorbar(pcm,ax=ax)
plt.show()


# ax2.plot_surface(dphobos[:,0],dphobos[:,1],dphobos[:,2],cmap='viridis')


# ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp[:,3])/abs(max(dexact[:,3])) * 100.0, "b", linewidth=3)
# ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dup[:,3])/abs(max(dexact[:,3])) * 100.0, "r", linewidth=3)
# ax[1,0].plot(dexact[:,0],abs(dexact[:,3] - dp2[:,3])/abs(max(dexact[:,3])) * 100.0, "g", linewidth=3)
# ax[1,0].legend(['Perturbed solution', 'Unperturbed solution', 'Perturb a'],fontsize = MEDIUM_SIZE)
# ax[1,0].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
# t = ax[1,0].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[1,0].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
# ax[1,0].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)


# ax[0,1].plot(dexact[:,c1],dexact[:,c11], "k", linewidth=3)
# ax[0,1].plot(dp[:, c1], dp[:, c11], "b", linewidth=3)
# ax[0,1].plot(dup[:, c1], dup[:, c11], "r--", linewidth=3)
# ax[0,1].plot(dp2[:, c1], dp2[:, c11], "g-.", linewidth=3)
# ax[0,1].legend([ 'Exact','Perturbed', 'Unperturbed','Perturb a'], fontsize = MEDIUM_SIZE)
# ax[0,1].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
# t = ax[0,1].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[0,1].set_ylabel("Potential", fontsize = BIGGER_SIZE)

# ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp[:,c11])/abs(max(dexact[:,c11])) * 100.0, "b", linewidth=3)
# ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dup[:,c11])/abs(max(dexact[:,c11])) * 100.0, "r", linewidth=3)
# ax[1,1].plot(dexact[:,c1],abs(dexact[:,c11] - dp2[:,c11])/abs(max(dexact[:,c11])) * 100.0, "g", linewidth=3)
# ax[1,1].legend(['Perturbed solution', 'Unperturbed solution','Perturb a'],fontsize = MEDIUM_SIZE)
# ax[1,1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
# t = ax[1,1].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[1,1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
# ax[1,1].set_ylabel("Rel diff (%)", fontsize = BIGGER_SIZE)
# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
#figure = plt.gcf()  
#figure.set_size_inches(9, 12) 
#plt.savefig("../../../Notes/New_coupling/figures/Benchmark_1.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
# plt.show()
# 
