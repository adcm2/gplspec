import numpy as np
import matplotlib.pyplot as plt
import math

# import data
path0 = "Phobos/Cartesian/MatrixSolution.out"

path_rp = "Phobos/Spherical/MatrixSolutionReferentialRotated.out"
path_pp = "Phobos/Spherical/MatrixSolutionPhysicalRotated.out"
path_rd = "Phobos/Spherical/ReferentialDensitySolutionRotated.out"
path_pd = "Phobos/Spherical/PhysicalDensitySolutionRotated.out"

dphobos = np.loadtxt(path0, delimiter=";")
dphobs_pd = np.loadtxt(path_pd,delimiter=";")
dphobs_rd = np.loadtxt(path_rd,delimiter=";")
dphobs_pp = np.loadtxt(path_pp,delimiter=";")
dphobs_rp = np.loadtxt(path_rp,delimiter=";")

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
        idxnextlayer = idx +2
lenln = l+1  

# number of longitude points
nlong = 2 * l
idxlow = 0

ncolat = l + 1
if (ncolat%2 == 0):
    idxlow = int(ncolat/2 - 1)
else:
    idxlow = int((ncolat - 1)/2)

nm = idxlow * nlong + 80
print('Theta: ', dphobs_pd[0,1 + nm * 5])
print('Phi: ', dphobs_pd[0,2 + nm * 5])
print(dphobs_pd[:,3 + nm * 5])

row2,col2 = dphobs_pd.shape
numpoints = nlong
print("numpoints: ",numpoints)
print("nlong: ",nlong)
xlen = idxnextlayer

# physical coordinates
x = np.zeros(xlen * numpoints)
y = np.zeros(xlen* numpoints)

xouter = np.zeros(numpoints+1)
youter = np.zeros(numpoints+1)
xouter2 = np.zeros([xlen,numpoints+1])
youter2 = np.zeros([xlen,numpoints+1])

# referential coordinates
xr = np.zeros(xlen * numpoints)
yr = np.zeros(xlen* numpoints)

xouterr = np.zeros(numpoints+1)
youterr = np.zeros(numpoints+1)
xouter2r = np.zeros([xlen,numpoints+1])
youter2r = np.zeros([xlen,numpoints+1])

# physical density
pdr = np.zeros(xlen* numpoints)
pdi = np.zeros(xlen* numpoints)

# referential density
rdr = np.zeros(xlen* numpoints)
rdi = np.zeros(xlen* numpoints)

# physical potential
ppr = np.zeros(xlen* numpoints)
ppi = np.zeros(xlen* numpoints)

# referential potential
rpr = np.zeros(xlen* numpoints)
rpi = np.zeros(xlen* numpoints)

idx = 0
idxfull = 0
if(ncolat%2==0):
    for idxr in range (0,xlen):
        idxlong = 0
        for idxt in range (idxlow,idxlow+2):
            for idxp in range (0,nlong):
                # indices
                idxfull = idxt * nlong + idxp
                idxlong = idxr * nlong + idxp
                
                # coordinates
                x[idxlong] = x[idxlong] + 0.5 * dphobs_pd[idxr,idxfull * 5] * np.sin(dphobs_pd[idxr,idxfull * 5 + 1] ) * np.cos(dphobs_pd[idxr,idxfull * 5 + 2] )
                y[idxlong] = y[idxlong] + 0.5 * dphobs_pd[idxr,idxfull * 5] * np.sin(dphobs_pd[idxr,idxfull * 5 + 1] ) * np.sin(dphobs_pd[idxr,idxfull * 5 + 2] )
                
                # coordinates
                xr[idxlong] = xr[idxlong] + 0.5 * dphobs_rd[idxr,idxfull * 5] * np.sin(dphobs_rd[idxr,idxfull * 5 + 1] ) * np.cos(dphobs_rd[idxr,idxfull * 5 + 2] )
                yr[idxlong] = yr[idxlong] + 0.5 * dphobs_rd[idxr,idxfull * 5] * np.sin(dphobs_rd[idxr,idxfull * 5 + 1] ) * np.sin(dphobs_rd[idxr,idxfull * 5 + 2] )

                # physical density
                pdr[idxlong] = pdr[idxlong] + 0.5 * dphobs_pd[idxr,idxfull * 5 + 3]
                pdi[idxlong] = pdi[idxlong] + 0.5 * dphobs_pd[idxr,idxfull * 5 + 4]
                
                # referential density
                rdr[idxlong] = rdr[idxlong] + 0.5 * dphobs_rd[idxr,idxfull * 5 + 3]
                rdi[idxlong] = rdi[idxlong] + 0.5 * dphobs_rd[idxr,idxfull * 5 + 4]
                
                # physical potential
                ppr[idxlong] = ppr[idxlong] + 0.5 * dphobs_pp[idxr,idxfull * 5 + 3]
                ppi[idxlong] = ppi[idxlong] + 0.5 * dphobs_pp[idxr,idxfull * 5 + 4]
                
                # referential potential
                rpr[idxlong] = rpr[idxlong] + 0.5 * dphobs_rp[idxr,idxfull * 5 + 3]
                rpi[idxlong] = rpi[idxlong] + 0.5 * dphobs_rp[idxr,idxfull * 5 + 4]

                xouter2[idxr,idxp] = x[idxlong]
                youter2[idxr,idxp] = y[idxlong]
                if (idxr == xlen-1):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]
                    
                xouter2r[idxr,idxp] = xr[idxlong]
                youter2r[idxr,idxp] = yr[idxlong]
                if (idxr == xlen-1):
                    xouterr[idxp] = xr[idxlong]
                    youterr[idxp] = yr[idxlong]
else:
    for idxr in range (0,xlen):
        idxlong = 0
        for idxt in range (idxlow,idxlow+1):
            for idxp in range (0,nlong):
                # indices
                idxfull = idxt * nlong + idxp
                idxlong = idxr * nlong + idxp
                
                # coordinates
                x[idxlong] =  dphobs_pd[idxr,idxfull * 5] * np.sin(dphobs_pd[idxr,idxfull * 5 + 1] ) * np.cos(dphobs_pd[idxr,idxfull * 5 + 2] )
                y[idxlong] =  dphobs_pd[idxr,idxfull * 5] * np.sin(dphobs_pd[idxr,idxfull * 5 + 1] ) * np.sin(dphobs_pd[idxr,idxfull * 5 + 2] )
                
                # coordinates
                xr[idxlong] =  dphobs_rd[idxr,idxfull * 5] * np.sin(dphobs_rd[idxr,idxfull * 5 + 1] ) * np.cos(dphobs_rd[idxr,idxfull * 5 + 2] )
                yr[idxlong] =  dphobs_rd[idxr,idxfull * 5] * np.sin(dphobs_rd[idxr,idxfull * 5 + 1] ) * np.sin(dphobs_rd[idxr,idxfull * 5 + 2] )
                
                # physical density
                pdr[idxlong] =  dphobs_pd[idxr,idxfull * 5 + 3]
                pdi[idxlong] =  dphobs_pd[idxr,idxfull * 5 + 4]
                
                # referential density
                rdr[idxlong] =  dphobs_rd[idxr,idxfull * 5 + 3]
                rdi[idxlong] =  dphobs_rd[idxr,idxfull * 5 + 4]
                
                # physical potential
                ppr[idxlong] =  dphobs_pp[idxr,idxfull * 5 + 3]
                ppi[idxlong] =  dphobs_pp[idxr,idxfull * 5 + 4]
                
                # referential potential
                rpr[idxlong] =  dphobs_rp[idxr,idxfull * 5 + 3]
                rpi[idxlong] =  dphobs_rp[idxr,idxfull * 5 + 4]

                xouter2[idxr,idxp] = x[idxlong]
                youter2[idxr,idxp] = y[idxlong]
                if (idxr == xlen-1):
                    xouter[idxp] = x[idxlong]
                    youter[idxp] = y[idxlong]
                
                xouter2r[idxr,idxp] = xr[idxlong]
                youter2r[idxr,idxp] = yr[idxlong]
                if (idxr == xlen-1):
                    xouterr[idxp] = xr[idxlong]
                    youterr[idxp] = yr[idxlong]

xouter[numpoints] = xouter[0]
youter[numpoints] = youter[0]
xouter2[:,numpoints] = xouter2[:,0]
youter2[:,numpoints] = youter2[:,0]

xouterr[numpoints] = xouterr[0]
youterr[numpoints] = youterr[0]
xouter2r[:,numpoints] = xouter2r[:,0]
youter2r[:,numpoints] = youter2r[:,0]

fig2, ax2 = plt.subplots(2,2)
ax2[0,0].set_aspect('equal','box')
cntr2 = ax2[0,0].tricontourf(x, y, pdr, levels=400, cmap="jet")
cb00 = fig2.colorbar(cntr2, format='%.1f')
cb00.set_label('Density $(kg/m^3)$', rotation=270, labelpad=30, fontsize=20)
cb00.ax.tick_params(axis='both', which='major', labelsize=15)
# ax2[0,0].set_xlabel('x $(m)$', fontsize=20)
ax2[0,0].set_ylabel('y $(m)$', fontsize=20)
ax2[0,0].tick_params(axis='both', which='major', labelsize=15)

ax2[0,1].set_aspect('equal','box')
cntr2 = ax2[0,1].tricontourf(xr, yr, rdr, levels=400, cmap="jet")
cb01 = fig2.colorbar(cntr2, format='%.0f')
cb01.set_label('Density $(kg/m^3)$', rotation=270, labelpad=30, fontsize=20)
cb01.ax.tick_params(axis='both', which='major', labelsize=15)
# ax2[0,1].set_xlabel('x $(m)$', fontsize=20)
# ax2[0,1].set_ylabel('y $(m)$', fontsize=20)
ax2[0,1].tick_params(axis='both', which='major', labelsize=15)

ax2[1,0].set_aspect('equal','box')
cntr2 = ax2[1,0].tricontourf(x, y, ppr, levels=400, cmap="jet")
cb10 = fig2.colorbar(cntr2, format='%.0f')
cb10.set_label('Potential $(m^2/s^2)$', rotation=270, labelpad=30, fontsize=20)
cb10.ax.tick_params(axis='both', which='major', labelsize=15)
ax2[1,0].set_xlabel('x $(m)$', fontsize=20)
ax2[1,0].set_ylabel('y $(m)$', fontsize=20)
ax2[1,0].tick_params(axis='both', which='major', labelsize=15)

ax2[1,1].set_aspect('equal','box')
cntr2 = ax2[1,1].tricontourf(xr, yr, rpr, levels=400, cmap="jet")
cb11 = fig2.colorbar(cntr2, format='%.0f')
cb11.set_label('Potential $(m^2/s^2)$', rotation=270, labelpad=30, fontsize=20)
cb11.ax.tick_params(axis='both', which='major', labelsize=15)
ax2[1,1].set_xlabel('x $(m)$', fontsize=20)
# ax2[1,1].set_ylabel('y $(m)$', fontsize=20)
ax2[1,1].tick_params(axis='both', which='major', labelsize=15)

ax2[0,0].plot(xouter,youter,'k',linewidth=1)
for idxr in range (0,xlen - 1):
    ax2[0,0].plot(xouter2[idxr,:],youter2[idxr,:],'k--',linewidth=1)
ax2[0,0].plot(xouter2[xlen-1,:],youter2[xlen-1,:],'k',linewidth=1)
ax2[0,0].fill_between(xouter2[xlen-1,:],youter2[xlen-1,:],1.5*youter2[xlen-1,:],color='white')
ax2[0,0].set(xlim=(-13000,13000),ylim=(-13000,13000))
ax2[0,0].text(-10000,10000,'A', fontsize=15,fontweight='bold')

ax2[0,1].plot(xouterr,youterr,'k',linewidth=1)
for idxr in range (0,xlen - 1):
    ax2[0,1].plot(xouter2r[idxr,:],youter2r[idxr,:],'k--',linewidth=1)
ax2[0,1].plot(xouter2r[xlen-1,:],youter2r[xlen-1,:],'k',linewidth=1)
ax2[0,1].fill_between(xouterr,youterr,1.2*youterr,color='white')
ax2[0,1].set(xlim=(-13000,13000),ylim=(-13000,13000))
ax2[0,1].text(-10000,10000,'B', fontsize=15,fontweight='bold')

ax2[1,0].plot(xouter,youter,'k',linewidth=1)
for idxr in range (0,xlen - 1):
    ax2[1,0].plot(xouter2[idxr,:],youter2[idxr,:],'k--',linewidth=1)
ax2[1,0].plot(xouter2[xlen-1,:],youter2[xlen-1,:],'k',linewidth=1)
ax2[1,0].fill_between(xouter2[xlen-1,:],youter2[xlen-1,:],1.5*youter2[xlen-1,:],color='white')
ax2[1,0].set(xlim=(-13000,13000),ylim=(-13000,13000))
ax2[1,0].text(-10000,10000,'C', fontsize=15,fontweight='bold')

ax2[1,1].plot(xouterr,youterr,'k',linewidth=1)
for idxr in range (0,xlen - 1):
    ax2[1,1].plot(xouter2r[idxr,:],youter2r[idxr,:],'k--',linewidth=1)
ax2[1,1].plot(xouter2r[xlen-1,:],youter2r[xlen-1,:],'k',linewidth=1)
ax2[1,1].fill_between(xouterr,youterr,1.2*youterr,color='white')
ax2[1,1].set(xlim=(-13000,13000),ylim=(-13000,13000))
ax2[1,1].text(-10000,10000,'D', fontsize=15,fontweight='bold')

plt.show()
plt.subplots_adjust(left=0.034, bottom=0.045, right=0.9, top=0.974, wspace=0.2, hspace=0.128)


