#%%
import numpy as np
import matplotlib.pyplot as plt

# d = np.loadtxt("Lagrange.out", delimiter=";")
d = np.loadtxt("gravitycheck.out")
dp = np.loadtxt("ellipticity2.out")
# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(1,2)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()



SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 30
MLWIDTH = 2.0

# density
# ax[0].plot((d[-1,0]-d[:, 0])/1000, d[:,3], "b", linewidth = MLWIDTH)
ax[0].plot(d[:,0], d[:,3], "b", linewidth = MLWIDTH)
ax[0].plot(dp[:,0],dp[:,4],"r",linewidth = MLWIDTH)
ax[0].plot(dp[:,0],dp[:,6],"k",linewidth = MLWIDTH)
# ax[0].plot((d[-1,0]-d[:, 0])/1000, d[:,4], "r", linewidth = MLWIDTH)
# ax[0].plot(d[:, 0], d[:,2], "r--", linewidth = MLWIDTH)
ax[0].legend(['FD', 'Clairaut','Potential FE'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)

ax[0].set_xlabel("Depth ($km$)", fontsize = MEDIUM_SIZE)
ax[0].set_ylabel("Ellipticity", fontsize = MEDIUM_SIZE)


ax[1].plot(dp[:,0], 100*abs(np.interp(dp[:,0],d[:,0],d[:,3])-dp[:,4])/max(dp[:,4]), "b", linewidth = MLWIDTH)
ax[1].plot(dp[:,0], 100*abs(dp[:,5]-dp[:,4])/max(dp[:,4]), "r", linewidth = MLWIDTH)
ax[1].plot(dp[:,0], 100*abs(dp[:,6]-dp[:,4])/max(dp[:,4]), "k", linewidth = MLWIDTH)
# ax[1].plot(dp[:,0],dp[:,6],"k",linewidth = MLWIDTH)
# ax[1].plot(dp[:,0],dp[:,5],"k",linewidth = MLWIDTH)
# ax[1].plot(dp[:,0],abs(dp[:,4] -dp[:,5]),"k",linewidth = MLWIDTH)
# ax[1].plot(dp[:,0],dp[:,3]-dp[:,4],"r",linewidth = MLWIDTH)
# ax[1].plot(dp[:,0],dp[:,5],"r")

# ax[1].plot(dp[:,0],100 * (dp[:,2]-dp[:,3])/max(dp[:,2]),"k")
# ax[1].plot(d[:,0],(d[:,1]-d[:,2])/max(abs(d[:,1])) * 100,"b",linewidth=MLWIDTH)
ax[1].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
# ax[1].set_ylabel("", fontsize = BIGGER_SIZE)
ax[1].set_xlabel("Radius", fontsize = BIGGER_SIZE)
ax[1].legend(['FD Rel. Err. (%)', 'CLFE2 (%)', 'POTFE (%)'], fontsize = MEDIUM_SIZE)


plt.show()

