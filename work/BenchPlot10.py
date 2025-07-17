#%%
import numpy as np
import matplotlib.pyplot as plt

# d = np.loadtxt("Lagrange.out", delimiter=";")
d = np.loadtxt("bench9.out", delimiter=";")
# d2 = np.loadtxt('time_test.STAT1.Z')

# initialise plot
fig, ax = plt.subplots(2,1)
# plt.rcParams.update({'font.size': 18})
fig.tight_layout()

SMALL_SIZE = 13
MEDIUM_SIZE = 20
BIGGER_SIZE = 20


# potential 
ax[0].plot(d[:,0],d[:,1], "b", linewidth=3)
ax[0].plot(d[:, 0], d[:, 2], "r--", linewidth=3)
ax[0].plot(d[:, 0], d[:, 6], "g-.", linewidth=3)
ax[0].plot(d[:, 0], d[:, 1], "ko")
ax[0].plot(d[:, 0], d[:, 4], "k-.", linewidth=3)
ax[0].legend(['Integral','Matrix-free pseudospectral', 'BP Theory','Nodes', 'Spherical'], fontsize = MEDIUM_SIZE)
ax[0].tick_params(axis='both', which='major', labelsize=SMALL_SIZE)
t = ax[0].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[0].set_ylabel("Potential", fontsize = BIGGER_SIZE)

# error in potential
# ax[1].plot(d[:,0], abs(d[:,2] - d[:,4])/abs(max(d[:,2])) * 100.0, "b", linewidth=3)
ax[1].plot(d[:,0], abs(d[:,2] - d[:,4])/max(abs(d[:,2]))* 100.0, "b", linewidth=3)
ax[1].plot(d[:,0], abs(d[:,1] - d[:,2])/max(abs(d[:,2])) * 100.0, "r", linewidth=3)
# ax[1].plot(d[:,0], abs(d[:,1] - d[:,6])/max(abs(d[:,2])) * 100.0, "k", linewidth=3)
ax[1].plot(d[:,0], abs(d[:,1] - d[:,6])/abs(d[:,2]) * 100.0, "k", linewidth=3)
# ax[1].plot(d[:,0],d[:,0]*d[:,0]*0.144,"g--",linewidth=3)
ax[1].legend(['Rel. err. nonadvec. (%)', 'Rel. err. nonlin', 'Rel. err. bp'],fontsize = MEDIUM_SIZE)
ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
t = ax[1].yaxis.get_offset_text()
t.set_size(MEDIUM_SIZE)
ax[1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)
# ax[1].plot(d[:,0], np.log(abs(d[:,2] - d[:,4])/abs(max(d[:,2])) * 100.0), "b", linewidth=3)
# ax[1].plot(d[:,0], np.log(abs(d[:,1] - d[:,2])/abs(max(d[:,2])) * 100.0), "r", linewidth=3)
# ax[1].plot(d[:,0], np.log(abs(d[:,1] - d[:,6])/abs(max(d[:,2])) * 100.0), "k", linewidth=3)
# ax[1].legend(['Relative error (%)', 'Rel err of nonlin', 'Rel err of sph'],fontsize = MEDIUM_SIZE)
# ax[1].tick_params(axis='both', which='major', labelsize=MEDIUM_SIZE)
# t = ax[1].yaxis.get_offset_text()
# t.set_size(MEDIUM_SIZE)
# ax[1].set_xlabel("Radius (scaled to Earth's)", fontsize = BIGGER_SIZE)


# plt.ylim(-0.00001, 0.00001)
# manager = plt.get_current_fig_manager()
# manager.full_screen_toggle()
# figure = plt.gcf()  
# figure.set_size_inches(9, 9) 
# figure = plt.gcf()  
# figure.set_size_inches(9, 12) 
# plt.savefig("../../../Notes/New_coupling/figures/Benchmark_2.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
# plt.savefig("Big_SolASPH.pdf",format = "pdf", bbox_inches='tight', pad_inches = 0)
plt.show()
