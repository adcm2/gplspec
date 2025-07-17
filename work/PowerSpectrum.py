 #%%
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
import math

# import data
# path0 = "Phobos/Cartesian/MatrixSolution.out"
# path1 = "Phobos/Spherical/MatrixSolution.out"
# path2 = "Phobos/Rotated/MatrixSolution.out"

path = "Phobos/Spherical/PowerSpectrum.out"
dphobos = np.loadtxt(path, delimiter=";")

row,col = dphobos.shape
print("Num rows is: ", row)
print("Num cols is: ", col)

lval = np.arange(0,col)

width = 0.45
fig, ax = plt.subplots(1,1)
ax.bar(lval,dphobos[60,:],width,color='b')
ax.bar(lval+width,dphobos[row-1,:],width,color='r')
ax.set_yscale('log')
ax.tick_params(axis='both',labelsize=15)
ax.set_xlabel('l',fontsize=20)
ax.set_ylabel('Normalized power',fontsize=20)
plt.show()
