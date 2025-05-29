# Postprocessing script for the lid-driven cavity

import numpy as np
import matplotlib.pyplot as plt 

data1 = np.genfromtxt('ekin.csv', delimiter=',')

plt.figure()
plt.plot(data1[:,0],data1[:,1],label=r'$Re=4000$')
plt.xlabel(r'time $t$')
plt.ylabel(r'kin. energy $E$')
plt.legend()
plt.show()

plt.figure()
plt.plot(data1[:,0],data1[:,2],label=r'$Re=4000$')
plt.xlabel(r'time $t$')
plt.ylabel(r'enstrophy $En$')
plt.legend()
plt.show()
