# Simple turbulent channel case at Re_tau=180 with 5832 elements and N=5 
Simulation of a turbulent channel in a domain of size $(4\pi\delta,2\delta,4/3\pi\delta)$  where $\delta$ is the height of half then channel. The bulk Reynolds number is specified ($U_b=1$) with constant mass flux s.t. $Re_b=\delta U_b/\nu=2800$. The mesh is deformed and refined closer to the wall in turb_channel.f90 (the original box.nmsh is of size (4, 2, 1.5) and centered in (2, 0, 0.75)).

The channel starts with a slightly perturbed initial condition which becomes turbulent after around 10 time units. 

This case can be run if one has 4GB of DRAM and a relatively modern CPU. On a Macbook Air 2024 M3 running to $T=100\delta/U_b$ takes around 1 hour. On one RTX2060 GPU it takes less than 10 min instead. In single precision the case can be run faster and the memory requirement is lower as well (enabled by configuring with --enable-real=sp).

# Enabled features in turb_channel.case 

- The force exerted on the walls of the channel is calculated and normalized by the channel area to obtain the wall-shear stress. This is printed in the log for the top and bottom wall every 10 time steps. This can be modified in the case file as well. One can for example extract the stress in the stream-wise direction as (assuming one saves the log in my_log) and put it into a textfile out:
```bash
  awk '/forcex/ {print($1,$2,$3)} ' my_log > out #take the time step, time and stress and put into out
```
One can then plot the evolution of $Re_\tau$ in python with something simple like:
```python
import numpy as np
import matplotlib.pyplot as plt
Re_b = 2800 #Because we are non-dimensional this is ok
dat = np.genfromtxt('out') # Read out into a numpy array
plt.plot(dat[::2,1],np.sqrt(dat[::2,2])*Re_b,label='Bottom wall') 
plt.plot(dat[1::2,1],np.sqrt(dat[1::2,2])*Re_b,'--',label='Top wall') 
plt.ylabel(r'$Re_\tau$')
plt.xlabel(r'Time, $\delta/U_b$')
plt.legend()
plt.show()
```

- All statistics required to compute the budgets for the Reynolds stress components are sampled from $T=60\delta/U_b$ and averaged in the homogenous directions. The output is then written into a CSV file called fluid_stats.csv. This file can then be easily read into python or preffered postprocessing tool of choice. For example, if one wants compute the mean velocity profile, one ccould load `fluid_stats.csv` into python accordin to the following:

```python
import numpy as np
import matplotlib.pyplot as plt
# dat = [output time, coordinate, <p>, <u>, <v>, <w>, <pp>, <uu>, <vv>, <ww>, <uv>, <uw>,... 
dat = np.genfromtxt('fluid_stats.csv',delimiter=',')
#Time of batch of interest (our batchsize is 5 convective time units) and sampling started at 60
#First batch of average between 60 and 65 is therefore written out at T=65 (OBSERVE you need to have run beyond T=65)
# Lets see how the mean profile looks
batch_time = 65 # Can be changed 
#Extract <u>
U_vel = dat[np.abs(dat[:,0]-65)<0.1,3]
#coordinates
y_coords = dat[np.abs(dat[:,0]-65)<0.1,1]
#plot the profile
plt.plot(y_coords,U_vel)
#If you run longer you can compute better averages by adding many batches together, for somewhat converged statistics you will need at least 100 time units.
plt.ylabel(r'$u$ ($U_b$)')
plt.xlabel(r'$y$ ($\delta$)')
plt.legend()
plt.show()
```
Details on all the averages that are computed can be found in the documentation under statistics_guide.

- Lambda 2 is also calculated. It can be nicely visualized in paraview or visit by creating an isosurfaces of the "temperature" field.

