# This is the script for plotting the spin
# texture and spin currents
# Here we will calculate the series expansion
# for sin(x) up to an arbitrary order N.
#
# Tian-Qi Chen, 18/11/2016

#from scipy.special import jn
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy as np

# Here we plot the spin superfluid current

Y, X = np.mgrid[-2:2:20j, -2:2:20j]  # X, Y are the plot starting points
U = -Y/(X**2+Y**2)*np.exp(-0.5*(X**2+Y**2))
V = X/(X**2+Y**2)*np.exp(-0.5*(X**2+Y**2))
speed = np.sqrt(U**2 + V**2)
UN = U/speed
VN = V/speed
plot1 = plt.figure()
plt.quiver(X, Y, UN, VN,        # data
           speed,                   # colour the arrows based on this array
           cmap=cm.seismic,     # colour map
           headlength=6)        # length of the arrows

plt.colorbar()                  # adds the colour bar

#plt.title('Spin superfluid current')
plt.show(plot1)                 # display the plot

# Here we add some random noice

