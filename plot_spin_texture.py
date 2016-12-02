from matplotlib.pyplot import cm
# from sympy import Symbol, diff
import matplotlib.pyplot as plt
import numpy as np


# Here we plot the spin superfluid current with noice

# First we try to test how to perform partial derivative
# in python using sympy
#X = Symbol('X')
#Y = Symbol('Y')
#f_X = diff(1/(X**2+Y**2), X)
#g_Y = diff(1/(X**2+Y**2), Y)

#print(f_X)

# Parameter define
gamma = 0.4  # the strength of the noise
sigma_2 = 1    # the sigma square in the function of noise
#sigma_1 = 1    # the sigma square in the initial state Chi(r)
x_0 = np.array([1, -4])
y_0 = np.array([0, 0])     # the coordinates of the centers of the noise
N = len(x_0)  # the length of the coordinates array
m = 1 # the index of magnetic quantum no.

#Define the functions to calculate the summation
def phasecal(x, y):
    tmp_phi = 0
    for j in range(N):
        tmp_phi = tmp_phi + np.exp(-(x-x_0[j])**2/sigma_2-(y-y_0[j])**2/sigma_2)
    return np.arctan(y/x) * (1 + gamma * tmp_phi)


# Plot the density hotmap
Y, X = np.mgrid[-4:4:70j, -4:4:70j]
Z = -m * phasecal(X, Y)
plot1 = plt.figure()

plt.pcolormesh(X, Y, Z)

plt.colorbar()                  # adds the colour bar

plt.xlabel('X')
plt.ylabel('Y')
plt.show(plot1)
