# This is the script for plotting the spin
# texture and spin currents.
# Here we will add some noise to the order
# parameter
#
# Tian-Qi Chen, 20/11/2016

# from scipy.special import jn
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
gamma = 0.0  # the strength of the noise
sigma_2 = 1    # the sigma square in the function of noise
sigma_1 = 1    # the sigma square in the initial state Chi(r)
x_0 = np.array([1, -2])
y_0 = np.array([0, 0])     # the coordinates of the centers of the noise
N = len(x_0)  # the length of the coordinates array

# Define functions to calculate partial derivatives
def devx(x, y):
    tmp_1 = 0
    tmp_2 = 0  # initialization
    for i in range(N):
        tmp_1 = tmp_1+np.exp(-((x-x_0[i])**2+(y-y_0[i])**2)/(2*sigma_2))
    for j in range(N):
        tmp_2 = tmp_2+(-1)*(x-x_0[j])/sigma_2 *np.exp(-((x-x_0[j])**2+(y-y_0[j])**2)/(2*sigma_2))
    return (-y/(x**2+y**2))*(1+gamma*tmp_1)+np.arctan(y/x)*gamma*tmp_2

def devy(x, y):
    tmp_1 = 0
    tmp_2 = 0
    for i in range(N):
        tmp_1 = tmp_1+np.exp(-((x-x_0[i])**2+(y-y_0[i])**2)/(2*sigma_2))
    for j in range(N):
        tmp_2 = tmp_2+(-1)*(y-y_0[j])/sigma_2 *np.exp(-((x-x_0[j])**2+(y-y_0[j])**2)/(2*sigma_2))
    return (x/(x**2+y**2))*(1+gamma*tmp_1)+np.arctan(y/x)*gamma*tmp_2

# Define the y-dependent density function
def rho(y,sigma_square):
    return 1/np.sqrt(2*np.pi*sigma_square) *np.exp(-y**2/(2*sigma_square))

# Define the Hermite polynomial function for n=2
def Her(x, y, a, b, s):
    rst = 0
    if s == 0:
        rst = 1
    elif s == 1:
        rst = a*(x**2+y**2) - b
    else:
        print('The parameter s is wrong.')
    return rst



# Plot the spin superfluid current density
Y, X = np.mgrid[-3:3:20j, -3:3:20j] # X, Y are the plot starting points
U = devx(X, Y)*Her(X, Y, 1, 1, 0)*np.exp(-(X**2+Y**2)/(2*sigma_1)) * rho(Y, 1)
V = devy(X, Y)*Her(X, Y, 1, 1, 0)*np.exp(-(X**2+Y**2)/(2*sigma_1)) * rho(Y, 1)
#U = devx(X, Y)   # The situation where Chi(r) is an constant
#V = devy(X, Y)

speed = np.sqrt(U**2 + V**2)
#speed = U*2+V
UN = U/speed
VN = V/speed
plot1 = plt.figure()
plt.quiver(X, Y, UN, VN,        # data
           speed,                   # colour the arrows based on this array
           cmap=cm.seismic,     # colour map
           headlength=5)        # length of the arrows

plt.colorbar()                  # adds the colour bar

#plt.title('Spin superfluid current with noise')
plt.xlabel('X')
plt.ylabel('Y')
plt.show(plot1)                 # display the plot



