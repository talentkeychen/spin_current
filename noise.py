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
from pylab import rcParams
rcParams['figure.figsize'] = 10, 8


# Here we plot the spin superfluid current with noice

# First we try to test how to perform partial derivative
# in python using sympy
#X = Symbol('X')
#Y = Symbol('Y')
#f_X = diff(1/(X**2+Y**2), X)
#g_Y = diff(1/(X**2+Y**2), Y)

#print(f_X)

# Parameter define
gamma = 1.00  # the strength of the noise
sigma_2 = 100                             # the sigma squared in the function of noise
sigma_1 = 10000                             # the sigma squared in the exponential of all Chi's
sigma_3 = 10000                             # the sigma squared in Chi_1 and Chi_{-1}
sigma_4 = 10000                             # the sigma squared in Chi_0
C_1 = 1/np.sqrt(2*np.pi*sigma_3)
C_2 = 1/np.sqrt(2*np.pi*sigma_4)
N_f = 1                                     # N_f winding number
N_g = 2                                     # N_g winding number
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
def rho(y):
    return np.sqrt((2*C_1**2+C_2**2)*np.exp(-y**2/sigma_1))

def KdotF1(y):
    return 2*np.sqrt(2/3)*C_1*C_2*np.exp(-y**2/(1*sigma_1))

def KdotF2(y):
    return 2/3 * (2*C_1**2+C_2**2) * np.exp(-y**2/(1*sigma_1))

# Define the Hermite polynomial function for n=2
#def Her(x, y, a, b, s):
#    rst = 0
#    if s == 0:
#        rst = 1
#    elif s == 1:
#        rst = a*(x**2+y**2) - b
#    else:
#        print('The parameter s is wrong.')
#    return rst



# Plot the spin superfluid current density
Y, X = np.mgrid[-5:5:20j, -5:5:20j] # X, Y are the plot starting points
#U = (devx(X, Y)*N_f + KdotF1(Y)*devx(X, Y)*N_g) * rho(Y)       # the conventional current
#V = (devy(X, Y)*N_f + KdotF1(Y)*devy(X, Y)*N_g) * rho(Y)

#U = (devx(X, Y)*N_f*KdotF1(Y) + KdotF2(Y)*devx(X, Y)*N_g) * rho(Y)       # the first-order spin current
#V = (devy(X, Y)*N_f*KdotF1(Y) + KdotF2(Y)*devy(X, Y)*N_g) * rho(Y)

U = (devx(X, Y)*N_f*KdotF2(Y) + KdotF1(Y)*devx(X, Y)*N_g) * rho(Y)       # the second-order spin current
V = (devy(X, Y)*N_f*KdotF2(Y) + KdotF1(Y)*devy(X, Y)*N_g) * rho(Y)


#U = devx(X, Y)   # The situation where Chi(r) is an constant
#V = devy(X, Y)

speed = np.sqrt(U**2 + V**2)
#speed = U*2+V
UN = U/speed
VN = V/speed
plot1 = plt.figure()
plt.quiver(X, Y, UN, VN,        # data
           speed*100000,                   # colour the arrows based on this array
           cmap=cm.seismic,     # colour map
           headlength=5)        # length of the arrows

plt.colorbar()                  # adds the colour bar

#plt.title('Spin superfluid current with noise')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$X$', fontsize=20)
plt.ylabel(r'$Y$', fontsize=20)
plt.savefig('Js2Plot(Noise='+str(gamma)+').eps', format='eps')
plt.show(plot1)                 # display the plot



