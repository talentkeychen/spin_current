# This is the script for plotting the spin
# texture and spin currents.
# Here we will add some noise to the order
# parameter
#
# Tian-Qi Chen, 20/11/2016


import numpy as np
from scipy.integrate import dblquad, quad
from scipy.stats import linregress
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 8



# Parameter define
gamma = 0                                # the strength of the noise
sigma_2 = 10                             # the sigma squared in the function of noise
sigma_1 = 10000                             # the sigma squared in the exponential of all Chi's
sigma_3 = 10000                             # the sigma squared in Chi_1 and Chi_{-1}
sigma_4 = 10000                             # the sigma squared in Chi_0
C_1 = 1/np.sqrt(2*np.pi*sigma_3)
C_2 = 1/np.sqrt(2*np.pi*sigma_4)
x_0 = np.array([2, -2])
y_0 = np.array([0, 0])                     # the coordinates of the centers of the noise
N = len(x_0)                               # the length of the coordinates array
N_f = 3                                    # the winding number of f(x,y)
N_g = 2                                    # the winding number of g(x,y)
NN = (N_f**2-N_g**2)/N_g                   # the big N
L = 200                                   # the length of the density window
yl_i = 0                                   # the initial value of yl
yl_f = L/2                                 # the final value of yl
yl = np.linspace(yl_i, yl_f, 10)           # the y position of placing single vortex when measuring
Nl = len(yl)                               # the length of yl
h_const = 1                                # h_bar over m
N_rho = 100                                # the N for the density function rho
delta = 0                                  # the length of the laser beam

# Define functions


def devx(x, y, t):                         # Derivative with respect to x
    tmp_1 = 0
    tmp_2 = 0                              # initialization
    rslt = 0
    if t == 'f':
        rslt = -y/(x**2+y**2)
    elif t == 'g':
        for i in range(N):
            tmp_1 = tmp_1+np.exp(-((x-x_0[i])**2+(y-y_0[i])**2)/(2*sigma_2))
        for j in range(N):
            tmp_2 = tmp_2+(-1)*(x-x_0[j])/sigma_2 *np.exp(-((x-x_0[j])**2+(y-y_0[j])**2)/(2*sigma_2))
        rslt = (-y/(x**2+y**2))*(1+gamma*tmp_1)+np.arctan(y/x)*gamma*tmp_2
    return rslt

def rho(y, yk):
    return np.sqrt((2*C_1**2+C_2**2)*np.exp(-(y+yk)**2/sigma_1))


def KdotF1(y, yk):
    return 2*np.sqrt(2/3)*C_1*C_2*np.exp(-(y+yk)**2/(1*sigma_1))

def KdotF2(y, yk):
    return 2/3 * (2*C_1**2+C_2**2) * np.exp(-(y+yk)**2/(1*sigma_1))

def integralResult(yl):
    jc1, err_jc1 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') + devx(x, y, 'g') * KdotF1(y, yl)), -np.Inf, 0, lambda x: -np.Inf, lambda x: np.Inf)
    jc2, err_jc2 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') + devx(x, y, 'g') * KdotF1(y, yl)), 0, np.Inf, lambda x: -np.Inf, lambda x: np.Inf)
    j11, err_j11 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') * KdotF1(y, yl) + devx(x, y, 'g') * KdotF2(y, yl)), -np.Inf, 0, lambda x: -np.Inf, lambda x: np.Inf)
    j12, err_j12 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') * KdotF1(y, yl) + devx(x, y, 'g') * KdotF2(y, yl)), 0, np.Inf, lambda x: -np.Inf, lambda x: np.Inf)
    j21, err_j21 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') * KdotF2(y, yl) + devx(x, y, 'g') * KdotF2(y, yl)), -np.Inf, 0, lambda x: -np.Inf, lambda x: np.Inf)
    j22, err_j22 = dblquad(lambda x,y: rho(y, yl) * (devx(x, y, 'f') * KdotF2(y, yl) + devx(x, y, 'g') * KdotF2(y, yl)), 0, np.Inf, lambda x: -np.Inf, lambda x: np.Inf)
    return (NN * (jc1+jc2) - (N_f * (j11+j12) - N_g * (j21+j22)))/np.pi


def integralAY(yl):
    val1, abserr1 = quad(lambda y: rho(y, 0), -np.Inf, yl)
    val2, abserr2 = quad(lambda y: rho(y, 0), yl, np.Inf)
    return val1 - val2


# Plot the Jx-Ay relation

Jtot = []    # the total current measured
Ay = []   # the A_y

print('Start the integration')
for j in range(Nl):
    #print('Start the integration')
    Jtot.append(integralResult(yl[j]))
    print('Integration completed for yl=', yl[j])
    Ay.append(integralAY(yl[j]))


# Linear fitting

z = np.polyfit(Ay, Jtot, 1)

plt.plot(Ay, Jtot, 'bo')

# Plot the fitted relation

xx = np.linspace(min(Ay), max(Ay), 50)
yy = z[0]*xx+z[1]
plt.plot(xx, yy, 'r')
#print(xx)

# Correlation coefficients calculation


def rsquared(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return slope, r_value**2

rsl1 = rsquared(Ay, Jtot)
print('The correlation coefficient is: ', rsl1[1])
print('The slope is: ', rsl1[0])

# Show the fits to data and fitted relation

#SIZE = 20
#plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
#plt.rc('axes', labelsize=SIZE)  # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
#plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
plt.xticks([0.0, 0.4, 0.8, 1.2], fontsize=20)
plt.yticks([1.0, 2.0, 3.0], fontsize=20)
plt.xlabel(r'$A_y$', fontsize=20)
plt.ylabel(r'$J_{tot}/G_{0}$', fontsize=20)
#plt.title(r'$\sigma_1=$'+str(sigma_1)+', '+r'$\sigma_4=$'+str(sigma_4)+', '+r'$\sigma_3=$'+str(sigma_3))
#plt.grid(True)
#plt.figure(figsize=(60,60))
plt.savefig('Jtot(Nf='+str(N_f)+', Ng='+str(N_g)+').eps', format='eps')
plt.show()

