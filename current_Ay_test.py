import numpy as np
from scipy.integrate import dblquad, quad
from scipy.stats import linregress
#from decimal import Decimal
import matplotlib.pyplot as plt


# Parameter define
gamma = 0.8    # the strength of the noise
sigma_2 = 1    # the sigma square in the function of noise
sigma_1 = 1    # the sigma square in the initial state Chi(r)
sigma_3 = 1000000    # the sigma square of the density rho(y)
sigma_4 = 1    # another sigma square in the initial state Chi(r)
x_0 = np.array([1, -2])
y_0 = np.array([0, 0])     # the coordinates of the centers of the noise
N = len(x_0)   # the length of the coordinates array
N_f = 6        # the winding number of f(x,y)
N_g = 1        # the winding number of g(x,y)
NN = (N_f**2-N_g**2)/N_g  # the big N
L = 2000          # the length of the density window
yl_i = 0      # the initial value of yl
yl_f = L/2     # the final value of yl
yl = np.linspace(yl_i, yl_f, 10)  # the y position of placing single vortex when measuring
Nl = len(yl)   # the length of yl
h_const = 1      # h_bar over m
N_rho = 100    # the N for the density function rho
delta = 0    # the length of the laser beam
# Define functions


def devx(x, y, t):     # Derivative with respect to x
    tmp_1 = 0
    tmp_2 = 0          # initialization
    rslt = 0
    if t == 'f':       # without noise
        rslt = -y/(x**2+y**2)
    elif t == 'g':     # with noise
        for i in range(N):
            tmp_1 = tmp_1+np.exp(-((x-x_0[i])**2+(y-y_0[i])**2)/(2*sigma_2))
        for j in range(N):
            tmp_2 = tmp_2+(-1)*(x-x_0[j])/sigma_2 *np.exp(-((x-x_0[j])**2+(y-y_0[j])**2)/(2*sigma_2))
        rslt = (-y/(x**2+y**2))*(1+gamma*tmp_1)+np.arctan(y/x)*gamma*tmp_2
    return rslt


def rho(y, t):# the density of order parameter
    if t == 0:
        return 1
    elif t == 1:
        return (N_rho/np.sqrt(2*np.pi*sigma_3)) * np.exp(-1*y**2/(2*sigma_3))
t =1


def rho1(y, yk):
    return (N_rho/np.sqrt(2*np.pi*sigma_3)) * np.exp(-1*(y+yk)**2/(2*sigma_3))

def initialchi(x, y):
    return 1/np.sqrt(2*np.pi*sigma_4)*np.exp(-x**2/(2*sigma_4)-y**2/(2*sigma_4))

#print(devx(100, 200, 'f'))

# Calculate the integral


def integralResult(N, yl):
    val1, abserr1 = dblquad(lambda x, y: rho1(y, yl) * N * devx(x, y, 'g'), -np.Inf, 0, lambda x: -np.Infinity,lambda x: np.Infinity)
    val2, abserr2 = dblquad(lambda x, y: rho1(y, yl) * N * devx(x, y, 'g'), 0, np.Inf, lambda x: -np.Infinity,lambda x: np.Infinity)
    return (val1 + val2) / np.pi

    #val1, abserr1 = dblquad(lambda x, y: rho(y, t) * N_f * devx(x, y, 'f'), -np.Inf, 0, lambda x: -np.Infinity, lambda x: np.Infinity)
    #val2, abserr2 = dblquad(lambda x, y: rho(y, t) * N_f * devx(x, y, 'f'), 0, yl-delta, lambda x: -np.Infinity, lambda x: np.Infinity)
    #val3, abserr3 = dblquad(lambda x, y: rho(y, t) * N_f * devx(x, y, 'f'), yl+delta, np.Inf, lambda x: -np.Infinity, lambda x: np.Infinity)


    #val1, abserr1 = dblquad(lambda x, y: rho1(y, yl) * N_f * devx(x, y, 'f'), -np.Inf, 0, lambda x: -np.Infinity, lambda x: np.Infinity)
    #val2, abserr2 = dblquad(lambda x, y: rho1(y, yl) * N_f * devx(x, y, 'f'), 0, np.Inf, lambda x: -np.Infinity, lambda x: np.Infinity)


    #val1, abserr1 = dblquad(lambda x, y: rho(y, t) * N_f * devx(x, y, 'f'), -L-yl, 0, lambda x: -np.Infinity, lambda x: np.Infinity)
    #val2, abserr2 = dblquad(lambda x, y: rho(y, t) * N_f * devx(x, y, 'f'), 0, L-yl, lambda x: -np.Infinity, lambda x: np.Infinity)


def integralAY(yl):
    val1, abserr1 = quad(lambda y: rho(y, t), -np.Inf, yl)
    val2, abserr2 = quad(lambda y: rho(y, t), yl, np.Inf)
    return val1 - val2

#############Test Area##############################
#a, b = dblquad(lambda x, y: rho(y) * N_f * devx(x, y, 'f'), -np.Inf, 0, lambda x: -np.Inf, lambda x: np.Inf)
#c, d = dblquad(lambda x, y: rho(y) * N_f * devx(x, y, 'f'), 0, 10, lambda x: -np.Inf, lambda x: np.Inf)

#print(a+c)


######################################################

# Plot the Jx-Ay relation
Jtot1 = []    # the total current measured
Jtot2 = []
Jtot3 = []
Ay = []   # the A_y
print('Start the integration')
for j in range(Nl):
    #print('Start the integration')
    Jtot1.append(integralResult(1, yl[j]))
    Jtot2.append(integralResult(2, yl[j]))
    Jtot3.append(integralResult(3, yl[j]))
    print('Integration completed for yl=', yl[j])
    Ay.append(integralAY(yl[j]))
#print(Jtot)
#print(Ay)

# Linear fitting
z1 = np.polyfit(Ay, Jtot1, 1)
z2 = np.polyfit(Ay, Jtot2, 1)
z3 = np.polyfit(Ay, Jtot3, 1)
#print(z1)
#print(yl)
plt.plot(Ay, Jtot1, 'ro')
plt.plot(Ay, Jtot2, 'go')
plt.plot(Ay, Jtot3, 'bo')

# Plot the fitted relation
xx = np.linspace(min(Ay), max(Ay), 50)
yy1 = z1[0]*xx+z1[1]
yy2 = z2[0]*xx+z2[1]
yy3 = z3[0]*xx+z3[1]
plt.plot(xx, yy1, 'r')
plt.plot(xx, yy2, 'g')
plt.plot(xx, yy3, 'b')

# Correlation coefficients calculation


def rsquared(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return slope, r_value**2

rsl1 = rsquared(Ay, Jtot1)
rsl2 = rsquared(Ay, Jtot2)
rsl3 = rsquared(Ay, Jtot3)
print('The correlation coefficient is: ', rsl1[1], rsl2[1], rsl3[1] )
print('The slope is: ', rsl1[0], rsl2[0], rsl3[0])

# Show the fits to data and fitted relation
plt.xlabel(r'$A_y$')
plt.ylabel(r'$J_x/G_{0}$')
plt.text(40, 150, 'N=3', color='blue')
plt.text(60, 100, 'N=2', color='green')
plt.text(60, 45, 'N=1', color='red')
plt.grid(True)
plt.show()