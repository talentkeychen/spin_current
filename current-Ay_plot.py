import numpy as np
from scipy.integrate import dblquad, quad
from scipy.stats import linregress
#from decimal import Decimal
import matplotlib.pyplot as plt


# Parameter define
gamma = 0.0    # the strength of the noise
sigma_2 = 1    # the sigma square in the function of noise
sigma_1 = 1    # the sigma square in the initial state Chi(r)
sigma_3 = 1000    # the sigma square of the density rho(y)
sigma_4 = 1    # another sigma square in the initial state Chi(r)
x_0 = np.array([1, -2])
y_0 = np.array([0, 0])     # the coordinates of the centers of the noise
N = len(x_0)   # the length of the coordinates array
N_f = 3        # the winding number of f(x,y)
N_g = 1        # the winding number of g(x,y)
NN = (N_f**2-N_g**2)/N_g  # the big N
yl_i = 35     # the initial value of yl
yl_f = 40     # the final value of yl
yl = np.linspace(yl_i, yl_f, 10)  # the y position of placing single vortex when measuring
Nl = len(yl)   # the length of yl
h_const=1      # h_bar over m
L = 40
# Define functions
def devx(x, y, t): # Derivative with respect to x
    tmp_1 = 0
    tmp_2 = 0  # initialization
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

def rho(y): # the density of order parameter
    return (1/np.sqrt(2*np.pi*sigma_3)) * np.exp(-y**2/(2*sigma_3))

def initialchi(x, y):
    return 1/np.sqrt(2*np.pi*sigma_4)*np.exp(-x**2/(2*sigma_4)-y**2/(2*sigma_4))

#print(devx(100, 200, 'f'))

# Calculate the integral

def integralResult(yl):
    #val, abserr = dblquad(lambda x, y: rho(y)*h_const*N_g*devx(x, y, 'f')*2*initialchi(x, y), yl_i-yl, yl_f-yl, lambda x: -2, lambda x: 2)
    val1, abserr1 = dblquad(lambda x, y: NN*rho(y) * N_f * devx(x, y, 'f'), -L, 0, lambda x: -np.Inf, lambda x: np.Inf)
    val2, abserr2 = dblquad(lambda x, y: NN*rho(y) * N_f * devx(x, y, 'f'), 0, L-yl, lambda x: -np.Inf, lambda x: np.Inf)
    return (val1+val2)/np.pi

def integralAY(yl):
    val1, abserr1 = quad(lambda y: rho(y), -np.Inf, yl)
    val2, abserr2 = quad(lambda y: rho(y), yl, np.Inf)
    return val1 - val2

Jtot =[]                               # the total current measured
Ay = []                                # the A_y
for j in range(Nl):
    Jtot.append(integralResult(yl[j]))
    Ay.append(integralAY(yl[j]))
print(Jtot)
print(Ay)

# Linear fitting
z1 = np.polyfit(Ay, Jtot, 1)
print(z1)
plt.plot(Ay, Jtot, 'ro')

# Plot the fitted relation
xx=np.linspace(min(Ay),max(Ay),50)
yy=z1[0]*xx+z1[1]
plt.plot(xx, yy)

# Show the fits to data and fitted relation
plt.show()

# Correlation coefficients calculation
def rsquared(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return r_value**2

print(rsquared(Ay, Jtot))
