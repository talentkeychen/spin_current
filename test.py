import numpy as np
from scipy.integrate import dblquad, quad
def f(a):
    return dblquad(lambda y, x: x*y*a, 0, 1, lambda x: 0, lambda x: 2)
print(f(1))
