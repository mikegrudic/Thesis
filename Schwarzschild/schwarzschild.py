import numpy as np
import math
import scipy
from scipy import special, integrate, misc, optimize, interpolate
import time

pi=np.pi

def bmin(e):
    return math.sqrt((8 - 36*e**2 + 27*e**4 + e*(9*e**2 - 8)**(3.0/2.0))/2)/(e**2 - 1)

def roots(e, b):
    coeffs = (e**2 - 1, 2, -b**2*(e**2-1), 2*b**2*(e**2-1))
    return np.sort(np.roots(coeffs))

def deflection(e, b):
    r1, r2, r3 = roots(e,b)
    m = r3*(r2-r1)/(r2*(r3-r1))

    p2 = lambda t: 1/math.sqrt((1-t*t)*(1-m*t*t))
    t=time.time()
    ellipf2 = integrate.quad(p2, 0, math.sqrt(r2/r3))[0]

    return 4*b*ellipf2/math.sqrt(r2*(r3-r1)) - np.pi

def SolveForAngle(e, theta):
    f = lambda b: deflection(e,b) - theta
    return optimize.newton(f, 10*bmin(e))

grid = np.linspace(0, pi, 1000)
    
E = 1.1

b = bmin(E) + 1/np.logspace(-5,10,10000)

deflections = np.array([deflection(E, B) for B in b])
deflections_modpi = deflections % pi
ddef_db = np.gradient(deflections)/np.gradient(b)
sigma = np.abs(b/np.sin(deflections)) * np.abs(1/ddef_db)

max_branch = np.max(deflections) - np.max(deflections) % pi
num_branches = int(max_branch/pi + 0.5)
branch_bounds = np.linspace(pi, max_branch, num_branches)

branch_index = [(deflections > bound - pi)*(deflections < bound) for bound in branch_bounds]

branch_sigma = np.array([np.interp(grid, deflections_modpi[i], sigma[i]) for i in branch_index])

summed_sigma = np.cumsum(branch_sigma, axis=0)

total_sigma = 2*pi*integrate.simps(summed_sigma*np.sin(grid), grid)
sigma_error = np.abs(total_sigma - total_sigma[-1])

#np.savetxt("deflections.dat", np.array((b, deflections)).T)
np.savetxt("sigma.dat", np.vstack((grid, summed_sigma)).T)
np.savetxt("error.dat", sigma_error.T)
