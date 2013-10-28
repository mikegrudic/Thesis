import numpy as np
import math
from optparse import OptionParser
import CarlsonR
import scipy
from scipy import special, integrate, misc, optimize, interpolate

pi=np.pi

def bmin(e):
    if e == "NULL":
        return 3*math.sqrt(3)
    else:
        return math.sqrt((8 - 36*e**2 + 27*e**4 + e*(9*e**2 - 8)**(3.0/2.0))/2)/(e**2 - 1)

def roots(e, b):
    if e == "NULL":
        coeffs = (1.0, 0, -b**2, 2*b**2)
    else:
        coeffs = (e**2 - 1, 2.0, -b*b*(e**2 - 1), 2*b*b*(e**2 - 1.0))
    return np.sort(np.roots(coeffs))

def deflection(e, b):
    r1, r2, r3 = roots(e,b)

    x = (r3 - r1)*(r3 - r2)
    y = r3*(r3 - r2)
    z = r3*(r3 - r1)
    ellipf = CarlsonR.BoostRF(x, y, z)
    return 4*b*ellipf - np.pi

def newton_deflection(e, b):
    return pi - 2*math.acos(1/math.sqrt((e**2 - 1)**2 * b**2 + 1))

def SolveForAngle(e, theta):
    f = lambda b: deflection(e,b) - theta
    return optimize.newton(f, 10*bmin(e))

p=OptionParser()
p.add_option("--E", help="Ratio between the total and rest energy of the particle. Let E=NULL for null geodesics.")
(opts,args)=p.parse_args()

grid = np.linspace(0, pi, 10000)
    
E = opts.E
if E != "NULL":
    E = float(E)

b = bmin(E) + 1/np.logspace(-5,10,100000)

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

np.savetxt("diff_sigma.dat", np.vstack((grid/pi, summed_sigma)).T)
np.savetxt("branch_sigma.dat", np.vstack((grid/pi, branch_sigma)).T)
np.savetxt("sigma.dat", (total_sigma[-1],))
np.savetxt("error.dat", sigma_error.T)
