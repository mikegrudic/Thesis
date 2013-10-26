import numpy as np
import math
import scipy
from scipy import special, integrate, misc, optimize
import time

pi=np.pi

def deflection(e, b):
    r3 = roots[-1]
#    m = r3*(r2-r1)/(r2*(r3-r1))

#    p2 = lambda t: 1/math.sqrt((1-t*t)*(1-m*t*t))
#    ellipf2 = integrate.quad(p2, 0, math.sqrt(r2/r3))[0]
    f = lambda r: (L-a*E)/math.sqrt(R(r)) + a/(r**2 + 2*r + a**2) * (1 - (E*(r**2 + a**2) - L * a)/math.sqrt(R(r)))
    return 2*integrate.quad(f, r3, np.inf)[0] - np.pi

def R(r):
#    print np.polyval(R_coeffs, r)
    return np.polyval(R_coeffs, r)

def dr_dphi(phi, r):
    sqrt_R = math.sqrt(np.abs(R(r)))
    return -sqrt_R/((L - a*E) + a*(sqrt_R - E*(r**2 + a**2) + L*a)/(r**2 + 2*r + a**2))

############################################################################3

E, a, L, C = 1.1, 0.99, 6, -7

R_coeffs = (E**2-1,2,-a**2+E**2*a**2-L**2-C,2*E**2*a**2-4*E*a*L+2*L**2+2*C,-a**2*C)

roots = np.sort(np.roots(R_coeffs))

#phi = np.linspace(0, (deflection(E, L/math.sqrt(E**2-1)) + np.pi), 1999)
r = 2*roots[-1] + 1.0/np.logspace(-3, 3, 1000)



#x, y = r*np.cos(phi), r*np.sin(phi)
#np.savetxt("traj.dat", np.array([x,y]).T)

