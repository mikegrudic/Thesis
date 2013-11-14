import numpy as np
import math
import scipy
from scipy import special, integrate, misc, optimize
from CarlsonR import *
import time

pi=np.pi

def deflection2(a, E, L):
    result = np.copy(L)
    R_coeffs = (E**2-1,2,-a**2+E**2*a**2-L**2,2*E**2*a**2-4*E*a*L+2*L**2,0.0)
    r_roots = np.sort(np.array([np.roots(c) for c in r_coeffs]), axis=1)
    
    if len(roots[np.nonzero(roots.imag)]) != 0:
        raise Exception("Particle falls into the black hole!")
    r1, r2, r3, r4 = roots
    C1 = math.sqrt(1-a**2)
    rplus, rminus = 1 + C1, 1 - C1

    int1 = InvSqrtQuartic(r1, r2, r3, r4, r4)
    int2 = TerribleIntegral(r1, r2, r3, r4, rplus, r4)
    int3 = TerribleIntegral(r1, r2, r3, r4, rminus, r4)

    part1 = (L - a*E)/math.sqrt(E**2 - 1) * int1
    part2 = -a*math.log((-rplus + r4)/(-rminus + r4))/(2.0*C1)
    part3 = - a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
    return 2*(part1 + part2 + part3) - pi

def deflection(a, E, L):
    R_coeffs = (E**2-1,2,-a**2+E**2*a**2-L**2,2*E**2*a**2-4*E*a*L+2*L**2,0.0)
    roots = np.sort(np.array([np.roots(c) for c in r_coeffs]), axis=1)
#    roots = np.sort(np.roots(R_coeffs))
    r3 = roots[:,-1]
#    m = r3*(r2-r1)/(r2*(r3-r1))

#    p2 = lambda t: 1/math.sqrt((1-t*t)*(1-m*t*t))
#    ellipf2 = integrate.quad(p2, 0, math.sqrt(r2/r3))[0]
    f = lambda r: (L-a*E)/math.sqrt(R(a, E, L, r))
    g = lambda r: a/(r**2 - 2*r + a**2)
    h = lambda r: -a/(r**2 - 2*r + a**2) * (E*(r**2 + a**2) - L * a)/math.sqrt(R(a, E, L, r)) 
    part1 = integrate.quad(f, r3, np.inf)[0]
    part2 = integrate.quad(g, r3, np.inf)[0]
    part3 = integrate.quad(h, r3, np.inf)[0]
#    print part1, part2, part3
    return 2*(part1 + part2 + part3) - pi

def R(a, E, L, r):
#    print np.polyval(R_coeffs, r)
    R_coeffs = (E**2-1,2,-a**2+E**2*a**2-L**2,2*E**2*a**2-4*E*a*L+2*L**2,0.0)
    return np.polyval(R_coeffs, r)

def dr_dphi(phi, r):
    sqrt_R = math.sqrt(np.abs(R(a, E, L, r)))
    return -sqrt_R/((L - a*E) + a*(sqrt_R - E*(r**2 + a**2) + L*a)/(r**2 + 2*r + a**2))

############################################################################3

#E, a, L, C = 1.1, 0.99, 6, -7

print deflection2(0.5, 1.1, 6)
print deflection(0.5, 1.1, 6)
