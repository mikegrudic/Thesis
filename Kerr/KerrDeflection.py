import numpy as np
import math
from scipy import weave
from CarlsonR import *

pi = np.pi

def CubicRoots(e, b):
    if type(b) != np.ndarray:
        b = np.array([b])
    r3 = np.zeros(b.shape)
    r1 = np.copy(r3)
    r2 = np.copy(r3)
    code = """                                                                                               
    int i;                                                                                                   
    double B, C, Q, R, Q3, theta, SQ;                                                                        
    double A = 2.0/(e + 1.0)/(e - 1.0);                                                                      
    double TAU = 2*3.141592653589793116;                                                                     
    for (i = 0; i < Nb[0]; ++i){                                                                             
        B = -b[i]*b[i];                                                                                      
        C = -2*B;                                                                                            
        Q = (A*A - 3*B)/9;                                                                                   
        R = (2*pow(A, 3) - 9*A*B + 27*C)/54.0;                                                               
        Q3 = pow(Q, 3);                                                                                      
        theta = acos(R/sqrt(Q3));                                                                            
        SQ = sqrt(Q);                                                                                        
        r1[i] = -2*SQ*cos(theta/3) - A/3;                                                                    
        r3[i] = -2*SQ*cos((theta + TAU)/3) - A/3;                                                            
        r2[i] = -2*SQ*cos((theta - TAU)/3) - A/3;                                                            
    }                                                                                                        
    """
    weave.inline(code, ['e', 'b', 'r3', 'r1', 'r2'])
    return r1, r2, r3

def EqDeflection(a, E, L):
    R_coeffs = (E**2-1,2,-a**2+E**2*a**2-L**2,2*E**2*a**2-4*E*a*L+2*L**2,0.0)
    roots = np.sort(np.roots(R_coeffs))
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



def KerrDeflection(a, theta, E, bx, by):
#    C1 = E**2 - 1
#    L = bx*C1*np.sin(theta)
#    Q = (E**2 - 1)*((bx**2 - a**2)*np.cos(theta) + by**2)
    
    if bx.shape != by.shape:
        print "Carter constants and angular momenta must have same dimension!"
        exit()
    
    eq_index = (theta==pi/2)*(Q==0.0)
    non_eq_index = np.invert(eq_index)
    ones = np.ones(L.shape)

    r_coeffs = np.array([C1*ones, 2*ones, a**2*C1 - L**2 - Q, 2*((-(a*E) + L)**2 + Q), -a**2*Q]).T
    mu_coeffs = np.array([C1*a**2*ones, C1*a**2 - L**2 - Q, Q]).T

    r_roots = np.sort(np.array([np.roots(c) for c in r_coeffs]), axis=1)
    mu_roots = np.sort(np.array([np.roots(c) for c in mu_coeffs]),axis=1)
    print mu_roots
    #Compute roots for mu
    a, b, c = C1*a**2, C1*a**2 - L**2 - Q, Q
    q = -0.5*(b + np.sign(b)*np.sqrt(b**2 - 4*a*c))
    mu1, mu2 = q/a, c/q
    print mu1
    print mu2
    exit()
 
    r1, r2, r3, r4 = r_roots[:,0], r_roots[:,1], r_roots[:,2], r_roots[:,3]
    r_integral = 2*InvSqrtQuartic(r1, r2, r3, r4, r4)
#    mu1, mu2 = np.sqrt(mu_roots[:,0]), np.sqrt(mu_roots[:,1])
#    print theta
#    print mu2
#    mu_big_integral = 2*InvSqrtQuartic(-mu1, -mu2, mu1, mu2, theta
