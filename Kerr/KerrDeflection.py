import numpy as np
import math
from scipy import weave, integrate
import CarlsonR
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
    part2 = -a*np.log((-rplus + r4)/(-rminus + r4))/(2.0*C1)
    part3 = - a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
    return 2*(part1 + part2 + part3) - pi

def PhiTerribleIntegral(r1, r2, r3, r4, a, E, L):
    C1 = math.sqrt(1-a**2)
    rplus, rminus = 1 + C1, 1 - C1
    int1 = InvSqrtQuartic(r1, r2, r3, r4, r4)
    int2 = TerribleIntegral(r1, r2, r3, r4, rplus, r4)
    int3 = TerribleIntegral(r1, r2, r3, r4, rminus, r4)
    part2 = -a*np.log((-rplus + r4)/(-rminus + r4))/(2.0*C1)
    part3 = - a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
    return 2*(part2 + part3)

def SchwarzDeflection(E, b):
    bm = bmin(e)
    result = np.zeros(b.shape)
    fall_in= (b<=bm)
    result[fall_in] = 0.0
    r1, r2, r3 = Roots2(e,b[b>bmin(e)])
    x = (r3 - r1)*(r3 - r2)
    y = r3*(r3 - r2)
    z = r3*(r3 - r1)
    ellipf = CarlsonR.BoostRF(x, y, z)
    result[np.invert(fall_in)] = 4*b[np.invert(fall_in)]*ellipf - np.pi
    return result

def KerrDeflection(a, theta, E, bx, by):
#    if a==0:
#        deflection = SchwarzDeflection(E, np.sqrt(bx**2 + by**2))
        
    mu0 = np.cos(theta)
    s_mu = np.zeros(bx.shape, dtype=np.int64)
    if theta != 0 and theta != pi:
        s_mu = np.sign(-(2*bx*mu0+(-1 + 3*np.cos(2*theta))*by)*np.sin(theta))
    elif theta == 0:
        s_mu = -np.ones(by.shape)
    else:
        s_mu = np.ones(by.shape)
    
    C1 = E**2 - 1
    C2 = math.sqrt(1-a**2)
    ones = np.ones(bx.shape)
    L = bx*math.sqrt(C1)*math.sin(theta)
    Q = (E**2 - 1)*((bx**2 - a**2)*mu0**2 + by**2)

#Solve r and mu polynomials    
    r_coeffs = np.array([C1*ones, 2*ones, a**2*C1 - L**2 - Q, 2*((-(a*E) + L)**2 + Q), -a**2*Q]).T
    mu_coeffs = np.array([C1*a**2*ones, C1*a**2 - L**2 - Q, Q]).T
# r quartic
    r_roots = np.sort(np.array([np.roots(c) for c in r_coeffs]), axis=1)
#mu biquadratic
    A = -a**2
    B = a**2 * (1 + mu0**2) - bx**2 - by**2
    C = by**2 + (bx**2 - a**2)*mu0**2
    q = -0.5*(B + np.sign(B)*np.sqrt(B**2 - 4*A*C))
    M1, M2 = np.sort((q/A, C/q),axis = 0)
    mu_max = np.sqrt(M2)
    mu_min = -np.sqrt(M2)
    kSqr = M2/(M2 - M1)
    xSqr = np.abs(1 - mu0**2/M2)
    n = M2/(1-M2)
    A = a*np.sqrt(M2 - M1)

#do complete integrals
    r1, r2, r3, r4 = r_roots[:,0], r_roots[:,1], r_roots[:,2], r_roots[:,3]
    r_integral = 2*InvSqrtQuartic(r1, r2, r3, r4, r4)
    mu_complete_integral = 2*CarlsonR.BoostRF(np.zeros(kSqr.shape), 1.0-kSqr, np.ones(kSqr.shape))/A

#do partial mu integral from mu0 to mu_max    
    case1 = (mu0 != mu_max)*(mu0 != mu_min)
    mu_initial_integral = np.empty(bx.shape)
    if mu0 > 0:
        U1 = (-mu0**4 + (M1[case1]+M2[case1])*mu0**2 + -M1[case1]*M2[case1] + (mu0**2 - mu_max[case1]**2)**2)/(mu_max[case1]-mu0)**2
        mu_initial_integral[case1] = 2*CarlsonR.BoostRF(U1, U1 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U1 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1]))/a
    else:
    #break into mu0 to 0 + 0 to mu_max
        U2 = (-M1[case1]*M2[case1] + mu_max[case1]**4)/mu_max[case1]**2
        U3 = (np.sqrt((-mu0**4 + (M1[case1]+M2[case1])*mu0**2 - M1[case1]*M2[case1]))+np.sqrt(-M1[case1]*M2[case1]))**2/mu0**2
        mu_initial_integral[case1] = 2*(CarlsonR.BoostRF(U2, U2 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U2 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1])) + CarlsonR.BoostRF(U3, U3 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U3 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1])))/a

    mu_initial_integral[np.invert(case1)] = 0.0 #mu_complete_integral[np.invert(case1)] 

    N = np.zeros(bx.shape, dtype=np.uint64)
    N = np.floor((r_integral + mu_initial_integral)/mu_complete_integral)
    integral_remainder = (r_integral - s_mu*mu_initial_integral - (2*np.floor((2*N + 3 - s_mu)/4) - 1)*mu_complete_integral)/(s_mu*(-1)**N)
    J = np.sqrt(M2-M1)*integral_remainder*a
    mu_final = mu_min*CarlsonR.JacobiCN(J, np.sqrt(kSqr))

#    xSqr = (M2 - mu0**2)/(M2-M1)

#    print xSqr

    mu_phi_int1 = BoostRF(1-xSqr, 1-kSqr*xSqr, np.ones(bx.shape))
    mu_phi_int2 = BoostRJ(1-xSqr, 1-kSqr*xSqr, np.ones(bx.shape), 1+n*xSqr)
    elliptic_pi = np.sqrt(xSqr)*mu_phi_int1 + n*np.sqrt(xSqr)**3*mu_phi_int2/3.0
    print elliptic_pi
    g = lambda mu: 1.0/((mu**2 - 1)*math.sqrt((M2[0] - mu**2)*(mu**2 - M1[0])))
    print integrate.quad(g,mu0,mu_final[0])
#(L/(a*np.sqrt(C1)))
    mu_phi_integral = L/np.sqrt(C1)/A/(1-M2)*elliptic_pi
    print mu_phi_integral

    r_phi_integral = PhiTerribleIntegral(r1, r2, r3, r4, a, E, L)
    
    phi = mu_phi_integral + r_phi_integral
#    part2 = -a*math.log((-rplus + r4)/(-rminus + r4))/(2.0*C1)
#    part3 = - a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)    

    return phi, mu_final
