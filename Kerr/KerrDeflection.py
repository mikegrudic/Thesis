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
#        raise Exception("Particle falls into the black hole!")
        return np.zeros(L.shape)
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
    part3 = -a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
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
    L = bx*math.sqrt(C1)*math.sin(theta)
    Q = (E**2 - 1)*((bx**2 - a**2)*mu0**2 + by**2)
    ones = np.ones(bx.shape)
    zeros = np.zeros(bx.shape)
#    if a==0:
#        deflection = SchwarzDeflection(E, np.sqrt(bx**2 + by**2))        
    mu0 = np.cos(theta)

    if mu0 == 0:
        return EqDeflection(a, E, L), np.pi/2
    s_mu = np.zeros(bx.shape, dtype=np.int64)
    if theta != 0 and theta != pi:
        s_mu = np.sign(-(2*bx*mu0+(-1 + 3*np.cos(2*theta))*by)*np.sin(theta))
    elif theta == 0:
        s_mu = -np.ones(by.shape)
    else:
        s_mu = np.ones(by.shape)
    
    C1 = E**2 - 1
    C2 = math.sqrt(1-a**2)

#Solve r and mu polynomials    
    r_coeffs = np.array([C1*ones, 2*ones, a**2*C1 - L**2 - Q, 2*((-(a*E) + L)**2 + Q), -a**2*Q]).T
    mu_coeffs = np.array([C1*a**2*ones, C1*a**2 - L**2 - Q, Q]).T
# r quartic
    r_roots = np.sort(np.array([np.roots(c) for c in r_coeffs]), axis=1)
#mu biquadratic
    discriminant = np.sqrt((bx**2 + by**2)**2 + 2*a**2*(bx - by)*(bx + by)*(-1 + mu0**2) + a**4*(-1 + mu0**2)**2)
    A = -a**2
    B = a**2 * (1 + mu0**2) - bx**2 - by**2
    C = by**2 + (bx**2 - a**2)*mu0**2
    q = -0.5*(B + np.sign(B)*np.sqrt(B**2 - 4*A*C))
    M1, M2 = np.sort((q/A, C/q),axis = 0)
    aSqrM2 = (-bx**2 - by**2 + discriminant + a**2*(1+mu0**2))/2
    mu_max = np.sqrt(M2)
    mu_min = -np.sqrt(M2)
    kSqr = M2/(M2 - M1)
    xSqr = np.abs(1 - mu0**2/M2)
    n = M2/(1-M2)
    A = a*np.sqrt(M2 - M1)

#do complete integrals
    r1, r2, r3, r4 = r_roots[:,0], r_roots[:,1], r_roots[:,2], r_roots[:,3]
    r_integral = 2*InvSqrtQuartic(r1, r2, r3, r4, r4)
#    mu_complete_integral = 2*CarlsonR.BoostRF(zeros, 1.0-kSqr, ones)/A
    mu_complete_integral = 2*CarlsonR.BoostRF(zeros, (bx**2+by**2+discriminant - a**2*(1+mu0**2))/2.0, discriminant)
#    print mu_complete_integral

#do partial mu integral from mu0 to mu_max    
    case1 = (mu0 != mu_max)*(mu0 != mu_min)
    mu_initial_integral = np.empty(bx.shape)
    if mu0 > 0:
        U1 = (np.abs(-mu0**4 + (M1[case1]+M2[case1])*mu0**2 + -M1[case1]*M2[case1] + (mu0**2 - mu_max[case1]**2)**2))/(mu_max[case1]-mu0)**2
#        U1 = (-a**4*mu0**4 + (aSqrM2 + by**2(mu0**2 - 1) - a**2*mu0**2*(mu0**2 + 1))**2)/(aSqrM2 - a**2*mu0**2)
        mu_initial_integral[case1] = 2*CarlsonR.BoostRF(U1, U1 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U1 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1]))/a
    else:
    #break into mu0 to 0 + 0 to mu_max
        U2 = (-M1[case1]*M2[case1] + mu_max[case1]**4)/mu_max[case1]**2
        U3 = (np.sqrt(np.abs(-mu0**4 + (M1[case1]+M2[case1])*mu0**2 - M1[case1]*M2[case1]))+np.sqrt(-M1[case1]*M2[case1]))**2/mu0**2 #unstable, fix this
        mu_initial_integral[case1] = 2*(CarlsonR.BoostRF(U2, U2 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U2 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1])) + CarlsonR.BoostRF(U3, U3 - (M1[case1]+M2[case1]) + 2*np.sqrt(-M1[case1]*M2[case1]), U3 - (M1[case1]+M2[case1]) - 2*np.sqrt(-M1[case1]*M2[case1])))/a

    mu_initial_integral[np.invert(case1)] = 0.0 #mu_complete_integral[np.invert(case1)] 

    N = np.zeros(bx.shape, dtype=np.uint64)
    N = np.floor((r_integral + mu_initial_integral)/mu_complete_integral)
    alpha1 = s_mu
    alpha2 = s_mu*(-1)**N
    alpha3 = 2*np.floor((2*N + 3 - s_mu)/4) - 1

    integral_remainder = (r_integral - alpha1*mu_initial_integral - alpha3*mu_complete_integral)/alpha2
    J = np.sqrt(M2-M1)*integral_remainder*a
    mu_final = mu_min*CarlsonR.JacobiCN(J, np.sqrt(kSqr))
# Do mu-integrals for phi deflection
    xSqr_init = np.abs(1 - mu0**2/M2)
    xSqr_final = np.abs(1 - mu_final**2/M2)
  
    pi_complete = 2*CarlsonR.LegendrePiComplete(-n, kSqr)
    pi_init = CarlsonR.LegendrePi(-n, xSqr_init, kSqr)
    pi_final = CarlsonR.LegendrePi(-n, xSqr_final, kSqr)

    if mu0 < 0.0:
        pi_init = pi_complete - pi_init
        pi_final[mu_final < 0] = pi_complete[mu_final<0] - pi_final[mu_final<0]
    else:
        pi_final[mu_final > 0] = pi_complete[mu_final>0] - pi_final[mu_final>0]

    P = 1/np.sqrt(M2 - M1)/(1-M2)        
    mu_phi_integral = (L*P/a*(alpha1*pi_init + alpha2*pi_final + alpha3*pi_complete) - a*E*r_integral)/np.sqrt(E**2 - 1)
    r_phi_integral = PhiTerribleIntegral(r1, r2, r3, r4, a, E, L)
    phi = mu_phi_integral + r_phi_integral

    return phi - pi, np.arccos(mu_final)
