import numpy as np
import math
from scipy import weave, integrate, linalg
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

def bmin(e):
    if e == "NULL":
        return 3*math.sqrt(3)
    else:
        return math.sqrt((8 - 36*e**2 + 27*e**4 + e*(9*e**2 - 8)**(3.0/2.0))/2)/(e**2 - 1)

def PhiTerribleIntegral(r1, r2, r3, r4, a, E, L):
    C1 = math.sqrt(1-a**2)
    rplus, rminus = 1 + C1, 1 - C1
    int1 = InvSqrtQuartic(r1, r2, r3, r4, r4)
    int2 = TerribleIntegral(r1, r2, r3, r4, rplus, r4)
    int3 = TerribleIntegral(r1, r2, r3, r4, rminus, r4)
    result = 2*a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
    return result

def SchwarzDeflection(E, b):
    bm = bmin(E)
    result = np.zeros(b.shape)
    fall_in= (b<=bm)
    result[fall_in] = np.NaN
    r1, r2, r3 = CubicRoots(E,b[b>bmin(E)])
    x = (r3 - r1)*(r3 - r2)
    y = r3*(r3 - r2)
    z = r3*(r3 - r1)
    ellipf = CarlsonR.BoostRF(x, y, z)
    result[np.invert(fall_in)] = 4*b[np.invert(fall_in)]*ellipf
    return result

def EquatorialDeflection(a, E, b, roots):
    C1 = math.sqrt(1-a**2)
    L = b*math.sqrt(E**2-1)
    rplus, rminus = 1 + C1, 1 - C1
    r1, r2, r3, r4 = roots
    int1 = InvSqrtQuartic(r1, r2, r3, r4, r4)
    int2 = TerribleIntegral(r1, r2, r3, r4, rplus, r4)
    int3 = TerribleIntegral(r1, r2, r3, r4, rminus, r4)
    part1 = (L - a*E)/math.sqrt(E**2 - 1) * int1
    part3 = a*E/math.sqrt(E**2 - 1) * (int1 + (a**2 -a*L/E + rplus**2)/(rplus - rminus)*int2  - (a**2 -a*L/E + rminus**2)/(rplus - rminus)*int3)
    phi_result = 2*(part1 + part3)
    return phi_result%(2*np.pi)

def KerrDeflection(a, theta, E, bx, by):
    phi_result = np.empty(bx.shape)
    mu_result = np.empty(bx.shape)
    
    mu0 = math.cos(theta)
    C1 = E**2 - 1
    C2 = math.sqrt(1-a**2)
    L = bx*math.sqrt(C1)*math.sin(theta)
    Q = (E**2 - 1)*((bx**2 - a**2)*mu0**2 + by**2)
    ones = np.ones(bx.shape)
    zeros = np.zeros(bx.shape)

    #Schwarzschild case
    if a==0.0:
        sch_def = SchwarzDeflection(E, np.sqrt(bx**2 + by**2))
        t1 = np.arctan2(by,bx)
        y = np.sin(t1)*np.sin(sch_def)
        x = np.cos(sch_def)*math.sin(theta) - np.sin(sch_def)*math.sin(theta)*np.sin(t1)
        phi_result = np.arctan2(y, x)
        theta_result = np.arccos(np.cos(sch_def)*math.cos(theta)+math.sin(theta)*np.sin(sch_def)*np.sin(t1))
        return phi_result%(2*np.pi), theta_result%(np.pi)

    s_mu = np.zeros(bx.shape, dtype=np.int64)
    if -1 < mu0 < 1:
        q1 = bx**2*math.cos(theta)
        s_mu[by > 0] = 1
        s_mu[by < 0] = -1
        s_mu[(by == 0)*(math.cos(theta)>0)] = -1
        s_mu[(by == 0)*(math.cos(theta)<0)] = 1
    elif mu0 == 1:
        s_mu = -np.ones(by.shape).astype(int)
    else:
        s_mu = np.ones(by.shape).astype(int)

    #solve r quartic
    r_coeffs = np.array([C1*ones, 2*ones, a**2*C1 - L**2 - Q, 2*((-(a*E) + L)**2 + Q), -a**2*Q])

    r_coeffs = r_coeffs/C1
    companion = np.swapaxes(np.array(
        [[zeros, zeros, zeros, -r_coeffs[4]],
        [ones, zeros, zeros, -r_coeffs[3]],
        [zeros, ones, zeros, -r_coeffs[2]],
        [zeros, zeros, ones, -r_coeffs[1]]]),0,2)
    r_roots = np.sort(np.linalg.eigvals(companion),axis=1)

    #Return NaN deflection for trajectories which go into the hole
    falls_in = (np.sum(np.abs(r_roots.imag),axis=1) > 0.0) + (np.max(r_roots.real,axis=1) < 1+C2)
    doesnt_fall_in = np.invert(falls_in)
    phi_result[falls_in] = np.NaN
    mu_result[falls_in] = np.NaN
    
    L, bx, by, r_roots, zeros, ones, s_mu = L[doesnt_fall_in], bx[doesnt_fall_in], by[doesnt_fall_in], r_roots[doesnt_fall_in].real, zeros[doesnt_fall_in], ones[doesnt_fall_in], s_mu[doesnt_fall_in]

#    print r_roots
    r1, r2, r3, r4 = r_roots.T
    
    #Equatorial case is easy
    if mu0-math.cos(np.pi/2) == 0.0 and np.all(by==0.0):
        phi_result[doesnt_fall_in] = EquatorialDeflection(a, E, bx, (r1,r2,r3,r4))
        return phi_result%(2*np.pi), ones*np.pi/2
    
    #mu biquadratic
    discriminant = np.sqrt((bx**2 + by**2)**2 + 2*a**2*(bx - by)*(bx + by)*(-1 + mu0**2) + a**4*(-1 + mu0**2)**2)
    A = -a**2
    B = a**2 * (1 + mu0**2) - bx**2 - by**2
    C = by**2 + (bx**2 - a**2)*mu0**2
    q = -0.5*(B + np.sign(B)*np.sqrt(B**2 - 4*A*C))
    M1, M2 = np.sort((q/A, C/q),axis = 0)
    aSqrM2 = (-bx**2 - by**2 + discriminant + a**2*(1+mu0**2))/2
    aSqrM1 = (-bx**2 - by**2 - discriminant + a**2*(1+mu0**2))/2
    mu_max = np.sqrt(M2)
    mu_min = -np.sqrt(M2)
    kSqr = M2/(M2 - M1)
    xSqr = np.abs(1 - mu0**2/M2)
    n = M2/(1-M2)

#do integrals
    r_integral = 2*InvSqrtQuartic(r1, r2, r3, r4, r4)

    mu_complete_integral = 2*CarlsonR.BoostRF(zeros, (bx**2+by**2+discriminant - a**2*(1+mu0**2))/2.0, discriminant)
    
    case1 = np.abs(np.abs(mu0)-mu_max) > 1e-16
    mu_initial_integral = np.empty(mu_complete_integral.shape)
    mu_initial_integral[case1] = (mu_complete_integral/2 - 1/np.sqrt(-M1*M2)*math.fabs(mu0)*CarlsonR.RF((M2-mu0**2)/M2, (M1-mu0**2)/M1, ones)/a)[case1]

    mu_initial_integral[np.invert(case1)] = mu_complete_integral[np.invert(case1)]

    A = np.sign(by)*np.sign(mu0) == -1    
    mu_initial_integral[A] = mu_complete_integral[A] - mu_initial_integral[A]

    N = np.floor((r_integral - mu_initial_integral)/mu_complete_integral)

    integral_remainder = r_integral - N*mu_complete_integral - mu_initial_integral
    
    alpha = s_mu*(-1)**N

    J = np.sqrt(M2-M1)*integral_remainder*a
    mu_final = mu_max*CarlsonR.JacobiCN(J, np.sqrt(kSqr))*alpha

# Do mu-integrals for phi deflection
    xSqr_init = 1 - mu0**2/M2
    xSqr_final = 1 - mu_final**2/M2

    P = 1/np.sqrt(M2 - M1)/(1-M2)

    pi_complete = P*2*CarlsonR.LegendrePiComplete(-n, kSqr)
    pi_init = P*CarlsonR.LegendrePi(-n, xSqr_init, kSqr)
    pi_final = P*CarlsonR.LegendrePi(-n, xSqr_final, kSqr)
    
    if mu0>0:
        pi_init[s_mu==-1] = pi_complete[s_mu==-1] - pi_init[s_mu==-1]
    else:
        pi_init[s_mu==1] = pi_complete[s_mu==1] - pi_init[s_mu==1]        

    A = integral_remainder > mu_complete_integral/2
    pi_final[A] = pi_complete[A] - pi_final[A]
    
    mu_phi_integral = ((pi_init + pi_final + N*pi_complete)*L/a - a*E*r_integral)/np.sqrt(C1)
    
    r_phi_integral = PhiTerribleIntegral(r1, r2, r3, r4, a, E, L)

    phi = mu_phi_integral + r_phi_integral
    
    phi_result[doesnt_fall_in] = phi
    mu_result[doesnt_fall_in] = mu_final

    return phi_result%(2*np.pi), np.arccos(mu_result)%(np.pi)
