import numpy as np
import math
from optparse import OptionParser
import CarlsonR
import scipy
from scipy import special, integrate, misc, optimize, interpolate, weave

pi=np.pi

def bmin(e):
    if e == "NULL":
        return 3*math.sqrt(3)
    else:
        return math.sqrt((8 - 36*e**2 + 27*e**4 + e*(9*e**2 - 8)**(3.0/2.0))/2)/(e**2 - 1)

def Roots2(e, b):
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

def ComputeDeflections(e, b):
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

def Sigma(e):
    deflections, diff_sigma = DiffSigma(e)
    deflections_modpi = deflections % pi
    max_branch = np.max(deflections) - np.max(deflections) % pi
    num_branches = int(max_branch/pi + 0.5)
    branch_bounds = np.linspace(pi, max_branch, num_branches)

    branch_index = [(deflections > bound - pi)*(deflections < bound) for bound in branch_bounds]
    branch_total_sigma = np.array([2*pi*integrate.simps(diff_sigma[i]*np.sin(deflections_modpi[i]), deflections_modpi[i]) for i in branch_index])
    print branch_total_sigma
    return np.sum(branch_total_sigma)

def DiffSigma(e):
    b = bmin(e) + 1/np.logspace(-10,2,10000)
    deflections = ComputeDeflections(e, b)
#    deflections = newton_deflection(e, b)
    ddef_db = np.gradient(deflections)/np.gradient(b)
    return deflections, np.abs(b/np.sin(deflections)) * np.abs(1.0/ddef_db)

def newton_deflection(e, b):
    L = b*math.sqrt(e**2 - 1.0)
    newton_e = 0.5 * (e**2 - 1.0)
    return pi - 2*np.arccos(1/np.sqrt(2*newton_e*L**2 + 1))

def SolveForAngle(e, theta):
    f = lambda b: deflection(e,b) - theta
    return optimize.brentq(f, 10*bmin(e))

#p=OptionParser()
#p.add_option("--E", help="Ratio between the total and rest energy of the particle. Let E=NULL for null geodesics.")
#(opts,args)=p.parse_args()
    
#E = opts.E
#if E != "NULL":
#    E = float(E)

E = 1+np.logspace(-5, 2, 1000)
#sigma = np.array([Sigma(e) for e in E])

#print DiffSigma(1.1)

#print E.shape, sigma.shape
#np.savetxt("sigma.dat", np.vstack((E, sigma)).T)

#branch_sigma = np.array([np.interp(grid, deflections_modpi[i], sigma[i]) for i in branch_index])

#summed_sigma = np.cumsum(branch_sigma, axis=0)

#total_sigma = 2*pi*integrate.simps(summed_sigma*np.sin(grid), grid)
#branch_total_sigma = 2*pi*integrate.simps(branch_sigma*np.sin(grid), grid)
#sigma_error = np.abs(total_sigma - total_sigma[-1])

#np.savetxt("diff_sigma.dat", np.vstack((grid/pi, summed_sigma)).T)
#np.savetxt("branch_sigma.dat", np.vstack((grid/pi, branch_sigma)).T)
#np.savetxt("sigma.dat", (total_sigma[-1],))
#np.savetxt("error.dat", sigma_error.T)
