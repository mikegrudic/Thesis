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

def roots(e, b):
    if e == "NULL":
        coeffs = (1.0, 0, -b**2, 2*b**2)
    else:
        coeffs = (e**2 - 1, 2.0, -b*b*(e**2 - 1), 2*b*b*(e**2 - 1.0))
    return np.sort(np.roots(coeffs))


def Roots2Newton(e, b):
    if type(b) != np.ndarray:
        b = np.array([b])
    rp= np.zeros(b.shape)
    r1 = np.copy(rp)
    r2 = np.copy(rp)
    code = """
    int i;
    double r0, rnew, val, error1, error2;
    for (i = 0; i < Nb[0]; ++i){
//        std::cout << b[i] << std::endl;
        double c[] = {e*e - 1.0, 2.0, -b[i]*b[i]*(e*e - 1.0), 2*b[i]*b[i]*(e*e - 1.0)};
        r0 = b[i] + e*e/(1.0 - e*e);
        val = c[0]*r0*r0*r0 + c[1]*r0*r0 + c[2]*r0 + c[3];
        error2 = 1e100;
        do {
            error1 = error2;
            rnew = r0 - val/(3.0*c[0]*r0*r0 + 2.0*c[1]*r0 + c[2]);
            val = c[0]*rnew*rnew*rnew + c[1]*rnew*rnew + c[2]*rnew + c[3];
            error2 = fabs(r0 - rnew);
            r0 = rnew;
//            std::cout << fabs(error1 - error2) << '\t';
        } while (fabs(error1-error2) > 1e-14);
        rp[i] = r0;
        double disc = sqrt(4 + 4*pow(b[i],2)*pow(-1 + pow(e,2),2) - (-1 + pow(e,2))*rp[i]*(4 + 3*(-1 + pow(e,2))*rp[i]));
        r1[i] = ((-2.0 - disc)/(e*e - 1) - rp[i])/2.0;
        r2[i] = ((-2.0 + disc)/(e*e - 1) - rp[i])/2.0;
    }
    """
    weave.inline(code, ['e', 'b', 'rp', 'r1', 'r2'])
    return r1, r2, rp

def ComputeDeflections(e, b):
    r1, r2, r3 = Roots2(e,b)
    x = (r3 - r1)*(r3 - r2)
    y = r3*(r3 - r2)
    z = r3*(r3 - r1)
    ellipf = CarlsonR.BoostRF(x, y, z)
    return 4*b*ellipf - np.pi

def newton_deflection(e, b):
    L = b*math.sqrt(e**2 - 1.0)
    newton_e = 0.5 * (e**2 - 1.0)
    return pi - 2*math.acos(1/math.sqrt(2*newton_e*L**2 + 1))

def SolveForAngle(e, theta):
    f = lambda b: deflection(e,b) - theta
    return optimize.brentq(f, 10*bmin(e))

p=OptionParser()
p.add_option("--E", help="Ratio between the total and rest energy of the particle. Let E=NULL for null geodesics.")
(opts,args)=p.parse_args()
    
E = opts.E
if E != "NULL":
    E = float(E)


grid = np.linspace(0, pi, 10000)
b = bmin(E) + 1/np.logspace(-5,10,10000)
#deflections = np.array([deflection(E, B) for B in b])
deflections = ComputeDeflections(E, b)

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
branch_total_sigma = 2*pi*integrate.simps(branch_sigma*np.sin(grid), grid)
sigma_error = np.abs(total_sigma - total_sigma[-1])

np.savetxt("diff_sigma.dat", np.vstack((grid/pi, summed_sigma)).T)
np.savetxt("branch_sigma.dat", np.vstack((grid/pi, branch_sigma)).T)
np.savetxt("sigma.dat", (total_sigma[-1],))
np.savetxt("error.dat", sigma_error.T)
