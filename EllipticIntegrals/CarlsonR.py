import numpy as np
import math
from scipy import weave

RFcode = """
    using namespace std;

    static const double ERRTOL = 0.0025, THIRD = 1.0/3.0, C1 = 1.0/24.0, C2=0.1, C3 = 3.0/44.0, C4 = 1.0/14.0;
    static const double TINY = 5.0*DBL_MIN, BIG= 0.2*DBL_MAX;
    double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;
    if (min(min(x,y),z) < 0.0 || min(min(x+y,x+z),y+z) < TINY || max(max(x,y),z) > BIG) throw("Invalid arguments in rf");
    xt = x;
    yt = y;
    zt = z;
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx*(sqrty + sqrtz) + sqrty*sqrtz;
        xt = 0.25*(xt + alamb);
        yt = 0.25*(yt + alamb);
        zt = 0.25*(zt + alamb);
        ave = THIRD*(xt + yt + zt);
        delx = (ave - xt)/ave;
        delz = (ave - zt)/ave;
        dely = (ave - yt)/ave;
    } while (max(max(abs(delx), abs(dely)),abs(delz)) > ERRTOL);
    e2 = delx*dely - delz*delz;
    e3 = delx*dely*delz;


    return_val = (1.0 + (C1*e2 - C2 - C3*e3)*e2 + C4*e3)/sqrt(ave);
    """

RJcode = """
using namespace std;

static const double ERRTOL = 0.0015, C1 = 3.0/4.0, C2 = 1.0/3.0, C3 = 3.0/22.0, C4 = 3.0/26.0, C5 = 0.75*C3, C6 = 1.5*C4, C7 = 0.5*C2, C8 = C3 + C3;
static const double TINY = pow(5.0*DBL_MIN, 1./3.) , BIG= 0.3*pow(0.2*DBL_MAX, 1./3.);

double a, alamb, alpha, ans, ave, b, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt, rcx, rho, sqrtx, sqrty, sum, tau, xt, yt, zt;

if (min(min(x,y),z) < 0.0 || min(min(x+y,x+z),min(y+z,abs(p))) < TINY || max(max(x,y),max(z,abs(p))) > BIG) throw("Invalid arguments in RJ");



"""


def RF(x,y,z):
    """ RF
    Calculates Carlson's elliptic integral of the first kind RF(x,y,z)
    
    C++ code comes directly from Numerical Recipes 3e by Press et al., section 6.12
    """

    return weave.inline(RFcode,['x','y','z'],headers=["<float.h>","<algorithm>"])

def BoostRF(x,y,z):
    RFcode = """
    return_val = boost::math::ellint_rf(x, y, z);
    """
    return weave.inline(RFcode, ['x','y','z'], headers = ["<boost/math/special_functions/ellint_rf.hpp>"])

def BoostRJ(x,y,z,p):
    RJcode = """
    return_val = boost::math::ellint_rj(x, y, z, p);
    """
    return weave.inline(RJcode, ['x','y','z','p'], headers = ["<boost/math/special_functions/ellint_rj.hpp>"])

def BoostRC(x,y):
    RCcode = """
    return_val = boost::math::ellint_rc(x, y);
    """
    return weave.inline(RCcode, ['x','y'], headers = ["<boost/math/special_functions/ellint_rc.hpp>"])

def InvSqrtQuartic(r1, r2, r3, r4, a, b=None):
    if b == None:
        U12 = math.sqrt((a - r1)*(a - r2)) + math.sqrt((a -r3)*(a -r4))
    else:
        U12 = (math.sqrt((b-r1)*(b-r2)*(a-r3)*(a-r4)) + math.sqrt((b-r4)*(b-r3)*(a-r2)*(a-r1)))/(b-a)

    U12squared = U12**2
    return 2*RF(U12squared, U12squared - (r4 - r1)*(r3 - r2), U12squared - (r3 - r1)*(r4 - r2))

def TerribleIntegral(r1,r2,r3,r4,r5, a, b=None):
    if b == None:
        U12 = math.sqrt((a - r1)*(a - r2)) + math.sqrt((a - r3)*(a - r4))
        U12sqr = U12*U12
        Wsqr = U12sqr - (r3 - r1)*(r4 - r1)*(r5 - r2)/(r5 - r1)
        Qsqr = (a - r5)/(a - r1)*Wsqr
    else:
        U12 = (math.sqrt((b-r1)*(b-r2)*(a-r3)*(a-r4)) + math.sqrt((b-r4)*(b-r3)*(a-r2)*(a-r1)))/(b-a)
        U12sqr = U12*U12
        Wsqr = U12sqr - (r3 - r1)*(r4 - r1)*(r5 - r2)/(r5 - r1)
        Qsqr = ((a - r5)*(b - r5)/(a - r1)/(b - r1))*Wsqr

    U13sqr = U12sqr - (r4 - r1)*(r3 - r2)
    U14sqr = U12sqr - (r3 - r1)*(r4 - r2)
    Psqr = Qsqr + (r5 - r2)*(r5 - r3)*(r5 - r4)/(r5 - r1)
    I3 = 2*(r2 - r1)*(r3 - r1)*(r4 - r1)/3.0/(r5-r1) * BoostRJ(U12sqr, U13sqr, U14sqr, Wsqr) + 2*BoostRC(Psqr, Qsqr)
    I1 = 2*BoostRF(U12sqr, U13sqr, U14sqr)

    return (I3 - I1)/(r5-r1)
