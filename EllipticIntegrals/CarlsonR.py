import numpy as np
from scipy import weave

def RF(x,y,z):
    """ RF
    Calculates Carlson's elliptic integral of the first kind RF(x,y,z)
    
    C++ code comes directly from Numerical Recipes 3e by Press et al., section 6.12
    """

    code = """
    using namespace std;

    const double ERRTOL = 0.0025, THIRD = 1.0/3.0, C1 = 1.0/24.0, C2=0.1, C3 = 3.0/44.0, C4 = 1.0/14.0;
    const double TINY = 5.0*DBL_MIN, BIG= 0.2*DBL_MAX;
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
    return weave.inline(code,['x','y','z'],headers=["<float.h>","<algorithm>"])
