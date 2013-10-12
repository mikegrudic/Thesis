import numpy as np
import scipy
import math
from scipy import weave, integrate


code = """
    // Initialize the PnVars derivative to zero
//    DTPNVARS = 0.0;

    // Decode the PnVars variable into the basic quantities
//    const double& mA = PNVARS1(0);
//    const double& mB = PNVARS1(1);
    const double& sAx = PNVARS1(0);
    const double& sAy = PNVARS1(1);
    const double& sAz = PNVARS1(2);
    const double& sBx = PNVARS1(3);
    const double& sBy = PNVARS1(4);
    const double& sBz = PNVARS1(5);
    const double& xAx = PNVARS1(6);
    const double& xAy = PNVARS1(7);
    const double& xAz = PNVARS1(8);
    const double& xBx = PNVARS1(9);
    const double& xBy = PNVARS1(10);
    const double& xBz = PNVARS1(11);
    const double& vAx = PNVARS1(12);
    const double& vAy = PNVARS1(13);
    const double& vAz = PNVARS1(14);
    const double& vBx = PNVARS1(15);
    const double& vBy = PNVARS1(16);
    const double& vBz = PNVARS1(17);

    // Compute some convenient quantities
    const double rAB = sqrt((xAx - xBx)*(xAx - xBx)+(xAy - xBy)*(xAy - xBy)+(xAz-xBz)*(xAz-xBz));
    const double vABsq = (vAx-vBx)*(vAx-vBx)+(vAy-vBy)*(vAy-vBy)+(vAz-vBz)*(vAz-vBz);
    const double nABvAB = ((xAx-xBx)*(vAx-vBx) + (xAy-xBy)*(vAy-vBy) +
                           (xAz-xBz)*(vAz-vBz))/rAB;
    const double vAsq = vAx*vAx + vAy*vAy + vAz*vAz;
    const double vBsq = vBx*vBx + vBy*vBy + vBz*vBz;
    const double vAvB = vAx*vBx + vAy*vBy + vAz*vBz;
    const double nABvA = ((xAx-xBx)*vAx + (xAy-xBy)*vAy + (xAz-xBz)*vAz)/rAB;
    const double nABvB = ((xAx-xBx)*vBx + (xAy-xBy)*vBy + (xAz-xBz)*vBz)/rAB;

    const double sAnABvAB = ((sAy*(xAz-xBz)-sAz*(xAy-xBy))*(vAx-vBx)
                            +(sAz*(xAx-xBx)-sAx*(xAz-xBz))*(vAy-vBy)
                            +(sAx*(xAy-xBy)-sAy*(xAx-xBx))*(vAz-vBz))/rAB;
    const double sBnABvAB = ((sBy*(xAz-xBz)-sBz*(xAy-xBy))*(vAx-vBx)
                            +(sBz*(xAx-xBx)-sBx*(xAz-xBz))*(vAy-vBy)
                            +(sBx*(xAy-xBy)-sBy*(xAx-xBx))*(vAz-vBz))/rAB;
    const double nABsA = ((xAx-xBx)*sAx + (xAy-xBy)*sAy + (xAz-xBz)*sAz)/rAB;
    const double nABsB = ((xAx-xBx)*sBx + (xAy-xBy)*sBy + (xAz-xBz)*sBz)/rAB;
//    const double sAsB = sAx*sBx + sAy*sBy + sAz*sBz;

    // Calculate derivatives of the basic quantities
//    double& dtmA = dtPnVars[0];
//    double& dtmB = dtPnVars[1];
    double& dtsAx = DTPNVARS1(0);
    double& dtsAy = DTPNVARS1(1);
    double& dtsAz = DTPNVARS1(2);
    double& dtsBx = DTPNVARS1(3);
    double& dtsBy = DTPNVARS1(4);
    double& dtsBz = DTPNVARS1(5);
    double& dtxAx = DTPNVARS1(6);
    double& dtxAy = DTPNVARS1(7);
    double& dtxAz = DTPNVARS1(8);
    double& dtxBx = DTPNVARS1(9);
    double& dtxBy = DTPNVARS1(10);
    double& dtxBz = DTPNVARS1(11);
    double& dtvAx = DTPNVARS1(12);
    double& dtvAy = DTPNVARS1(13);
    double& dtvAz = DTPNVARS1(14);
    double& dtvBx = DTPNVARS1(15);
    double& dtvBy = DTPNVARS1(16);
    double& dtvBz = DTPNVARS1(17);


    const double OmegaAx = -sBx/(rAB*rAB*rAB)
        + mB/(rAB*rAB*rAB)*((2.*vBy-1.5*vAy)*(xAz-xBz)-(2.*vBz-1.5*vAz)*(xAy-xBy))
        + 3./(rAB*rAB*rAB*rAB)*(nABsB + (mB/mA)*nABsA)*(xAx-xBx);
    const double OmegaAy = -sBy/rAB*rAB*rAB
        + mB/(rAB*rAB*rAB)*((2.*vBz-1.5*vAz)*(xAx-xBx)-(2.*vBx-1.5*vAx)*(xAz-xBz))
        + 3./(rAB*rAB*rAB*rAB)*(nABsB + (mB/mA)*nABsA)*(xAy-xBy);
    const double OmegaAz = -sBz/rAB*rAB*rAB
        + mB/(rAB*rAB*rAB)*((2.*vBx-1.5*vAx)*(xAy-xBy)-(2.*vBy-1.5*vAy)*(xAx-xBx))
        + 3./(rAB*rAB*rAB*rAB)*(nABsB + (mB/mA)*nABsA)*(xAz-xBz);

    const double OmegaBx = -sAx/rAB*rAB*rAB
        + mA/(rAB*rAB*rAB)*((2.*vAy-1.5*vBy)*(xBz-xAz)-(2.*vAz-1.5*vBz)*(xBy-xAy))
        - 3./(rAB*rAB*rAB*rAB)*(nABsA + (mA/mB)*nABsB)*(xBx-xAx);
    const double OmegaBy = -sAy/rAB*rAB*rAB
        + mA/(rAB*rAB*rAB)*((2.*vAz-1.5*vBz)*(xBx-xAx)-(2.*vAx-1.5*vBx)*(xBz-xAz))
        - 3./(rAB*rAB*rAB*rAB)*(nABsA + (mA/mB)*nABsB)*(xBy-xAy);
    const double OmegaBz = -sAz/rAB*rAB*rAB
        + mA/(rAB*rAB*rAB)*((2.*vAx-1.5*vBx)*(xBy-xAy)-(2.*vAy-1.5*vBy)*(xBx-xAx))
        - 3./(rAB*rAB*rAB*rAB)*(nABsA + (mA/mB)*nABsB)*(xBz-xAz);

    // Calculate the spin derivatives
    dtsAx += OmegaAy*sAz - OmegaAz*sAy;
    dtsAy += OmegaAz*sAx - OmegaAx*sAz;
    dtsAz += OmegaAx*sAy - OmegaAy*sAx;

    dtsBx += OmegaBy*sBz - OmegaBz*sBy;
    dtsBy += OmegaBz*sBx - OmegaBx*sBz;
    dtsBz += OmegaBx*sBy - OmegaBy*sBx;

    // Calculate the position derivatives (dt x is just the velocity)
    dtxAx += vAx;
    dtxAy += vAy;
    dtxAz += vAz;

    dtxBx += vBx;
    dtxBy += vBy;
    dtxBz += vBz;

    // Newtonian terms
    const double aA0PN = -mB/(rAB*rAB);
    const double aB0PN = -mA/(rAB*rAB);

//    std::cout << aA0PN << std::endl;
//    std::cout << rAB << std::endl;
//    std::cout << mB << std::endl;

    // 1 PN terms
    const double aA1PNx = (5.*mA*mB)/(rAB*rAB*rAB)
        + 4.*mB*mB/rAB*rAB*rAB
        + mB/rAB*rAB*(1.5*nABvB*nABvB - vAsq + 4*vAvB - 2.*vBsq);
    const double aB1PNx = (5.*mA*mB)/(rAB*rAB*rAB)
        + 4.*mA*mA/rAB*rAB*rAB
        + mA/rAB*rAB*(1.5*nABvA*nABvA - vBsq + 4*vAvB - 2.*vAsq);
    const double aA1PNv = mB/(rAB*rAB)*(4.*nABvA - 3.*nABvB);
    const double aB1PNv = mA/(rAB*rAB)*(3.*nABvA - 4.*nABvB);

    // 2 PN terms
    const double aA2PNx = -(57./4.)*mA*mA*mB/(rAB*rAB*rAB*rAB)
        - (69./2.)*mA*mB*mB/(rAB*rAB*rAB*rAB)
        - 9.*mB*mB*mB/(rAB*rAB*rAB*rAB)
        + mB/rAB*rAB*(-(15./8.)*nABvB*nABvB*nABvB*nABvB + (3./2.)*nABvB*nABvB*vAsq
            - 6.*nABvB*nABvB*vAvB - 2.*vAvB*vAvB
            + (9./2.)*nABvB*nABvB*vBsq + 4.*vAvB*vBsq - 2.*vBsq*vBsq)
        + mA*mB/(rAB*rAB*rAB)*((39./2.)*nABvA*nABvA - 39.*nABvA*nABvB
            + (17./2.)*nABvB*nABvB - (15./4.)*vAsq
            - (5./2.)*vAvB + (5./4.)*vBsq)
        + mB*mB/(rAB*rAB*rAB)*(2.*nABvA*nABvA - 4.*nABvA*nABvB
            - 6.*nABvB*nABvB - 8.*vAvB + 4.*vBsq);
    const double aB2PNx = -(57./4.)*mB*mB*mA/(rAB*rAB*rAB*rAB)
        - (69./2.)*mB*mA*mA/(rAB*rAB*rAB*rAB)
        - 9.*mA*mA*mA/(rAB*rAB*rAB*rAB)
        + mA/rAB*rAB*(-(15./8.)*nABvA*nABvA*nABvA*nABvA + (3./2.)*nABvA*nABvA*vBsq
            - 6.*nABvA*nABvA*vAvB - 2.*vAvB*vAvB
            + (9./2.)*nABvA*nABvA*vAsq + 4.*vAvB*vAsq - 2.*vAsq*vAsq)
        + mA*mB/(rAB*rAB*rAB)*((39./2.)*nABvB*nABvB - 39.*nABvB*nABvA
            + (17./2.)*nABvA*nABvA - (15./4.)*vBsq
            - (5./2.)*vAvB + (5./4.)*vAsq)
        + mA*mA/(rAB*rAB*rAB)*(2.*nABvB*nABvB - 4.*nABvB*nABvA
            - 6.*nABvA*nABvA - 8.*vAvB + 4*vAsq);

    const double aA2PNv = mB*mB/(rAB*rAB*rAB)*(-2.*nABvA - 2.*nABvB)
        + mA*mB/(rAB*rAB*rAB)*(-(63./4.)*nABvA + (55./4.)*nABvB)
        + mB/rAB*rAB*(-6.*nABvA*nABvB*nABvB + (9./2.)*nABvB*nABvB*nABvB + nABvB*vAsq
            - 4.*nABvA*vAvB + 4.*nABvB*vAvB + 4.*nABvA*vBsq
            - 5.*nABvB*vBsq);
    const double aB2PNv = mA*mA/(rAB*rAB*rAB)*(2.*nABvB + 2.*nABvA)
        + mA*mB/(rAB*rAB*rAB)*((63./4.)*nABvB - (55./4.)*nABvA)
        + mA/rAB*rAB*(6.*nABvB*nABvA*nABvA - (9./2.)*nABvA*nABvA*nABvA - nABvA*vBsq
            + 4.*nABvB*vAvB - 4.*nABvA*vAvB - 4.*nABvB*vAsq
            + 5.*nABvA*vAsq);

    // spin-orbit terms (enter between 1 and 2 PN)
    const double aASOx = 6.0/(rAB*rAB*rAB)*(sBnABvAB + (mB/mA)*sAnABvAB);
    const double aBSOx = 6.0/(rAB*rAB*rAB)*(sAnABvAB + (mA/mB)*sBnABvAB);

    // leading order RR terms (enter at 2.5 PN)
    const double aA0RRx = (208.*mA*mB*mB)/(15.*rAB*rAB*rAB*rAB)*nABvAB
                        - (24.*mA*mA*mB)/(5.*rAB*rAB*rAB*rAB)*nABvAB
                        + (12.*mA*mB)/(5.*rAB*rAB*rAB)*nABvAB*vABsq;
    const double aB0RRx = (208.*mB*mA*mA)/(15.*rAB*rAB*rAB*rAB)*nABvAB
                        - (24.*mB*mB*mA)/(5.*rAB*rAB*rAB*rAB)*nABvAB
                        + (12.*mA*mB)/(5.*rAB*rAB*rAB)*nABvAB*vABsq;
    const double aA0RRv = (8.*mA*mA*mB)/(5.*rAB*rAB*rAB*rAB)
                        - (32.*mA*mB*mB)/(5.*rAB*rAB*rAB*rAB)
                        - (4.*mA*mB)/(5.*rAB*rAB*rAB)*vABsq;
    const double aB0RRv = (8.*mB*mB*mA)/(5.*rAB*rAB*rAB*rAB)
                        - (32.*mB*mA*mA)/(5.*rAB*rAB*rAB*rAB)
                        - (4.*mA*mB)/(5.*rAB*rAB*rAB)*vABsq;

    // Initialize the meta-terms in the acceleration equations
    double aAx(0),aAv(0),aBx(0),aBv(0);     // scalar terms
    double cSAx(0),cSAy(0),cSAz(0);         // cross-product spin-orbit terms
    double cSBx(0),cSBy(0),cSBz(0);
    // Add the integer PostNewtonian terms (to the PN order specified)
    if (mPnOrder == 0) {
        aAx = aA0PN;
        aBx = aB0PN;
    } else if (mPnOrder == 1) {
        aAx = aA0PN + aA1PNx;
        aAv = aA1PNv;
        aBx = aB0PN + aB1PNx;
        aBv = aB1PNv;
    } else if (mPnOrder == 2) {
        aAx = aA0PN + aA1PNx + aA2PNx;
        aAv = aA1PNv + aA2PNv;
        aBx = aB0PN + aB1PNx + aB2PNx;
        aBv = aB1PNv + aB2PNv;
    } else {
        std::cout << "PnOrder = "<<mPnOrder<< " is not implemented." << std::endl;
//        return -1;
    }
    // Add the radiation reaction terms (2.5PN) unless conservative requested
    if (mRadiationReaction) {
        aAx += aA0RRx;
        aAv += aA0RRv;
        aBx += aB0RRx;
        aBv += aB0RRv;
    }
    // Add the Spin-Orbit terms
    if (mSpinOrbit) {
        // Scalar terms
        aAx += aASOx;
        aBx += aBSOx;
        // Cross-product terms
        cSAx =  4./(rAB*rAB*rAB)*(sBy*(vAz-vBz)-sBz*(vAy-vBy))
              + 3./(rAB*rAB*rAB)*(sAy*(vAz-vBz)-sAz*(vAy-vBy))*(mB/mA)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sBy*(xAz-xBz)-sBz*(xAy-xBy))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sAy*(xAz-xBz)-sAz*(xAy-xBy))*(mB/mA);
        cSAy =  4./(rAB*rAB*rAB)*(sBz*(vAx-vBx)-sBx*(vAz-vBz))
              + 3./(rAB*rAB*rAB)*(sAz*(vAx-vBx)-sAx*(vAz-vBz))*(mB/mA)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sBz*(xAx-xBx)-sBx*(xAz-xBz))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sAz*(xAx-xBx)-sAx*(xAz-xBz))*(mB/mA);
        cSAz =  4./(rAB*rAB*rAB)*(sBx*(vAy-vBy)-sBy*(vAx-vBx))
              + 3./(rAB*rAB*rAB)*(sAx*(vAy-vBy)-sAy*(vAx-vBx))*(mB/mA)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sBx*(xAy-xBy)-sBy*(xAx-xBx))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sAx*(xAy-xBy)-sAy*(xAx-xBx))*(mB/mA);
        cSBx =  4./(rAB*rAB*rAB)*(sAy*(vBz-vAz)-sAz*(vBy-vAy))
              + 3./(rAB*rAB*rAB)*(sBy*(vBz-vAz)-sBz*(vBy-vAy))*(mA/mB)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sAy*(xBz-xAz)-sAz*(xBy-xAy))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sBy*(xBz-xAz)-sBz*(xBy-xAy))*(mA/mB);
        cSBy =  4./(rAB*rAB*rAB)*(sAz*(vBx-vAx)-sAx*(vBz-vAz))
              + 3./(rAB*rAB*rAB)*(sBz*(vBx-vAx)-sBx*(vBz-vAz))*(mA/mB)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sAz*(xBx-xAx)-sAx*(xBz-xAz))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sBz*(xBx-xAx)-sBx*(xBz-xAz))*(mA/mB);
        cSBz =  4./(rAB*rAB*rAB)*(sAx*(vBy-vAy)-sAy*(vBx-vAx))
              + 3./(rAB*rAB*rAB)*(sBx*(vBy-vAy)-sBy*(vBx-vAx))*(mA/mB)
              - 6./(rAB*rAB*rAB*rAB)*nABvAB*(sAx*(xBy-xAy)-sAy*(xBx-xAx))
              - 3./(rAB*rAB*rAB*rAB)*nABvAB*(sBx*(xBy-xAy)-sBy*(xBx-xAx))*(mA/mB);
    }
    
    //-------------------------------------------------------
    // Calculate accelerations (for any implemented PN order)
    //-------------------------------------------------------
    dtvAx += aAx*(xAx-xBx)/rAB
            +aAv*(vAx-vBx) + cSAx;
    dtvAy += aAx*(xAy-xBy)/rAB
            +aAv*(vAy-vBy) + cSAy;
    dtvAz += aAx*(xAz-xBz)/rAB
            +aAv*(vAz-vBz) + cSAz;
            
    dtvBx += aBx*(xBx-xAx)/rAB
            +aBv*(vBx-vAx) + cSBx;
    dtvBy += aBx*(xBy-xAy)/rAB
            +aBv*(vBy-vAy) + cSBy;
    dtvBz += aBx*(xBz-xAz)/rAB
            +aBv*(vBz-vAz) + cSBz;
"""

def PostNewtonianAccel(PnVars, t, mA=1, mB=1, mPnOrder=2, mSpinOrbit=1, mRadiationReaction=0):
    dtPnVars = np.zeros(18)
    weave.inline(code, ['mA', 'mB', 'PnVars', 'dtPnVars', 'mPnOrder', 'mSpinOrbit', 'mRadiationReaction'])
    return dtPnVars

a = 1000

x = np.array([0,0,0,0,0,0,-a,0,0,a,0,0,0,math.sqrt(0.25/a),0,0,-math.sqrt(0.25/a),0])
print PostNewtonianAccel(x, 0)
T = 2*np.pi*math.sqrt(4*a**3)

t= np.linspace(0, T, 1000)

xx = integrate.odeint(PostNewtonianAccel,x, t)
print xx[:5]
np.savetxt("traj.dat", xx[:,6:12])

#print PostNewtonianAccel(x, 0, mPnOrder=0)
