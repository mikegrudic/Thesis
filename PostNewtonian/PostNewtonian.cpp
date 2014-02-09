#include <iostream>
using namespace std;

inline void PN0(double x, double y, double z, double &ax, double &ay, double &az){
  double r3 = pow(x*x + y*y + z*z, 1.5);
  ax = -x/r3;
  ay = -y/r3;
  az = -z/r3;
}

inline void PN1(double eta, double x, double y, double z, double vx, double vy, double vz, double &ax, double &ay, double &az){
  double r = sqrt(x*x + y*y + z*z);
  double r2 = r*r;
  double rdot = (x*vx + y*vy + z*vz)/r;
  double vSqr = vx*vx + vy*vy + vz*vz;
  double a_n = (1+3*eta)*vSqr - 2*(2+eta)/r - 1.5*eta*rdot*rdot;
  double a_v = -2*(2-eta)*rdot;
  ax = -(a_n*x/r + a_v*vx)/r2;
  ay = -(a_n*y/r + a_v*vy)/r2;
  az = -(a_n*z/r + a_v*vz)/r2;
}
