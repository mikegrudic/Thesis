const double PI = 3.141592653589793;

inline void PNDeriv(double x, double y, double z, double vx, double vy, double vz, double &ax, double &ay, double &az, double eta, int order, int radiation){
  double rSqr = x*x + y*y + z*z;
  double r = sqrt(rSqr);
  double vSqr = vx*vx + vy*vy + vz*vz;
  double v4 = vSqr*vSqr;
  double v6 = v4*vSqr;
  double rdot = (x*vx + y*vy + z*vz)/r;
  double rdotSqr = rdot*rdot;
  double rdot4 = rdotSqr*rdotSqr;
  double rdot6 = pow(rdotSqr,3);

  double A1 = -3.0/2.0*rdotSqr*eta + (1+3*eta)*vSqr -(4+2*eta)/r;
  double A2 = 0.75*(12 + 29*eta)/rSqr + eta*(3-4*eta)*vSqr*vSqr + 15./8*eta*(1-3*eta)*rdotSqr*rdotSqr 
    - 1.5*eta*(3-4*eta)*vSqr*rdotSqr - 0.5*eta*(13 - 4*eta)*vSqr/r - (2 + 25*eta + 2*eta*eta)/r*rdotSqr;
  double A25 = -8.0/5.0*rdot*eta*(3*vSqr/r - 17.0/3.0/rSqr);
  double A3 = 1.0/16.0*eta*rdot6*(-35 + 175*(eta - eta*eta)) + rdot4*vSqr*eta*(15.0/2.0 - 135.0/4.0*eta + 255.0/8.0*eta*eta) + rdotSqr*v4*eta*(-15.0/2.0 + 237.0/8.0*eta - 45.0/2.0*eta*eta) +eta*v6*(11.0/4.0 -49.0/4.0*eta + 13.0*eta*eta) + 1.0/r*(rdot4*eta*(79 - 69.0*eta/2.0 - 30*rdot4*eta*eta) + rdotSqr*vSqr*eta*(-121 + 12*eta + 20*eta*eta) + v4*eta*(75.0/4.0 + 8*eta - 10*eta*eta)) + 1.0/rSqr * (rdotSqr*(1 + 22717.0/168.0*eta + 11.0/8.0*eta*eta - 7.0*eta*eta*eta + 615.0/64.0*eta*PI*PI) + vSqr*(-20827.0/840.0*eta + eta*eta*eta - 123.0/64.0*eta*PI*PI)) + 1.0/rSqr/r*(-16 - 1399.0/12.0*eta - 71.0/2.0*eta*eta + 41.0/16.0*PI*PI*eta);

  double B1 = -2*(2-eta)*rdot;
  double B2 = -0.5*rdot*(eta*(15.0 + 4.0*eta)*vSqr - (4.0 + 41.0*eta + 8.0*eta*eta)/r - 3.0*eta*(3.0+2.0*eta)*rdotSqr);
  double B25 = -8.0/5.0*eta*(6*vSqr - 2/r - 15*rdotSqr);
  double B3 = rdot4*rdot*eta*(-45.0/8.0 + 15*eta +15.0/4.0*eta*eta) + rdotSqr*rdot*vSqr*eta*(12 - 111.0/4.0*eta -12*eta*eta) + rdot*v4*eta*(65.0/8.0 + 19*eta + 6*eta*eta) + 1.0/r*(eta*rdot*rdotSqr*(329.0/6.0 + 59.0/2.0*eta + 18*eta*eta) + rdot*vSqr*eta*(-15 -27*eta - 10*eta*eta)) + 1.0/rSqr*(rdot*(-4 - 5849.0/840.0*eta + 25.0*eta*eta + 8*eta*eta*eta - 123.0/32.0*eta*PI*PI));
  
  double A = 0, B = 0;
  if (order > 0) {
    A += A1;
    B += B1;
  }
  if (order > 1) {
    A += A2;
    B += B2;
    if (radiation) {
      A += A25;
      B += B25;
    }
  }
  if (order > 2) {
    A += A3;
    B += B3;
  }

  ax = -((1+A)*x/r + B*vx)/rSqr;
  ay = -((1+A)*y/r + B*vy)/rSqr;
  az = -((1+A)*z/r + B*vz)/rSqr;
}
