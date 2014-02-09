inline double TerribleIntegral(double r1, double r2, double r3, double r4, double r5){
  double U12sqr = (r4 - r1)*(r4 - r2),
    U13sqr = U12sqr - (r4 - r1)*(r3-r2),
    U14sqr = U12sqr - (r3-r1)*(r4-r2);

  double Wsqr = U12sqr - (r3-r1)*(r4-r1)*(r5-r2)/(r5-r1),
    Qsqr = (r4 - r5)/(r4-r1)*Wsqr,
    Psqr = Qsqr + (r5 - r2)*(r5 - r3)*(r5-r4)/(r5-r1),
    rc = acosh(sqrt(((r1 - r5)*(-r4 + r5))/(r2*r3 - r1*r4 + (r1 - r2 - r3 + r4)*r5)))/sqrt(-((-r2 + r5)*(-r3 + r5)*(-r4 + r5))/(r1 - r5));
    return (2*(r2 - r1)*(r3 - r1)*(r4 - r1)/3.0/(r5-r1) * boost::math::ellint_rj(U12sqr, U13sqr, U14sqr, Wsqr) + 2*rc)/(r5 - r1);
}

inline void KerrDeflection(double a, double E, double theta, double bx, double by, double &theta_result, double &phi_result){
  const double C1 = E*E - 1,
    C2 = sqrt(1-a*a),
    aSqr = a*a,
    rplus = 1 + C2,
    rminus = 1 - C2;

  double mu0 = cos(theta);
  if (mu0 < -1.0) mu0 = -1.0;
  if (mu0 > 1.0) mu0 = 1.0;
  const double mu0Sqr = mu0*mu0;

  const double bxSqr = bx*bx,
    bySqr = by*by,
    L = bx*sqrt(C1)*sin(theta),
    Q = C1*((bxSqr - aSqr)*mu0Sqr + bySqr);

  double s_mu;
  if (-1.0 < mu0 && mu0 < 1.0){
    if (by == 0.0) s_mu = -copysign(1.0, mu0);
    else s_mu = copysign(1.0, by);
  } else {
    s_mu = -copysign(1.0, mu0);
  }

  double rootsr[4], rootsi[4];
  double coeffs[] = {-aSqr*Q, 2*(pow(L - a*E, 2) + Q), aSqr*C1 - L*L - Q, 2.0, C1};

  int info;
  quartic(coeffs, rootsr, rootsi, &info);
    
  std::sort(rootsr, rootsr + 4);

  //Check whether capture trajectory, if so return NaN
  double root_sum = fabs(rootsi[0])+fabs(rootsi[1])+fabs(rootsi[2])+fabs(rootsi[3]);
  if (root_sum != 0.0 || rootsr[3] < rplus){
    phi_result = NAN;
    theta_result = NAN;
    return;
  }
  
  double r1 = rootsr[0],
    r2 = rootsr[1],
    r3 = rootsr[2],
    r4 = rootsr[3];

  //Solve biquadratic polynomial equation M(mu) == 0
  double disc = sqrt(pow(bxSqr + bySqr, 2) + 2*aSqr*(bxSqr - bySqr)*(mu0Sqr - 1.0) + aSqr*aSqr*pow(mu0Sqr - 1.0, 2)),
    A = -aSqr,
    B = aSqr*(mu0Sqr + 1.0) - bxSqr - bySqr,
    C = Q/C1;
  double q = -0.5*(B + copysign(1.0, B)*sqrt(B*B - 4*A*C));
  double M1 = std::min(q/A, C/q),
    M2 = std::max(q/A, C/q),
    k = sqrt(M2/(M2-M1)),
    n = M2/(1-M2);

  //protect against roundoff error
  if (M2 > 1.0) M2 = 1.0;

  double U12sqr = (r4 - r1)*(r4 - r2),
    U13sqr = U12sqr - (r4 - r1)*(r3-r2),
    U14sqr = U12sqr - (r3-r1)*(r4-r2);

  double r_integral = 4*boost::math::ellint_rf(U12sqr, U13sqr, U14sqr);
  double mu_complete_integral = 2*boost::math::ellint_rf(0.0, (disc - B)/2.0, disc);
  double mu_initial_integral;
  //  double mu_initial_integral2;

  if(fabs(M2 - mu0Sqr)/M2 > 1e-15){
    mu_initial_integral = boost::math::ellint_rf(mu0Sqr, M2*(mu0Sqr - M1)/(M2-M1), M2)*sqrt(fabs((M2-mu0Sqr)/(M2-M1)))/fabs(a);
    //    mu_initial_integral2 = mu_complete_integral/2 - mu0*boost::math::ellint_rf(fabs(M1*(mu0Sqr - M2)), fabs(M2*(mu0Sqr - M1)), fabs(-M2*M1))/fabs(a);
  } else {
    mu_initial_integral = mu_complete_integral;
    //mu_initial_integral2 = mu_complete_integral;
  }

  if (mu0*by < 0.0) mu_initial_integral = mu_complete_integral - mu_initial_integral;
  //if (by < 0) mu_initial_integral2 = mu_complete_integral - mu_initial_integral2;

  int N = int((r_integral - mu_initial_integral)/mu_complete_integral);

  double integral_remainder = r_integral - N*mu_complete_integral - mu_initial_integral;

  double alpha = s_mu*pow(-1.0, N);
  
  double cn, dn;
  boost::math::jacobi_elliptic(k, sqrt(M2-M1)*integral_remainder*fabs(a), &cn, &dn);
  double mu_final = sqrt(M2)*cn*alpha;

  theta_result = acos(mu_final);

  double xSqr_init = 1 - mu0Sqr/M2,
    xSqr_final = 1- mu_final*mu_final/M2,
    P = 1/sqrt(M2 - M1)/(1-M2);

  //Insurance against roundoff error: both of these guys should always be non-negative
  xSqr_init = std::max(0.0, xSqr_init);
  xSqr_final = std::max(0.0, xSqr_final);
  
  double pi_complete, pi_init, pi_final;
  if (fabs(bx/by) > 1e-8){
    pi_complete = 2*boost::math::ellint_3(k, -n);
    pi_init = boost::math::ellint_3(k, -n, asin(sqrt(xSqr_init)));
    pi_final = boost::math::ellint_3(k, -n, asin(sqrt(xSqr_final)));
  } else {
    //Zero angular momentum orbits which pass over the poles
    pi_complete = 3.141592653589793;
    pi_init = 0;
    pi_final = 0;
  }

  if (mu0*s_mu < 0.0) pi_init = pi_complete - pi_init;  
  if (integral_remainder > mu_complete_integral/2) pi_final = pi_complete - pi_final;

  double mu_phi_integral;
  if (fabs(bx/by) > 1e-8){
    mu_phi_integral = P*(pi_init + pi_final + N*pi_complete)*L/fabs(a)/sqrt(C1);
  } else {
    mu_phi_integral = pi_init + pi_final + N*pi_complete;
  }

  //Evaluate the Terrible Integrals
  double Tminus = TerribleIntegral(r1, r2, r3, r4, rminus) - r_integral/2/(rminus-r1),
    Tplus = TerribleIntegral(r1, r2, r3, r4, rplus)- r_integral/2/(rplus-r1);

  double r_phi_integral = 2*a*E/sqrt(C1) * ((aSqr -a*L/E + rplus*rplus)/(rplus - rminus)*Tplus  - (aSqr -a*L/E + rminus*rminus)/(rplus - rminus)*Tminus);
  
  phi_result = mu_phi_integral + r_phi_integral;
}
