inline void KerrDeflection(double a, double E, double theta, double bx, double by, double &theta_result, double &phi_result){
  const double mu0 = cos(theta),
    C1 = E*E - 1,
    C2 = sqrt(1-a*a),
    mu0Sqr = mu0*mu0,
    aSqr = a*a,
    rplus = 1 + C2,
    rminus = 1 - C2;
  double bxSqr = bx*bx,
    bySqr = by*by,
    L = bx*sqrt(C1)*sin(theta),
    Q = C1*((bxSqr - aSqr)*mu0Sqr + bySqr);

  double s_mu;
  if (mu0 > -1.0 && mu0 < 1.0){
    if (by == 0.0) s_mu = -copysign(1.0, mu0);
    else s_mu = copysign(1.0, by);
  } else {
    phi_result = NAN;
    theta_result = NAN;
    return;
    //    throw("Polar orbits not implemented.");
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
  //
  double r_integral = 4*boost::math::ellint_rf(U12sqr, U13sqr, U14sqr);

  double mu_complete_integral = 2*boost::math::ellint_rf(0.0, (disc - B)/2.0, disc);
  
  double mu_initial_integral;
  if(fabs(M2 - mu0Sqr) > 1e-15){
    mu_initial_integral = boost::math::ellint_rf(mu0Sqr, M2*(mu0Sqr - M1)/(M2-M1), M2)*sqrt(fabs((M2-mu0Sqr)/(M2-M1)))/fabs(a);
  } else {
    mu_initial_integral = mu_complete_integral;
  }

  if (mu0*by < 0.0) mu_initial_integral = mu_complete_integral - mu_initial_integral;

  int N = int((r_integral - mu_initial_integral)/mu_complete_integral);
  
  double integral_remainder = r_integral - N*mu_complete_integral - mu_initial_integral;
  
  double alpha = s_mu*pow(-1.0, N);
  
  double cn, dn;
  boost::math::jacobi_elliptic(k, sqrt(M2-M1)*integral_remainder*a, &cn, &dn);
  double mu_final = sqrt(M2)*cn*alpha;

  theta_result = acos(mu_final);

  //Phi deflection
  if (fabs(bx/by) < 1e-8){
    std::cout << "OH SHIT" << std::endl;
    phi_result = NAN;
    return;
  }


  double xSqr_init = 1 - mu0Sqr/M2,
    xSqr_final = 1- mu_final*mu_final/M2,
    P = 1/sqrt(M2 - M1)/(1-M2);
  
  double pi_complete = 2*boost::math::ellint_3(k, -n),
    pi_init = boost::math::ellint_3(k, -n, asin(sqrt(xSqr_init))),
    pi_final = boost::math::ellint_3(k, -n, asin(sqrt(xSqr_final)));

  if (mu0*s_mu < 0.0) pi_init = pi_complete - pi_init;
  
  if (integral_remainder > mu_complete_integral/2) pi_final = pi_complete - pi_final;

  double mu_phi_integral = (P*(pi_init + pi_final + N*pi_complete)*L/fabs(a) - a*E*r_integral)/sqrt(C1);

  //Evaluate the Terrible Integrals
  // r_-
  double Wsqr = U12sqr - (r3-r1)*(r4-r1)*(rminus-r2)/(rminus-r1),
    Qsqr = (r4 - rminus)/(r4-r1)*Wsqr,
    Psqr = Qsqr + (rminus - r2)*(rminus - r3)*(rminus-r4)/(rminus-r1),
    int2 = (2*(r2 - r1)*(r3 - r1)*(r4 - r1)/3.0/(rminus-r1) * boost::math::ellint_rj(U12sqr, U13sqr, U14sqr, Wsqr) + 2*boost::math::ellint_rc(Psqr, Qsqr) - r_integral/2)/(rminus - r1);

  // r_+
  Wsqr = U12sqr - (r3-r1)*(r4-r1)*(rplus-r2)/(rplus-r1);
  Qsqr = (r4 - rplus)/(r4-r1)*Wsqr;
  Psqr = Qsqr + (rplus - r2)*(rplus - r3)*(rplus-r4)/(rplus-r1);
  double int3 = (2*(r2 - r1)*(r3 - r1)*(r4 - r1)/3.0/(rplus-r1) * boost::math::ellint_rj(U12sqr, U13sqr, U14sqr, Wsqr) + 2*boost::math::ellint_rc(Psqr, Qsqr) - r_integral/2)/(rplus - r1);

  //put all the integrals together
  double r_phi_integral = 2*a*E/sqrt(C1) * (r_integral/2 + (aSqr -a*L/E + rplus*rplus)/(rplus - rminus)*int3  - (aSqr -a*L/E + rminus*rminus)/(rplus - rminus)*int2);
  
  phi_result = mu_phi_integral + r_phi_integral;

}
