int nn = Nbx[0];

#pragma omp parallel for
for (int i = 0; i < nn; i++)
{
  KerrDeflection(a, E, theta, bx[i], by[i], theta_result[i], phi_result[i]);
}
