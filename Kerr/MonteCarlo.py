import numpy as np
import KerrDeflection

params = [(0.0001, np.pi/2, 1000000.0, 400, 400,100.0),
          ]

for p in params:
    print p
    counts, data = KerrDeflection.MonteCarloCrossSection(p[0], p[1], p[2], p[3], p[4], p[5])

    np.savetxt("bins_a%g_th%gpi_E%g_bm%g.dat"%(p[0],p[1]/np.pi,p[2],p[5]),data)
    with open("info_a%g_th%gpi_E%g_bm%g.dat"%(p[0],p[1]/np.pi,p[2],p[5]), 'w') as f:
        f.write("""Capture counts: %d
Scattered counts: %d
Estimated capture cross section: %g
"""%(counts, np.sum(data), float(counts)/(np.sum(data)+counts)*np.pi*p[5]**2))
