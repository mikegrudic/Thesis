import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import animation
import sys

a = 1.0

r, phi, theta = np.loadtxt(sys.argv[2]).T
bound = float(sys.argv[1])

#r, phi, theta = r[rbound], phi[rbound], theta[rbound]
x = r*np.sin(phi)*np.sin(theta)
y = r*np.cos(phi)*np.sin(theta)
z = r*np.cos(theta)

rbound = np.max(np.abs((x,y,z)),axis=0) < bound

x, y, z = x[rbound], y[rbound], z[rbound]

fig = plt.figure()
ax = fig.gca(projection='3d',aspect='equal')
ax._axis3don = True

ax.plot(x,y,z)

for direction in (-1, 1):
    for point in np.diag(direction * bound * np.array([1,1,1])):
        ax.plot([point[0]], [point[1]], [point[2]], 'w')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
u, v = np.meshgrid(u,v)
rEH = 1 + math.sqrt(1-a**2)
rErg = 1+np.sqrt(1-a**2*np.cos(v)**2)
xs = rEH * np.cos(u)*np.sin(v)
ys = rEH * np.sin(u)*np.sin(v)
zs = rEH * np.cos(v)

xErg = rErg*np.cos(u)*np.sin(v)
yErg = rErg*np.sin(u)*np.sin(v)
zErg = rErg*np.cos(v)
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='black',linewidth=0, alpha=1.0)
#ax.plot_surface(xErg, yErg, zErg,  rstride=4, cstride=4, color='r',alpha=0.2,linewidth=0)

plt.savefig("trajectory.png", bbox_inches='tight')
