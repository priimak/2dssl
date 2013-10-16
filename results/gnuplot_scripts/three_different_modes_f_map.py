#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
from matplotlib.mlab      import griddata
from mpl_toolkits.mplot3d import axes3d
import numpy as np

fig  = plt.figure(figsize=(18, 6))
fig.subplots_adjust(left=0.0, right=1, top=1, bottom=0, wspace=0.08, hspace=0.09)

ax = fig.add_subplot(131)
data = np.genfromtxt('B=4/E_dc=6_mu=116_alpha=0.9496_B=4.0_frame.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-4.5, 0.6, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-4.5,0.6)
ax.text(-2.6, 0, '(a)', color='white', fontsize=28)

ax = fig.add_subplot(132)
data = np.genfromtxt('B=4/E_dc=8_mu=116_alpha=0.9496_B=4.0_frame.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-4.5, 0.6, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-4.5,0.6)
ax.text(-2.6, 0, '(b)', color='white', fontsize=28)

ax = fig.add_subplot(133)
data = np.genfromtxt('B=4/E_dc=10_mu=116_alpha=0.9496_B=4.0_frame.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-4.5, 0.6, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-4.5,0.6)
ax.text(-2.6, 0, '(c)', color='white', fontsize=28)

#plt.show()
outfile="plots/three_types_of_f.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")
