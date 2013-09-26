#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
from matplotlib.mlab      import griddata
from mpl_toolkits.mplot3d import axes3d
import numpy as np

fig  = plt.figure(figsize=(12, 6))
fig.subplots_adjust(left=0.0, right=1, top=1, bottom=0, wspace=0.08, hspace=0.09)

ax = fig.add_subplot(231)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f0.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '0', color='white', fontsize=18)

ax = fig.add_subplot(232)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f1.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '1', color='white', fontsize=18)

ax = fig.add_subplot(233)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f2.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '2', color='white', fontsize=18)

ax = fig.add_subplot(234)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f3.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '3', color='white', fontsize=18)

ax = fig.add_subplot(235)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f4.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '4', color='white', fontsize=18)

ax = fig.add_subplot(236)
data = np.genfromtxt('E_dc=7.0/frame_E_dc=7_E_omega=0_mu=11.6_alpha=0.9496_B=4.0_f5.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=400)
yi = np.linspace(-3.9, 1, num=400)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.91)
ax.text(-2.6, 0.2, '5', color='white', fontsize=18)

#plt.show()
outfile="plots/example_of_time_evolution_of_f.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")
