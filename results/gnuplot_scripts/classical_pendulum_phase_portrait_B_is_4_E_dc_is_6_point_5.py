#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
from matplotlib.mlab      import griddata
from mpl_toolkits.mplot3d import axes3d
import numpy as np

fig  = plt.figure(figsize=(12, 6))
fig.subplots_adjust(left=0.05, right=0.95, top=0.93, wspace=0.13, hspace=0.27)
ax = fig.add_subplot(121)

data_new = np.genfromtxt('B=4/classical_pendilum_data_B_is_4_E_dc_is_6_point_5.data', delimiter=' ', comments='#',
                         names=['t', 'phi_x', 'phi_y'])
sx_new = np.genfromtxt('B=4/classical_pendilum_separatrix_B_is_4_E_dc_is_6_point_5.data', delimiter=' ', comments='#',
                         names=['t', 'phi_x', 'phi_y'])
#p1 = ax.plot(data_new['phi_x'], data_new['phi_y'], ',', lw=3, color='black', label='')
ax.scatter(data_new['phi_x'], data_new['phi_y'], marker='.', s=1, color='black')
ax.scatter(sx_new['phi_x'], sx_new['phi_y'], marker='.', s=1, color='red')
ax.plot([-3.142,3.142], [0,0], color='black', lw=1.2, ls='dashed')
ax.plot([0,0], [-3.9,1], color='black', lw=1.2, ls='dashed')

ax.set_xlabel('$\phi_x$')
ax.set_ylabel('$\phi_y$', labelpad=-6)

ax.xaxis.set_ticks([-3.14,0,3.14])
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = '$-\pi$'
labels[1] = '$0$'
labels[2] = '$\pi$'
ax.set_xticklabels(labels)

ax.yaxis.set_ticks([0])
ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.9,1)

ax = fig.add_subplot(122)

data = np.genfromtxt('B=4/f_B=4_E_dc=6.5_mu=5_alpha=0.9496.data',
                     delimiter=' ', 
                     names=['phi_x', 'phi_y', 'f'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(-3.141, 3.141, num=800)
yi = np.linspace(-3.9, 1, num=800)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['phi_x'], data['phi_y'], data['f'], xi, yi)
plt.pcolor(X, Y, Z, cmap='hot', vmin=0)

ax.xaxis.set_ticks([-3.14,0,3.14])
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = '$-\pi$'
labels[1] = '$0$'
labels[2] = '$\pi$'
ax.set_xticklabels(labels)

ax.yaxis.set_ticks([])

ax.set_xlim(-3.142,3.142)
ax.set_ylim(-3.9,1)

ax.set_xlabel('$\phi_x$')

#plt.show()
outfile="plots/classical_pendulum_phase_portrait_B_is_4_E_dc_is_6_point_5.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")
