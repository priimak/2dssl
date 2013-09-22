#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab      import griddata

fig  = plt.figure()
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, wspace=0.13, hspace=0.27)

data = \
np.genfromtxt('omega=0.5/B_is_4_E_omega_is_0.1_alpha_0.9496_omega_is_0.5_vary_E_dc_and_omega.data',
              delimiter=' ', comments='#',
              names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 
                     'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(111)
ax.set_title('$\omega=0.5$')
ax.xaxis.set_ticks([0,7,10])
ax.yaxis.set_ticks([0,1])

index = np.asarray([row['mu'] > 1.15999995 and row['mu'] < 1.15999999 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['E_dc'], fdata['v_dr'], color='#ff0000', lw=2)

index = np.asarray([row['mu'] > 2.9 and row['mu'] < 3.1 for row in data])
fdata = data[index]
p1_1 = ax.plot(fdata['E_dc'], fdata['v_dr'], color='#999900', lw=2)

index = np.asarray([row['mu'] > 11 and row['mu'] < 12 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['E_dc'], fdata['v_dr'], color='#33aaaa', lw=2)

index = np.asarray([row['mu'] > 115 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['E_dc'], fdata['v_dr'], color='#00bbff', lw=2)

ax.legend(['$\mu=1.16$', '$\mu=3$', '$\mu=11.6$', '$\mu=116$'], loc='best')

p5 = ax.plot([0,10], [0,0], color='black', lw=1.2, ls='dashed')
p5 = ax.plot([7,7], [-1,1], color='black', lw=1.2, ls='dashed')

ax.set_xlim(0,10)
ax.set_ylim(0,1)
ax.set_ylabel('$<v_{dr}>$')
ax.set_xlabel('$E_{dc}$')

plt.show()
#outfile="plots/A_and_ASIN_of_omega_E_dc_is_6_point_5_vary_mu.pdf"
#print "Writing {o}".format(o=outfile)
#plt.savefig(outfile, format="pdf")
