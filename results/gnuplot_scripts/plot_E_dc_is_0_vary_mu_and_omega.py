#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab      import griddata

fig  = plt.figure()
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, wspace=0.13, hspace=0.27)

data = \
np.genfromtxt('B=3/B_is_3_E_omega_is_0.1_alpha_0.9496_E_dc_is_0.0_vary_mu_and_omega.data',
              delimiter=' ', comments='#',
              names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 
                     'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(121)
ax.set_title('$A$')
ax.xaxis.set_ticks([0,11])
ax.yaxis.set_ticks([0])

index = np.asarray([row['mu'] > 0.29 and row['mu'] < 0.31 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#ff0000', lw=2)

index = np.asarray([row['mu'] > 0.79 and row['mu'] < 0.81 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#00aa22', lw=2)

index = np.asarray([row['mu'] > 1.9 and row['mu'] < 2.1 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#00bbff', lw=2)

index = np.asarray([row['mu'] > 115.9 and row['mu'] < 116.1 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#0000ff', lw=2, alpha=0.9)

ax.legend(['$\mu=0.3$', '$\mu=0.8$', '$\mu=2.0$', '$\mu=116$'], loc='best')

ax.set_xlim(0,11)
ax.set_ylim(0,0.055)
ax.set_xlabel('$\\tilde{\omega}$')

ax = fig.add_subplot(122)
ax.set_title('$ASIN$')
ax.xaxis.set_ticks([0,11])
ax.yaxis.set_ticks([0])

index = np.asarray([row['mu'] > 0.29 and row['mu'] < 0.31 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#ff0000', lw=2)

index = np.asarray([row['mu'] > 0.79 and row['mu'] < 0.81 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#00aa22', lw=2)

index = np.asarray([row['mu'] > 1.9 and row['mu'] < 2.1 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#00bbff', lw=2)

index = np.asarray([row['mu'] > 115.9 and row['mu'] < 116.1 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#0000ff', lw=2, alpha=0.9)

ax.legend(['$\mu=0.3$', '$\mu=0.8$', '$\mu=2.0$', '$\mu=116$'], loc='best')

p5 = ax.plot([0,11], [0,0], color='black', lw=1.2, ls='dashed')

ax.set_xlim(0,11)
ax.set_ylim(-0.018,0.035)
ax.set_xlabel('$\\tilde{\omega}$')

plt.show()
#outfile="plots/A_and_ASIN_of_omega_E_dc_is_6_point_5_vary_mu.pdf"
#print "Writing {o}".format(o=outfile)
#plt.savefig(outfile, format="pdf")



