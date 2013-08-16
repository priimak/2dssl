#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab      import griddata

fig  = plt.figure()
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, wspace=0.13, hspace=0.27)

data = \
np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_alpha_0.9496_E_dc_is_7.8_vary_mu_and_omega.data',
              delimiter=' ', comments='#',
              names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 
                     'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(121)
ax.set_title('$A$')
ax.xaxis.set_ticks([0,10])
ax.yaxis.set_ticks([0])

index = np.asarray([row['mu'] > 1.15999995 and row['mu'] < 1.15999999 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#ff0000', lw=2)

index = np.asarray([row['mu'] > 2.9 and row['mu'] < 3.1 for row in data])
fdata = data[index]
p1_1 = ax.plot(fdata['omega'], fdata['A'], color='#999900', lw=2)

index = np.asarray([row['mu'] > 11 and row['mu'] < 12 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#33aaaa', lw=2)

index = np.asarray([row['mu'] > 31 and row['mu'] < 33 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#00bbff', lw=2)

index = np.asarray([row['mu'] > 110 and row['mu'] < 120 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['A'], color='#0000ff', lw=2)

ax.legend(['$\mu=1.16$', '$\mu=3$', '$\mu=11.16$', '$\mu=32.16$', '$\mu=116$'], loc='best')

p5 = ax.plot([0,10], [0,0], color='black', lw=1.2, ls='dashed')

ax.set_xlim(0,10)
ax.set_ylim(-0.04,0.023)
ax.set_xlabel('$\\tilde{\omega}$')

ax = fig.add_subplot(122)
ax.set_title('$ASIN$')
ax.xaxis.set_ticks([0,10])
ax.yaxis.set_ticks([0])

index = np.asarray([row['mu'] > 1.15999995 and row['mu'] < 1.15999999 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#ff0000', lw=2)

index = np.asarray([row['mu'] > 2.9 and row['mu'] < 3.1 for row in data])
fdata = data[index]
p1_1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#999900', lw=2)

index = np.asarray([row['mu'] > 11 and row['mu'] < 12 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#33aaaa', lw=2)

index = np.asarray([row['mu'] > 31 and row['mu'] < 33 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#00bbff', lw=2)

index = np.asarray([row['mu'] > 110 and row['mu'] < 120 for row in data])
fdata = data[index]
p1 = ax.plot(fdata['omega'], fdata['ASIN'], color='#0000ff', lw=2)

ax.legend(['$\mu=1.16$', '$\mu=3$', '$\mu=11.16$', '$\mu=32.16$', '$\mu=116$'], loc='best')

p5 = ax.plot([0,10], [0,0], color='black', lw=1.2, ls='dashed')

ax.set_xlim(0,10)
ax.set_ylim(-0.046,0.006)
ax.set_xlabel('$\\tilde{\omega}$')

#plt.show()
outfile="plots/A_and_ASIN_of_omega_E_dc_is_7_point_8_vary_mu.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")



