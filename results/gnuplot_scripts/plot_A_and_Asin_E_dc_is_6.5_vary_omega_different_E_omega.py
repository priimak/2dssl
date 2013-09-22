#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special        import ellipk
from matplotlib.mlab      import griddata

fig  = plt.figure()
fig.subplots_adjust(left=0.05, right=0.97, top=0.93, wspace=0.13, hspace=0.27)
ax   = fig.add_subplot(121)
ax.yaxis.set_ticks([])

data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_6.5_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax.set_title('A')
p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_2.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='g', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_3.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='b', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_4.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='black', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_5.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='orange', lw=2)

p5 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')

ax.set_xlabel('$\\tilde{\omega}$')
ax.set_ylim(-0.075,0.08)
ax.set_xlim(0,10)

ax   = fig.add_subplot(122)

ax.set_title('ASIN')
ax.yaxis.set_ticks([])
data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_6.5_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='r', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_2.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='g', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_3.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='b', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_4.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='black', lw=2)

data_mu_116 = np.genfromtxt('B=4/result_omega_Eomega_is_0_point_5.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])
p2 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='orange', lw=2)
ax.legend(['$E_\omega=0.1$', '$E_\omega=0.2$', '$E_\omega=0.3$', '$E_\omega=0.4$', '$E_\omega=0.5$'], loc='best')

p5 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')

ax.set_xlabel('$\\tilde{\omega}$')
ax.set_xlim(0,10)
ax.set_ylim(-0.02,0.1)

ax.grid(False)

plt.show()
#outfile="plots/A_and_ASIN_of_omega_E_dc_is_6_point_5_vary_mu_different_E_omega.pdf"
#print "Writing {o}".format(o=outfile)
#plt.savefig(outfile, format="pdf")
