#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk

tilde_B=4.0
E_dc=0.0

print math.pi*tilde_B
print 2*ellipk(E_dc/(2*tilde_B))
print (math.pi*tilde_B/(2*ellipk(E_dc/(2*tilde_B))))

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

omega_resonance_1 = omega_c_low(E_dc)
omega_1 = [omega_resonance_1, omega_resonance_1]
line_value = [-0.02, 0.06]

fig  = plt.figure()
fig.subplots_adjust(right=0.87, top=0.93, wspace=0.13, hspace=0.27)

data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_0_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

data_mu_infty = np.genfromtxt('mu=infty/E_dc=0.csv', comments='#', delimiter=' ', names=['omega', 'A'])

ax = fig.add_subplot(111)
ax.set_title('Absoprtion $A=f(\omega)$. $\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\omega}=0.1$, $E_{dc}=0$')
ax.xaxis.set_ticks([0,3,6,9,12])
ax.yaxis.set_ticks([0])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2, label='A $\mu=116$')
p1 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='blue', lw=2, label='ASIN $\mu=116$')
p2 = ax.plot(data_mu_infty['omega'], data_mu_infty['A']/10, 'bs', color='black', label='$\mu=\infty$')
p2[0].set_markersize(4)
p3 = ax.plot(omega_1, line_value, color='black', lw=1.2, ls='dashed')
p3 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')
ax.annotate('$\omega=\Omega$', xy=(4.1,0.055), xycoords='data', fontsize='large')


ax.set_xlabel('$\omega$')
ax.set_ylabel('$A$ and $ASIN$')

ax.legend()
ax.set_xlim(0,12)
ax.set_ylim(-0.02, 0.06)
ax.grid(False)
plt.show()
