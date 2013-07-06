#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk

tilde_B=4.0
E_dc=6.0

print math.pi*tilde_B
print 2*ellipk(E_dc/(2*tilde_B))
print (math.pi*tilde_B/(2*ellipk(E_dc/(2*tilde_B))))

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

omega_resonance_1 = omega_c_low(E_dc)
omega_1 = [omega_resonance_1, omega_resonance_1]
omega_2 = [3*omega_resonance_1, 3*omega_resonance_1]
line_value = [-0.01, 0.03]

fig  = plt.figure()
fig.subplots_adjust(right=0.87, top=0.93, wspace=0.13, hspace=0.27)

data_mu_116 = np.genfromtxt('B=4/absorption_E_dc=6_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

data_mu_infty = np.genfromtxt('mu=infty/E_dc=6.csv', comments='#', delimiter=' ', names=['omega', 'A'])

ax = fig.add_subplot(111)
ax.set_title('Absoprtion $A=f(\omega)$. $\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\omega}=0.1$, $E_{dc}=6$')
ax.xaxis.set_ticks([0,3,6,9,12])
ax.yaxis.set_ticks([0])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2, label='$\mu=116$')
p2 = ax.plot(data_mu_infty['omega'], data_mu_infty['A']/10, 'bs', color='black', label='$\mu=\infty$')
p3 = ax.plot(omega_1, line_value, color='black', lw=1.2, ls='dashed')
p4 = ax.plot(omega_2, line_value, color='black', lw=1.2, ls='dashed')
p5 = ax.plot([0,12], [0,0], color='black', lw=1, ls='dashed')

ax.annotate('$\omega=\Omega$', xy=(3.4,0.025), xycoords='data', fontsize='large')
ax.annotate('$\omega=3\Omega$', xy=(10,0.02), xycoords='data', fontsize='large')
ax.set_xlabel('$\omega$')
ax.set_ylabel('$Absorption$')

ax.legend()
ax.set_xlim(0,12)
ax.set_ylim(-0.01, 0.03)
ax.grid(False)
plt.show()
