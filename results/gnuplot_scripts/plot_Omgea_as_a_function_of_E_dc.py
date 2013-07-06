#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk

#data = np.genfromtxt('mu=infty/Omega_of_E_dc_B=4.data', delimiter=' ', names=['omega', 'OOmega'])
data = np.genfromtxt('results/mu=infty/PRL1.csv', delimiter=',', names=['omega', 'OOmega'])

tilde_B=4.0
E_dc=6.0

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

E_dc_cyclotron=np.linspace(0, 2*tilde_B-0.01, num=200)
E_dc_bloch=np.linspace(2*tilde_B+0.01, 12, num=100)
OmegaCLow=map(omega_c_low, E_dc_cyclotron)
OmegaCHigh=map(omega_c_high, E_dc_bloch)
E_dc_cyclotron = np.append(E_dc_cyclotron, 2*tilde_B)
E_dc_bloch = np.insert(E_dc_bloch, 0, 2*tilde_B)
OmegaCHigh = np.insert(OmegaCHigh, 0, 0)
OmegaCLow = np.append(OmegaCLow, 0)

fig  = plt.figure()
fig.subplots_adjust(right=0.87, top=0.93, wspace=0.13, hspace=0.27)
ax = fig.add_subplot(111)

p1 = ax.plot(E_dc_cyclotron, OmegaCLow, color='black', lw=2, label='$\mu=116$')
p2 = ax.plot(E_dc_bloch, OmegaCHigh, color='black', lw=2)
p3 = ax.plot(data['omega'], data['OOmega'], 'ro', label='From PRL Fig.1(a) $\omega_{c}\\tau = 4$')
ax.set_ylim(0, 10)
ax.set_xlim(0,12)
ax.grid(True)
ax.set_xlabel('$\omega$')
ax.set_ylabel('$\Omega$')

ax.legend(loc=0)

plt.show()
