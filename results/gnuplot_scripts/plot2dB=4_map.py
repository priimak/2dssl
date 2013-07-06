#!/usr/bin/env python

from matplotlib           import cm
from matplotlib.mlab      import griddata
from mpl_toolkits.mplot3d import axes3d
from scipy.special        import ellipk

import matplotlib.pyplot as plt
import math
import numpy as np

tilde_B=4

def omega_c_low(E):  return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

fig  = plt.figure()
fig.subplots_adjust(left=0.1, top=0.95)
ax   = fig.add_subplot(111)
ax.set_title("$Absorption$ $map$ $for$ $\\tilde{B}=4$, $\mu=116$, $\\alpha=0.0496$, $E_{\omega}=0.1$")

plt.yticks([0.8, 2, 4, 6, 8, 10])

data = np.genfromtxt('B=4/absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(10, 5, num=100)
yi = np.linspace(0.8, 10, num=300)

E_dc_cyclotron = np.linspace(5, 2*tilde_B-0.01, num=200)
E_dc_bloch     = np.linspace(2*tilde_B+0.01, 10, num=100)
OmegaCLow      = map(omega_c_low, E_dc_cyclotron)
OmegaCHigh     = map(omega_c_high, E_dc_bloch)
E_dc_cyclotron = np.append(E_dc_cyclotron, 2*tilde_B)
E_dc_bloch     = np.insert(E_dc_bloch, 0, 2*tilde_B)
OmegaCHigh     = np.insert(OmegaCHigh, 0, 0.8)
OmegaCLow      = np.append(OmegaCLow, 0.8)

X, Y = np.meshgrid(xi, yi)
Z = griddata(data['E_dc'], data['omega'], data['m_eff'], xi, yi)
A = griddata(data['E_dc'], data['omega'], data['A'], xi, yi)
S = griddata(data['E_dc'], data['omega'], [0]*len(data['E_dc']), xi, yi)

from pylab import *
cdict = {'red': ((0.0, 0.0, 0.0),
                  (0.5, 1.0, 0.7),
                  (1.0, 1.0, 1.0)),
          'green': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 0.0),
                    (1.0, 1.0, 1.0)),
          'blue': ((0.0, 0.0, 0.0),
                   (0.5, 1.0, 0.0),
                   (1.0, 0.5, 1.0))}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

p = ax.contourf(X, Y, A, zdir='z', cmap=my_cmap, linewidth=0, antialiased=True, hold='on',
                levels=[-0.05,-0.02,-0.015,-0.01,-0.005,0,0.01,0.02,0.03,0.04,0.05])
p = ax.contour(X, Y, Z, zdir='z', cmap=my_cmap, linewidth=0, antialiased=True, hold='on', 
               levels=[0,0.9], linestyles='dashed')
plt.clabel(p, fmt='$m/m_x$ = %1.3f', inline=1, fontsize=11)

p = ax.plot(E_dc_cyclotron, OmegaCLow,  lw=1.5, color='black')
p = ax.plot(E_dc_bloch,     OmegaCHigh, lw=1.5, color='black')
ax.annotate('$\Omega$', xy=(7.8, 7), xytext=(5.5, 3.0),  fontsize=15)
ax.annotate('$\Omega$', xy=(7.8, 7), xytext=(8.83, 5.7), fontsize=15)

ax.set_xlabel('$E_{dc}$')
ax.set_ylabel('$\omega$')
plt.show()
exit

