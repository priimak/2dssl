#!/usr/bin/env python

from matplotlib import cm
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


fig  = plt.figure()
#fig.subplots_adjust(right=0.87, top=0.93, wspace=0.13, hspace=0.27)
ax   = fig.add_subplot(111, projection='3d')

ax.set_title("$Absorption$ $for$ $\\tilde{B}=8$, $\mu=116$, $\\alpha=0.0496$, $E_{\omega}=0.1$")

data = np.genfromtxt('B=8/absorption_B_is_8_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(18, 5, num=100)
yi = np.linspace(0.8, 12, num=300)
X, Y = np.meshgrid(xi, yi)

Z = griddata(data['E_dc'], data['omega'], data['A'], xi, yi)
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

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

ax.plot_wireframe(X, Y, Z, rstride=5, color='black')
p = ax.contourf(X, Y, Z, zdir='z', offset=-0.05, cmap=my_cmap, linewidth=0, antialiased=True, hold='on', 
                levels=[-0.05,-0.02,-0.015,-0.01,-0.005,0,0.01,0.02,0.03,0.04,0.05])
ax.set_xlim(5, 18)
ax.set_xlabel('$E_{dc}$')
ax.set_ylabel('$\omega$')
ax.set_zlabel('$A$')
plt.show()
exit
