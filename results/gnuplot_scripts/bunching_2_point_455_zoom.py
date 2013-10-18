#!/usr/bin/env python

from matplotlib           import cm
from matplotlib.mlab      import griddata
from mpl_toolkits.mplot3d import axes3d
from scipy.special        import ellipk
from matplotlib.ticker    import MaxNLocator
from matplotlib.colors    import BoundaryNorm

import matplotlib.pyplot as plt
import math
import numpy as np

data = np.genfromtxt('E_dc=7.8/B=4/omega=2.455/f_x_time.data.bz2', 
                     delimiter=' ', 
                     names=['t', 'phi_x', 'f'])
data2=np.genfromtxt('E_dc=7.8/B=4/E_dc=7.8,E_omega=0.1,omega=2.455,mu=116,alpha=0.9496,B=4.time_evolution.data', 
                     delimiter=' ',
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'A_inst', 't'])

fig  = plt.figure(figsize=(12, 6))
fig.subplots_adjust(left=0.08, top=0.97, right=1.02)
ax   = fig.add_subplot(111)

xi = np.linspace(50, 90, num=100)
yi = np.linspace(1.987, 3.141, num=200)
X, Y = np.meshgrid(xi, yi)
F = griddata(2.455*data['t'], data['phi_x'], data['f'], xi, yi)

df=(F.max()-F.min())/150
levels = np.arange(F.min()-df,F.max()+df,df)
fmap = ax.contourf(X, Y, F, antialiased=False, levels=levels)
ax.set_xlim([50,90])
ax.set_ylim([1.8,3.12])
ax.yaxis.set_ticks([2,3])
ax.grid()

t = np.arange(50, 90, 0.01)
e_omega=0.089*np.cos(t)+1.8935
ax.plot(t, e_omega, lw=2, color='black')
ax.plot([50,90], [1.987,1.987], lw=1.5, color='black')
ax.plot([50,90], [1.8935,1.8935], '--', lw=1, color='black')
ax.plot(data2['t']*2.455, 1.63*data2['v_dr_inst']+1.08, lw=2, color='red')
for x in range(0, 7):
    ax.plot([54.96+x*5,54.96+x*5], [1.987,1.996], lw=1.8, color='black')

ax.set_xlabel('$\\tilde{\omega} t$', fontsize='x-large')
ax.text(-0.05, 0.59, '$f(\phi_x)$',
         horizontalalignment='center',
         fontsize='x-large',
         rotation='vertical',
         transform = ax.transAxes)
cbar = plt.colorbar(fmap, ticks=[0.1275, 0.474])
cbar.ax.set_yticklabels(['Low', 'High'])

#plt.show()
outfile="plots/bunching_2_point_455_zoom.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")
