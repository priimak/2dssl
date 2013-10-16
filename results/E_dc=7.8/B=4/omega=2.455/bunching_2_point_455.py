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

data = np.genfromtxt('f_x_time.data.bz2', 
                     delimiter=' ', 
                     names=['t', 'phi_x', 'f'])

fig  = plt.figure()
ax   = fig.add_subplot(111)

xi = np.linspace(40, 60, num=100)
yi = np.linspace(-3.141, 3.141, num=100)
X, Y = np.meshgrid(xi, yi)
F = griddata(2.455*data['t'], data['phi_x'], data['f'], xi, yi)

#plt.pcolor(X, Y, F, cmap='RdBu', vmin=0)
#plt.colorbar()
#ax.set_ylim([-3.141,3.141])


#levels = MaxNLocator(nbins=30)#.tick_values(F.min(), F.max())

df=(F.max()-F.min())/100
levels = np.arange(F.min()-df,F.max()+df,df)
#p = ax.contourf(X, Y, F, linewidth=0, antialiased=True, levels=levels)
p = ax.contourf(X, Y, F, antialiased=False, levels=levels)
#plt.colorbar()

plt.show()

