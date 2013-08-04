#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special        import ellipk
from mpl_toolkits.mplot3d import axes3d
from matplotlib.mlab      import griddata

fig  = plt.figure()
fig.subplots_adjust(left=0.06, right=0.96, top=0.94, bottom=0.05, wspace=0.13, hspace=0.27)

ax   = fig.add_subplot(221)
ax.set_title("(a)")
plt.yticks([0.5, 10])
plt.xticks([5, 8, 10])

data = np.genfromtxt('B=4/absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

X, Y, Z = axes3d.get_test_data(0.05)

xi = np.linspace(10, 5, num=100)
yi = np.linspace(0.5, 10, num=300)

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
p = ax.plot([6.5, 6.5], [0.5,10], lw=1.5, color='black', linestyle='--')
p = ax.plot([7.8, 7.8], [0.5,10], lw=1.5, color='black', linestyle='--')
p = ax.plot([8.2, 8.2], [0.5,10], lw=1.5, color='black', linestyle='--')
p = ax.plot([9.5, 9.5], [0.5,10], lw=1.5, color='black', linestyle='--')
ax.annotate('(b)', xy=(6.55,7.7), xycoords='data', fontsize='medium')
ax.annotate('(e)', xy=(7.45,6.7), xycoords='data', fontsize='medium')
ax.annotate('(c)', xy=(8.25,7.7), xycoords='data', fontsize='medium')
ax.annotate('(d)', xy=(9.55,7.7), xycoords='data', fontsize='medium')
ax.set_xlabel('$\\tilde{E}_{dc}$', labelpad=-10)
ax.set_ylabel('$\\tilde{\omega}$', labelpad=-20)
ax.set_xlim(5,10)

##################################################### E_dc=6.5 ##################################################### 
tilde_B=4.0

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

omega_resonance_1 = omega_c_low(6.5)
omega_1 = [omega_resonance_1, omega_resonance_1]
omega_2 = [omega_resonance_1*3, omega_resonance_1*3]
omega_3 = [omega_resonance_1*3, omega_resonance_1*3]
line_value = [-0.02, 0.06]

data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_6.5_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(222)
#ax.set_title('$\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\omega}=0.1$, $E_{dc}=6.5$')
ax.set_title('(b)')
ax.xaxis.set_ticks([0,10])
ax.yaxis.set_ticks([0])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2)
p1 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='blue', lw=2)
p3 = ax.plot(omega_1, line_value, color='black', lw=1.2, ls='dashed')
p4 = ax.plot(omega_2, line_value, color='black', lw=1.2, ls='dashed')
p5 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')
ax.annotate('$A$', xy=(1.3,0.015), xycoords='data', fontsize='medium')
ax.annotate('$ASIN$', xy=(3.5,0.019), xycoords='data', fontsize='medium')
ax.annotate('$\Omega$', xy=(2.59,-0.015), xycoords='data', fontsize='medium')
ax.annotate('$2\Omega$', xy=(8.53,-0.015), xycoords='data', fontsize='medium')

ax.set_xlabel('$\\tilde{\omega}$', labelpad=-10)
#ax.set_xlabel('$\omega$')
#ax.set_ylabel('$A$ and $ASIN$')

ax.set_xlim(0,10)
ax.set_ylim(-0.018, 0.033)
ax.grid(False)

##################################################### E_dc=8.2 ##################################################### 
tilde_B=4.0

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

omega_resonance_1 = omega_c_high(8.2)

omega_1 = [omega_resonance_1, omega_resonance_1]
line_value = [-0.02, 0.06]

data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_8.2_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(223)
#ax.set_title('$\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\omega}=0.1$, $E_{dc}=6.5$')
ax.set_title('(c)')
ax.xaxis.set_ticks([0,10])
ax.yaxis.set_ticks([0])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2)
p1 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='blue', lw=2)
p3 = ax.plot(omega_1, [-0.038, 0.03], color='black', lw=1.2, ls='dashed')
p5 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')
ax.annotate('$A$', xy=(6.5,0.015), xycoords='data', fontsize='medium')
ax.annotate('$ASIN$', xy=(5.5,-0.02), xycoords='data', fontsize='medium')
ax.annotate('$\Omega$', xy=(3.8,0.02), xycoords='data', fontsize='medium')

ax.set_xlabel('$\\tilde{\omega}$', labelpad=-10)

ax.set_xlim(0,10)
ax.set_ylim(-0.038, 0.03)
ax.grid(False)

##################################################### E_dc=9.5 ##################################################### 
tilde_B=4.0

def omega_c_low(E): return math.pi*tilde_B/(2*ellipk((E/(2*tilde_B))**2))
def omega_c_high(E): return math.pi*E/(2*ellipk((2*tilde_B/E)**2))

omega_resonance_1 = omega_c_high(9.2)

omega_1 = [omega_resonance_1, omega_resonance_1]
omega_2 = [omega_resonance_1*3, omega_resonance_1*3]
omega_3 = [omega_resonance_1*3, omega_resonance_1*3]
line_value = [-0.02, 0.06]

data_mu_116 = np.genfromtxt('B=4/B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_9.5_vary_omega.data',
                            delimiter=' ', comments='#',
                            names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff', 'ASIN'])

ax = fig.add_subplot(224)
#ax.set_title('$\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\omega}=0.1$, $E_{dc}=6.5$')
ax.set_title('(d)')
ax.xaxis.set_ticks([0,10])
ax.yaxis.set_ticks([0])

p1 = ax.plot(data_mu_116['omega'], data_mu_116['A'], color='r', lw=2)
p1 = ax.plot(data_mu_116['omega'], data_mu_116['ASIN'], color='blue', lw=2)
p3 = ax.plot(omega_1, [-0.038, 0.03], color='black', lw=1.2, ls='dashed')
p4 = ax.plot(omega_2, [-0.038, 0.03], color='black', lw=1.2, ls='dashed')
p5 = ax.plot([0,12], [0,0], color='black', lw=1.2, ls='dashed')
ax.annotate('$A$', xy=(3.9,-0.015), xycoords='data', fontsize='medium')
ax.annotate('$ASIN$', xy=(8.0,-0.02), xycoords='data', fontsize='medium')
ax.annotate('$\Omega$', xy=(6.1,0.02), xycoords='data', fontsize='medium')

ax.set_xlabel('$\\tilde{\omega}$', labelpad=-10)

ax.set_xlim(0,10)
ax.set_ylim(-0.042, 0.03)
ax.grid(False)

outfile="plots/B=4_map_and_3_A_and_Asin_plots.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")
#plt.show()
