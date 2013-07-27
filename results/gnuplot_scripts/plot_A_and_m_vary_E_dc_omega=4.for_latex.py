#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

fig  = plt.figure()
fig.subplots_adjust(right=0.87, top=0.93, wspace=0.13, hspace=0.27)
ax1   = fig.add_subplot(221)
ax2 = ax1.twinx()

data = np.genfromtxt('B=4/absorption_omega=4_E_omega=0.1_B=4_mu=216_alpha=0.9496_vary_E_dc.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

ax1.plot(data['E_dc'], data['A'], color='r', lw=2, label='Absorption')
ax1.fill_between(data['E_dc'], data['A'], 0, interpolate=True, where=data['A']<0, color='#ff0000', alpha=0.3)
ax1.set_title('(a) $\\mu=216$')
ax1.annotate('$A$', xy=(6.65,0.0215), xycoords='data', fontsize='large')
ax1.set_xlim(0,11)
ax1.set_ylim(-0.04,0.055)
ax1.yaxis.set_ticks([-0.04, 0, 0.04])
ax1.plot([0,11],[0,0], ls='dotted', color='black')
ax1.set_ylabel('$A$', fontsize='x-large')

ax2.set_ylim(-0.4, 1)
ax2.yaxis.set_ticks([])
ax2.plot(data['E_dc'], data['m_eff'], color='b', lw=2, label='m_{x}', ls='dashed')
ax2.plot([0,11],[0,0], ls='dotted', color='black')
ax2.annotate('$\\frac{m_*}{m_x}$', xy=(5,0.59), xycoords='data', fontsize='large')
ax2.fill_between(data['E_dc'], data['m_eff'], 0, interpolate=True, where=data['m_eff']<0, color='#0000ff', alpha=0.3)

ax1   = fig.add_subplot(222)
ax2 = ax1.twinx()

# second plot 

data = np.genfromtxt('B=4/absorption_omega=4_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_E_dc.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

ax1.plot(data['E_dc'], data['A'], color='r', lw=2, label='Absorption')
ax1.fill_between(data['E_dc'], data['A'], 0, interpolate=True, where=data['A']<0, color='#ff0000', alpha=0.3)
ax1.annotate('$A$', xy=(6.54,0.0215), xycoords='data', fontsize='large')
ax1.set_title('(b) $\\mu=116$')
ax1.set_xlim(0,11)
ax1.set_ylim(-0.04,0.055)
ax1.yaxis.set_ticks([])
ax1.plot([0,11],[0,0], ls='dotted', color='black')
ax1.grid(False)

ax2.set_ylabel('$\\frac{m_*}{m_x}$', fontsize='x-large')
ax2.set_ylim(-0.4, 1)
ax2.yaxis.set_ticks([-0.4, 0, 1])
ax2.plot(data['E_dc'], data['m_eff'], color='b', lw=2, label='m_{x}', ls='dashed')
ax2.plot([0,11],[0,0], ls='dotted', color='black')
ax2.annotate('$\\frac{m_*}{m_x}$', xy=(5,0.59), xycoords='data', fontsize='large')
ax2.fill_between(data['E_dc'], data['m_eff'], 0, interpolate=True, where=data['m_eff']<0, color='#0000ff', alpha=0.3)

# third plot

ax1   = fig.add_subplot(223)
ax2 = ax1.twinx()

data = np.genfromtxt('B=4/absorption_omega=4_E_omega=0.1_B=4_mu=11.6_alpha=0.9496_vary_E_dc.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

ax1.plot(data['E_dc'], data['A'], color='r', lw=2, label='Absorption')
ax1.fill_between(data['E_dc'], data['A'], 0, interpolate=True, where=data['A']<0, color='#ff0000', alpha=0.3)
ax1.set_title('(c) $\\mu=11.6$')
ax1.annotate('$A$', xy=(5.7,0.018), xycoords='data', fontsize='large')
ax1.set_xlim(0,11)
ax1.set_ylim(-0.04,0.055)
ax1.yaxis.set_ticks([-0.04, 0, 0.04])
ax1.plot([0,11],[0,0], ls='dotted', color='black')
ax1.set_ylabel('$A$', fontsize='x-large')
ax1.set_xlabel('$E_{dc}/E_{*}$', fontsize='large')

ax2.set_ylim(-0.4, 1)
ax2.yaxis.set_ticks([])
ax2.plot(data['E_dc'], data['m_eff'], color='b', lw=2, label='m_{x}', ls='dashed')
ax2.plot([0,11],[0,0], ls='dotted', color='black')
ax2.annotate('$\\frac{m_*}{m_x}$', xy=(5,0.2), xycoords='data', fontsize='large')
ax2.fill_between(data['E_dc'], data['m_eff'], 0, interpolate=True, where=data['m_eff']<0, color='#0000ff', alpha=0.3)

# forth plot 
ax1   = fig.add_subplot(224)
ax2 = ax1.twinx()

data = np.genfromtxt('B=4/absorption_omega=4_E_omega=0.1_B=4_mu=1.16_alpha=0.9496_vary_E_dc.data', 
                     delimiter=' ', 
                     names=['E_dc', 'E_omega', 'omega', 'mu', 'v_dr_inst', 'A', 'NORM', 'v_y_inst', 'm_eff_inst', 'v_dr', 'v_y', 'm_eff'])

ax1.plot(data['E_dc'], data['A'], color='r', lw=2, label='Absorption')
ax1.fill_between(data['E_dc'], data['A'], 0, interpolate=True, where=data['A']<0, color='#ff0000', alpha=0.3)
ax1.annotate('$A$', xy=(5.7,0.01), xycoords='data', fontsize='large')
ax1.set_title('(d) $\\mu=1.16$')
ax1.set_xlim(0,11)
ax1.set_ylim(-0.04,0.055)
ax1.yaxis.set_ticks([])
ax1.plot([0,11],[0,0], ls='dotted', color='black')
ax1.set_xlabel('$E_{dc}/E_{*}$', fontsize='large')
ax1.grid(False)

ax2.set_ylabel('$\\frac{m_*}{m_x}$', fontsize='x-large')
ax2.set_ylim(-0.4, 1)
ax2.yaxis.set_ticks([-0.4, 0, 1])
ax2.plot(data['E_dc'], data['m_eff'], color='b', lw=2, label='m_{x}', ls='dashed')
ax2.plot([0,11],[0,0], ls='dotted', color='black')
ax2.annotate('$\\frac{m_*}{m_x}$', xy=(4.8,0.1), xycoords='data', fontsize='large')
ax2.fill_between(data['E_dc'], data['m_eff'], 0, interpolate=True, where=data['m_eff']<0, color='#0000ff', alpha=0.3)

#plt.show()
outfile="plots/A_and_m_vary_E_dc_and_mu_omega=4.pdf"
print "Writing {o}".format(o=outfile)
plt.savefig(outfile, format="pdf")




