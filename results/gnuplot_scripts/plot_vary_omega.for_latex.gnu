set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/combined_B=0_and_4_mu=116_and_1_16_E_omega=0_1.tex"
print "Writing gnuplot_scripts/plot_vary_omega.for_latex.gnu"

set title "Absoprtion $A$ as a function of $\\omega$. $\\mu=116$ $\\alpha=0.9496$ $E_{\\omega}=1$ $E_{dc}=6.4$ $B=4$"

set xrange [0:12]
set yrange [-0.015:0.0245]
set grid
set xlabel "$\\omega$"
set ylabel "$Absorption$"

#plot "absorption_E_dc=6.4_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#ee0000" t '$\mu=116, B=4$', "absorption_E_dc=6.4_E_omega=0.1_B=4_mu=1.16_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#0000ff" t '$\mu=1.16, B=4$', "absorption_E_dc=6.4_E_omega=0.1_B=0_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#000000" t '$\mu=116, B=0$ Bloch gain'

plot "B=4/absorption_E_dc=6.4_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#ee0000" t '$\mu=116$', "B=4/absorption_E_dc=6.4_E_omega=0.1_B=4_mu=1.16_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#0000ff" t '$\mu=1.16$'

#plot "E_dc=6.csv" u 1:($2/10) w l lw 2 lc rgb "#ee0000" t '$\mu=\infty$', "absorption_E_dc=6.4_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#008800" t '$\mu=116$', "absorption_E_dc=6.4_E_omega=0.1_B=4_mu=1.16_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#0000ff" t '$\mu=1.16$'