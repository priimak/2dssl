set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/combined_B=4_mu=116_E_omega=0_1_E_dc=6_vs_6_point_4.tex"
print "Writing plots/combined_B=4_mu=116_E_omega=0_1_E_dc=6_vs_6_point_4.tex"

set title "Absoprtion $A$ as a function of $\\omega$. $\\mu=116$, $\\alpha=0.9496$, $E_{\\omega}=0.1$, $B=4$, $E_{dc}=[6,6.4]$"

set xrange [0:12]
set yrange [-0.015:0.03]
set grid
set xlabel "$\\omega$"
set ylabel "$Absorption$"

plot "B=4/absorption_E_dc=6.4_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#ee0000" t '$E_{dc}=6.4$', "B=4/absorption_E_dc=6_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#0000ff" t '$E_{dc}=6$'