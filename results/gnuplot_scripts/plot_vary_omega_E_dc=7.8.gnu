set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/vary_omega_B=4_mu=116_E_omega=0_1_E_dc=7_point_8.tex"
print "Writing plots/vary_omega_B=4_mu=116_E_omega=0_1_E_dc=7_point_8.tex"

set title "Absoprtion $A$ as a function of $\\omega$. $\\mu=116$, $\\alpha=0.9496$, $E_{\\omega}=0.1$, $B=4$, $E_{dc}=7.8$"

set xrange [0:12]
set yrange [-0.04:0.025]
set grid
set xlabel "$\\omega$"
set ylabel "$Absorption$"

plot "B=4/absorption_E_dc=7.8_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#ee0000" t '', '< echo "5.1 0.020614"' w p pt 7 lw 2 lc rgb '#000000' t '', '< echo "2.455 -0.0380679"' w p pt 7 lw 2 lc rgb '#000000' t ''