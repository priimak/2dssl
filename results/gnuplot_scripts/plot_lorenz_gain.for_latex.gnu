set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/lorentz_gain.tex"

set title "Absoprtion $A$ as a function of $\\omega$. $\\mu=116$, $\\alpha=0.9496$, $B=4$, $E_{\\omega}=0.1$, $E_{dc}=0$"

set xrange [0:12]
set yrange [0:0.053]
set grid
set xlabel "$\\omega$"
set ylabel "$Absorption$"

plot "B=4/absorption_E_dc=0_E_omega=0.1_B=4_mu=116_alpha=0.9496_vary_omega.data" u 3:6 w l lw 2 lc rgb "#ee0000" t '$\mu=116$', "mu=infty/E_dc=0.csv" u 1:($2/10) w p pt 5 lc rgb '#000000' t '$\mu=\infty$'