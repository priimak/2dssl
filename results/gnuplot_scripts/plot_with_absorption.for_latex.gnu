set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "combined_B=4_mu=116_omega=5_E_omega=1.tex"

set title "$\\mu=116$ $\\alpha=0.9496$ $E_{\\omega}=1$ $\\omega=5$ $B=4$"

set xrange [0:12]
set yrange [-0.26:0.92]
set grid
set xlabel "$E_{dc}/E_{*}$"

plot "v_dr_B=4_alpha=0.9496_mu=116_E_omega=1_omega=5_vary_E_dc.m_e_GaAs_scaled.data" u 1:6 w l lw 2 lc rgb "#ee0000" t '$A$', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=1_omega=5_vary_E_dc.m_e_GaAs_scaled.data" u 1:12 w l lw 2 lc rgb '#0000ee' t '$<m/m_{x,k}>$', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=1_omega=5_vary_E_dc.m_e_GaAs_scaled.data" u 1:10 w l lw 2 lc rgb '#00dd00' t '$<v_{dr}/v_{p}>$', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l lw 2 lc rgb "#000000" t '$v_{dr}/v_{p}$ for $\omega=0$'
