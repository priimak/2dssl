set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/v_dr_of_e_dc_B=diff.tex"

set title "$v_{dr}/v_{p}=f(E_{dc}/E_{*})$ $\\mu=116$ $\\alpha=0.9496$ $B=[0, 4, 8]$ $E_{\\omega}=0$"

set xrange [0:18.4]
set grid
set xlabel "$E_{dc}/E_{*}$"
set ylabel "$v_{dr}/v_{p}$"

plot "B=8/v_dr_B=8_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new2.data" u 1:5 w l t 'B=8' lw 2 lc rgb '#0000ee', "B=4/v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=4' lw 2 lc rgb '#ee0000', "B=0/v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=0 Esaki-Tsu' lw 2 lc rgb '#000000'

set output "plots/m_over_m_xk_of_e_dc_B=diff.tex"

set title "$m/m_{x,k}=f(E_{dc}/E_{*})$ $\\mu=116$ $\\alpha=0.9496$ $B=[0, 4, 8]$ $E_{\\omega}=0$"
set xrange [0:18.4]
set yrange [-0.4:1.1]
set grid
set xlabel "$E_{dc}/E_{*}$"
set ylabel "$m/m_{x,k}$"

plot "B=8/v_dr_B=8_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new2.data" u 1:9 w l t 'B=8' lw 2 lc rgb '#0000ee', "B=4/v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'B=4' lw 2 lc rgb '#ee0000', "B=0/v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'B=0' lw 2 lc rgb '#000000'
