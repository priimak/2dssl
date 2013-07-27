set terminal epslatex 8 solid color colortext input lw 2
set size 1.4,1.2
set output "plots/v_dr_of_e_dc_B=4_mu=11_dot_6_116.tex"
print "Writing plots/v_dr_of_e_dc_B=4_mu=11_dot_6_116.tex"

set title "$v_{dr}/v_{p}=f(E_{dc}/E_{*})$ $\\mu=[1.16, 11.6, 116, 232]$ $\\alpha=0.9496$ $B=4$ $E_{\\omega}=0$"

set xrange [0:12]
set grid
set xlabel "$E_{dc}/E_{*}$"
set ylabel "$v_{dr}/v_{p}$"

plot "B=4/v_dr_B=4_alpha=0.9496_mu=232_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l lw 2 lc rgb "#daa520" t '$\mu=232$', "B=4/v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t '$\mu=116$' lw 2 lc rgb '#0000ee', "B=4/v_dr_B=4_alpha=0.9496_mu=11.6_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t '$\mu=11.6$' lw 2 lc rgb '#ee0000', "B=4/v_dr_B=4_alpha=0.9496_mu=1.16_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t '$\mu=1.16$' lw 2 lc rgb '#00dd00', "B=0/v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=0 Esaki-Tsu' lw 2 lc rgb '#000000'

set output "plots/m_over_m_xk_of_e_dc_B=4_mu=11_dot_6_116.tex"

set title "$m/m_{x,k}=f(E_{dc}/E_{*})$ $\\mu=[1.16, 11.6, 116, 232]$ $\\alpha=0.9496$ $B=4$ $E_{\\omega}=0$"
set xrange [0:12]
set yrange [-0.4:1.1]
set grid
set xlabel "$E_{dc}/E_{*}$"
set ylabel "$m/m_{x,k}$"

plot "B=4/v_dr_B=4_alpha=0.9496_mu=232_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t '$\mu=232$' lw 2 lc rgb "#daa520", "B=4/v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t '$\mu=116$' lw 2 lc rgb '#0000ee', "B=4/v_dr_B=4_alpha=0.9496_mu=11.6_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t '$\mu=11.6$' lw 2 lc rgb '#ee0000', "B=4/v_dr_B=4_alpha=0.9496_mu=1.16_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t '$\mu=1.16$' lw 2 lc rgb '#00dd00', "B=0/v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t '$\mu=116 B=0$' lw 2 lc rgb '#000000'
