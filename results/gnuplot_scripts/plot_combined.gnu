set terminal pngcairo size 1024,600 enhanced font 'Verdana,12' 
set output "v_dr_of_e_dc_B=4_mu=11_dot_6_116.png"

set title "v_{dr}/v_{p}=f(E_{dc}/E_{*}) μ=\\{1.16, 11.6, 116, 232\\} α=0.9496 B=4"

set xrange [0:12]
set grid
set xlabel "E_{dc}/E_{*}"
set ylabel "v_{dr}/v_{p}"

plot "v_dr_B=4_alpha=0.9496_mu=232_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l lw 2 lc rgb "#daa520" t 'μ=232', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'μ=116' lw 2 lc rgb '#0000ee', "v_dr_B=4_alpha=0.9496_mu=11.6_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'μ=11.6' lw 2 lc rgb '#ee0000', "v_dr_B=4_alpha=0.9496_mu=1.16_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'μ=1.16' lw 2 lc rgb '#00dd00', "v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=0 Esaki-Tsu' lw 2 lc rgb '#000000'

set output "m_over_m_xk_of_e_dc_B=4_mu=11_dot_6_116.png"

set title "m/m_{x,k}=f(E_{dc}/E_{*}) μ=\\{1.16, 11.6, 116, 232\\} α=0.9496 B=4"
set xrange [0:12]
set yrange [-0.4:1.1]
set grid
set xlabel "E_{dc}/E_{*}"
set ylabel "m/m_{x,k}"

plot "v_dr_B=4_alpha=0.9496_mu=232_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'μ=232' lw 2 lc rgb "#daa520", "v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'μ=116' lw 2 lc rgb '#0000ee', "v_dr_B=4_alpha=0.9496_mu=11.6_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'μ=11.6' lw 2 lc rgb '#ee0000', "v_dr_B=4_alpha=0.9496_mu=1.16_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'μ=1.16' lw 2 lc rgb '#00dd00', "v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'μ=116 B=0' lw 2 lc rgb '#000000'
