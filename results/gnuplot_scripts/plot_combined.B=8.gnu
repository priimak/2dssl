set terminal epslatex 8 solid color colortext input lw 2
#set size 0.75,0.75

#set terminal epslatex size 2,2 standalone color colortext 10
set output 'sample1.tex'
#set terminal pngcairo size 1024,600 enhanced font 'Verdana,12' 

#set output "v_dr_of_e_dc_B=8_mu=set.png"
#set output "v_dr_of_e_dc_B=8_mu=set.tex"

set title "v_{dr}/v_{p}=f(E_{dc}/E_{*}) μ=116 α=0.9496 B=\\{0, 4, 8\\}"

set xrange [0:18.4]
set grid
set xlabel "E_{dc}/E_{*}"
set ylabel "v_{dr}/v_{p}"

plot "v_dr_B=8_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new2.data" u 1:5 w l t 'B=8' lw 2 lc rgb '#0000ee', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=4' lw 2 lc rgb '#ee0000', "v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:5 w l t 'B=0 Esaki-Tsu' lw 2 lc rgb '#000000'

#set output "m_over_m_xk_of_e_dc_B=8_mu=set.png"
set output "m_over_m_xk_of_e_dc_B=8_mu=set.tex"

set title "m/m_{x,k}=f(E_{dc}/E_{*}) μ=116 α=0.9496 B=\\{0, 4, 8\\}"
set xrange [0:18.4]
set yrange [-0.4:1.1]
set grid
set xlabel "E_{dc}/E_{*}"
set ylabel "m/m_{x,k}"

plot "v_dr_B=8_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new2.data" u 1:9 w l t 'B=8' lw 2 lc rgb '#0000ee', "v_dr_B=4_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'B=4' lw 2 lc rgb '#ee0000', "v_dr_B=0_alpha=0.9496_mu=116_E_omega=0_vary_E_dc.m_e_GaAs_scaled.new.data" u 1:9 w l t 'B=0' lw 2 lc rgb '#000000'
