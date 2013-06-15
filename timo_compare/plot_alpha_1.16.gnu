#!/bin/sh

gnuplot -persist <<EOF
#set terminal wxt size 1366,689 enhanced font 'Verdana,10'
set terminal pngcairo size 1024,600 enhanced font 'Verdana,10'
set output "priimak_vs_timo_take_2.png"

set ylabel "j_{dc}/j_{p}"
set xlabel "ω_{Β}τ"
set title "GaAs sample; Δ/(2k_{B}T)=1.16; Δ=60meV; d=6nm; T_{phys}=300Kelvin; α=0.9496; ω_{c}τ=4; B_{dmitri}=ω_{c}τ/sqrt(α)=4.10475"
plot "alpha=1.16_B_timo=4_jdc_vs_E_dc.data" u 1:(\$2/10) w l lw 2 lc rgb "#00bb00" t 'Timo', "../v_dr_B=4_alpha=0.9496_mu=1.16_E_omega=0_vary_E_dc.m_e_GaAs_scaled.data" u 1:5 w p pt 4 lw 2 lc rgb '#000000' t 'Dmitri'
EOF
