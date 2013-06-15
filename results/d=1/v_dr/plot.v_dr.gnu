set terminal pngcairo size 1024,600 enhanced font 'Verdana,12'
set output "v_dr_vary_E_dc_B=0:1:4_E_omega=0.png"

set xlabel "E_{dc}/E_{*}≡ω_{b}τ"
set ylabel "v_{dr}/v_{0}"
set xrange [0:7.7]
set yrange [0:0.45]
set ytics 0,0.1,0.5
unset grid
set label 1 "0" at 0.5,0.35 left front font 'Times Bold,13'
set label 2 "1" at 1.6,0.438 left front font 'Times Bold,13'
set label 3 "2" at 2.8,0.41 left front font 'Times Bold,13'
set label 4 "3" at 4,0.35 left front font 'Times Bold,13'
set label 5 "4" at 5.4,0.288 left front font 'Times Bold,13'
set label 6 "T=0.2 d=1" at 5.1,0.4 left front font 'Times,16'
plot "v_dr_B=0_vary_E_dc.cuda.t_max=20.data" u 1:5 w l lw 2 lc rgb "#ff0000" t '', "v_dr_B=1_vary_E_dc.cuda.data" u 1:5 w l lw 2 lc rgb "#0000aa" t '', "v_dr_B=2_vary_E_dc.cuda.data" u 1:5 w l lw 2 lc rgb "#0000aa" t '', "v_dr_B=3_vary_E_dc.cuda.data" u 1:5 w l lw 2 lc rgb "#0000aa" t '', "v_dr_B=4_vary_E_dc.cuda.data" u 1:5 w l lw 2 lc rgb "#0000aa" t ''
