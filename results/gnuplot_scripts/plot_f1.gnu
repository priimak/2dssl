set terminal pngcairo size 1024,600 enhanced font 'Verdana,10'
#set terminal pngcairo size 2048,1200 enhanced font 'Verdana,10'
set output "f1.png"

set lmargin at screen 0.04
set rmargin at screen 0.89
set bmargin at screen 0.08
set tmargin at screen 0.95

set title 'Frame 39'

set view map
set xrange [-3.14159:3.14159]
set yrange [-4.5:4.5]
splot "f1.data" u 1:2:3 w image
