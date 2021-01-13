#!/usr/bin/gnuplot
set ylabel 'Uz'
set xlabel 'z'
set grid
set style line 1 lt rgb "red" lw 2
set style line 2 lt rgb "black" lw 3
set title 'velocity along the bottom wall'
set style arrow 1 head filled size char 4,10,25
set arrow 1 from 0.45,0 to 0.65,0 heads linestyle 2

plot 'postProcessing/singleGraph/500/line_Uz.xy' using 1:2 with lines linestyle 1
unset key

set term png
set output "Uz_wall.png"
replot
set term x11
