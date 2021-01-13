#!/usr/bin/gnuplot
set logscale y
#set title "Residuals"
set format y "%.2e"
set ylabel 'Residual'
set xlabel 'Iteration'
set grid
set style line 1 lt rgb "red" lw 2
set style line 2 lt rgb "orange" lw 2
set style line 3 lt rgb "yellow" lw 2
set style line 4 lt rgb "green" lw 2
set style line 5 lt rgb "cyan" lw 2
set style line 6 lt rgb "blue" lw 2
set style line 7 lt rgb "violet" lw 2
plot "< cat log.simpleFoam | grep 'Solving for Ux' | cut -d' ' -f9 | tr -d ','" title 'Ux' ls 1 with lines,\
     "< cat log.simpleFoam | grep 'Solving for Uy' | cut -d' ' -f9 | tr -d ','" title 'Uy' ls 2 with lines,\
     "< cat log.simpleFoam | grep 'Solving for Uz' | cut -d' ' -f9 | tr -d ','" title 'Uz' ls 3 with lines,\
     "< cat log.simpleFoam | grep 'Solving for p' | awk 'NR % 3 == 0'| cut -d' ' -f9 | tr -d ','" title 'p' ls 4 with lines,\
     "< cat log.simpleFoam | grep 'Solving for omega' | cut -d' ' -f9 | tr -d ','" title 'omega' ls 5 with lines,\
     "< cat log.simpleFoam | grep 'Solving for k' | cut -d' ' -f9 | tr -d ','" title 'k' ls 6 with lines
#pause 10

set term png
set output "residuals.png"
replot
set term x11
