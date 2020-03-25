#!/bin/bash

writeout=$1


/usr/bin/gnuplot << EOF

#readin = "1d_xyz_dissipativ_site6_Jy1.2_Nbr5000_h0.out"
#writeout = "1d_xyz_dissipativ_site6_Jy1.2_Nbr5000_h0.png"

readin = "${writeout}.out"


set terminal pngcairo size 2048, 1024
set output "${writeout}.png"

plot readin using 1:4
MAX=GPVAL_DATA_X_MAX

set output "${writeout}.png"
set multiplot layout 2,2 rowsfirst

set xrange [0:MAX]
set title "Shift"
plot readin using 1:4 with lines lt rgb "#1E90FF" notitle

set xrange [0:MAX]
set title "Walker number - Diag. walk total"
plot readin using 1:2 with lines lt rgb "#1E90FF" notitle, readin using 1:3 with lines lt rgb "#FF0000" notitle

set xrange [0:MAX]
set title "Diagonal walk: Real - Imag"
plot readin using 1:9 with lines lt rgb "#1E90FF" notitle, readin using 1:10 with lines lt rgb "#FF0000" notitle

set xrange [0:MAX]
set yrange [-0.52:-0.46]
set title "Mz: Real - Imag"
plot readin using 1:11 with lines lt rgb "#1E90FF" notitle, readin using 1:12 with lines lt rgb "#FF0000" notitle

unset multiplot

pause mouse key

EOF

