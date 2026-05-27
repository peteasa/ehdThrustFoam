################################################################################
# once the simulation has been run view the system/functions/minMax results with
# gnuplot -p minmax.gp
################################################################################


set terminal qt size 900,1100

set datafile commentschars "#"

set title "Time Evolution of U, p min max values"
set xlabel "Time (s)"
set ylabel "Values"
set grid

set multiplot layout 4, 1 scale 1.0, 2.0 margins 0.1, 0.9, 0.1, 0.9 spacing 0.1

# Customize line styles
set linetype  1 lc rgb "dark-violet" lw 2
set linetype  2 lc rgb "#009e73" lw 2
set linetype  3 lc rgb "#56b4e9" lw 2
set linetype  4 lc rgb "#e69f00" lw 2
set linetype  5 lc rgb "#f0e442" lw 2
set linetype  6 lc rgb "#0072b2" lw 2
set linetype  7 lc rgb "#e51e10" lw 2
set linetype  8 lc rgb "black"   lw 2
set linetype  9 lc rgb "gray50"  lw 2
set linetype cycle  9

set xrange [0:0.0001]

set ylabel "U" offset 0.1
set ytics nomirror
set y2label "p" offset 0.2
set y2tics
set logscale y2

plot "postProcessing/0/minMax/minMax.dat" using 1:2 with lines ls 1 axes x1y1 title "U min", \
     "postProcessing/0/minMax/minMax.dat" using 1:3 with lines ls 2 axes x1y1 title "U max", \
     "postProcessing/0/minMax/minMax.dat" using 1:6 with lines ls 3 axes x1y2 title "p min", \
     "postProcessing/0/minMax/minMax.dat" using 1:7 with lines ls 4 axes x1y2 title "p max"

set yrange [*:*]
unset logscale y
unset y2label
unset y2tics

set title "Time Evolution of phi min max values"

set ylabel "phi" offset 0.1

plot "postProcessing/0/minMax/minMax.dat" using 1:4 with lines ls 1 axes x1y1 title "phi min", \
     "postProcessing/0/minMax/minMax.dat" using 1:5 with lines ls 2 axes x1y1 title "phi max"

set yrange [*:*]

set title "Time Evolution of ne, nN2p min max values"

set ylabel "max" offset 0.1
set logscale y

set y2label "min" offset 0.2
set y2tics
set logscale y2

plot "postProcessing/0/minMax/minMax.dat" using 1:8 with lines ls 1 axes x1y2 title "ne min", \
     "postProcessing/0/minMax/minMax.dat" using 1:9 with lines ls 2 axes x1y1 title "ne max", \
     "postProcessing/0/minMax/minMax.dat" using 1:10 with lines ls 3 axes x1y2 title "nN2p min", \
     "postProcessing/0/minMax/minMax.dat" using 1:11 with lines ls 4 axes x1y1 title "nN2p max"

unset logscale y
unset logscale y2

set title "Time Evolution of phiE, E min max values"

set ylabel "max" offset 0.1
set logscale y
set y2label "min" offset 0.2
set logscale y2

plot "postProcessing/0/minMax/minMax.dat" using 1:12 with lines ls 1 axes x1y2 title "phiE min", \
     "postProcessing/0/minMax/minMax.dat" using 1:13 with lines ls 2 axes x1y1 title "phiE max", \
     "postProcessing/0/minMax/minMax.dat" using 1:14 with lines ls 3 axes x1y2 title "E min", \
     "postProcessing/0/minMax/minMax.dat" using 1:15 with lines ls 4 axes x1y1 title "E max"

unset multiplot
