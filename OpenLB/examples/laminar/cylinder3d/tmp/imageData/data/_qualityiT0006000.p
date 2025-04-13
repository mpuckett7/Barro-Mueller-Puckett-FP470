if (strstrt(GPVAL_TERMINALS, 'jpeg') > 0) {set terminal jpeg size 5825,1000font ",25"
set output './tmp/imageData/_qualityiT0006000.jpeg'
} else {set terminal png size 5825,1000font ",25"
set output './tmp/imageData/_qualityiT0006000.png'
}
set pm3d map
unset key
set xtics out
set ytics out
set xtics nomirror
set ytics nomirror
set pm3d interpolate 0,0
set size ratio -1
set size 0.925,1.0
set colorbox vertical user origin 0.85,0.1 size 0.00429185 ,0.8
set xlabel "x-axis in m "
set ylabel "y-axis in m "
set cblabel offset 0.5 "Discretization(Rounding(refinementMetricKnudsen))"
set autoscale fix
set palette defined ( 0 "black", 1 "red", 2 "yellow")
splot './tmp/imageData/data/_qualityiT0006000.matrix' u ($1*0.0042+-0.0075):($2*0.0042+0.0025):3 matrix with pm3d
