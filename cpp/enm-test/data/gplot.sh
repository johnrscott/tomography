file=$1
gnuplot << EOF
set term dumb
unset key
set xlabel 'Purity'
set ylabel 'Error distance'
set grid
plot for [i=2:5] '$file' using 1:(column(i))
set term png
set output 'test.png'
replot
EOF

display test.png
