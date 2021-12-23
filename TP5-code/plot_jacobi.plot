set term png size 1900,1000

set output "image/jacobi_converg.png"

set grid

set ylabel "Itération"
set xlabel "Erreur"
set logscale y 3

plot "jacobi.dat" using 2:1 t "Convergence Jacobi" w lp;