set term png size 1900,1000

set output "image/LDLT_time.png"

set grid

set ylabel "Temps"
set xlabel "Taille matrice nxn"
set logscale y
set logscale x

plot "data/LDLT.dat" using 6:5 t "my_LU3b" w lp, "data/LDLT.dat" using 6:4 t "my_LU1b" w lp,"data/LDLT.dat" using 6:3 t "LU" w lp, "data/LDLT.dat" using 6:2 t "LDLT" w lp;

unset logscale y
unset logscale x

set output "image/LDLT_solo_time.png"
plot "data/LDLT.dat" using 6:2 t "LDLT" w lp;

set output "image/LDLT_err.png"
set ylabel "Erreur"

plot "data/LDLT.dat" using 6:1 t "LDLT" w lp;
