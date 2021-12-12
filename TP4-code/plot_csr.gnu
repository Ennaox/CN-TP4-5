set term png size 1900,1000

set output "image/CSR_time.png"

set grid

set ylabel "Temps"
set xlabel "Taille matrice nxn"

plot "data/csr.dat" using 3:1 t "CSR" w lp; 

set output "image/CSRt_time.png"
plot "data/csr.dat" using 3:2 t "CSRt" w lp;