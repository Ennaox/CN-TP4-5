exec CSR.sci
exec mat_vect.sci

ROW = 7;
COL = 5;
FILL = 0.35;
A = sprand(ROW,COL,FILL);
A = full(A);
v = rand(COL,1);
disp(v);
disp(A);
disp(A*v);
[AA,JA,IA,n,m,nnz_] = myCSR(A);
v = my_sparse_matvect(AA,JA,IA,n,m,nnz_,v);
disp(v);
