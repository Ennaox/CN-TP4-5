exec CSR.sci
exec mat_vect.sci

ROW = 5;
COL = 7;
FILL = 0.35;
A = sprand(ROW,COL,FILL);
A = full(A);
v = ones(COL,1);
disp(A*v)
[AA,JA,IA,n,m,nnz_] = myCSR(A);
v = my_sparse_matvect(AA,JA,IA,n,m,nnz_,v);
disp(v);
