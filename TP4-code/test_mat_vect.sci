exec CSR.sci
exec mat_vect.sci

NB_REP = 10;
MAX_SIZE = 505;
FILL = 0.35;
[fic, mod] = mopen("data/mat_vect.dat", "w");

for SIZE=5:10:MAX_SIZE
    error_ = 0;
    error_2 = 0;
    disp(string(SIZE)+"/"+string(MAX_SIZE));
    for REP=1:NB_REP
        A = sprand(SIZE,SIZE,FILL);
        A = full(A);
        v = ones(SIZE,1);
        [AA,JA,IA,li,col,nnz_] = myCSR(A);
        r = my_sparse_matvect(AA,JA,IA,li,col,nnz_,v);
        error_ = error_ + norm(r - (A*v));
        v = sprand(SIZE,1,0.35);
        v = full(v);
        r = my_sparse_matvect(AA,JA,IA,li,col,nnz_,v);
        error_2 = error_2 + norm(r - (A*v));
    end
    mfprintf(fic,"%.17lf %.17lf %d\n",error_/NB_REP, error_2/NB_REP, SIZE);
end
mclose(fic);
