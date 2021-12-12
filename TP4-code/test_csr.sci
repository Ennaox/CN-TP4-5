exec CSR.sci
exec CSRt.sci

ROW = 4;
COL = 5;
FILL = 0.35;
A = sprand(ROW,COL,FILL);
A = full(A);
disp(A);


[AA,JA,IA,col,li,nnz_]= myCSR(A);
disp("AA");
disp(AA(:));
disp("JA");
disp(JA(:));
disp("IA");
disp(IA(:));


disp(A');
[AAt,JAt,IAt,colt,lit,nnz_t] = myCSRt(AA,JA,IA,col,li,nnz_);
disp("AAt");
disp(AAt(:));
disp("JAt");
disp(JAt(:));
disp("IAt");
disp(IAt(:));
