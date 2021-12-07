exec CSR.sci
exec CSRt.sci

ROW = 7;
COL = 5;
FILL = 0.35;
A = sprand(ROW,COL,FILL);
A = full(A);
disp(A);


[AA,JA,IA,n,m,nnz_]= myCSR(A);
disp("AA");
disp(AA(:));
disp("JA");
disp(JA(:));
disp("IA");
disp(IA(:));


disp(A');
[AAt,JAt,IAt,nt,mt,nnz_t] = myCSRt(AA,JA,IA,n,m,nnz_);
disp("AAt");
disp(AAt(:));
disp("JAt");
disp(JAt(:));
disp("IAt");
disp(IAt(:));