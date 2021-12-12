exec CSR.sci
exec CSRt.sci

TAILLE_MAX = 250;
NB_REP = 10;
FILL = 0.35;

[fic, mod] = mopen("data/csr.dat", "w");
for TAILLE=10:10:TAILLE_MAX
    temps_csr = 0;
    temps_csrt = 0;
    disp(string(TAILLE)+"/"+string(TAILLE_MAX));
    for REP = 1:NB_REP;
        A = sprand(TAILLE,TAILLE,FILL);
        A = full(A);
        tic();
        [AA,JA,IA,col,li,nnz_]= myCSR(A);
        temps_csr = temps_csr + toc();
        tic();
        [AAt,JAt,IAt,colt,lit,nnz_t] = myCSRt(AA,JA,IA,col,li,nnz_);
        temps_csrt = temps_csrt + toc();
    end
    mfprintf(fic,"%.14lf %.14lf %d\n",temps_csr/NB_REP, temps_csrt/NB_REP, TAILLE);
end
mclose(fic);
