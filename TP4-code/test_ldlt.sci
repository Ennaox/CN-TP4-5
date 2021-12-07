exec LDLT.sci

TAILLE_MAX = 751;
NB_REP = 10;

[fic, mod] = mopen("data/LDLT.dat", "w");
for taille = 1:10:TAILLE_MAX
    dif = 0;
    disp(string(taille)+"/"+string(TAILLE_MAX));
    time_ldlt = 0;
    time_lu = 0;
    for rep = 1 : NB_REP
        B = tril(rand(taille,taille)+ones(taille,taille));
        A = zeros(taille,taille);
        A = B + B';
        Abis = A;
        tic();
        A = LDLT(A);
        time_ldlt =+ toc();
        tic()
        [L1,U1] = lu(Abis);
        time_lu =+ toc();
        L = tril(A);
        for k = 1 : taille
            L(k,k) = 1;
        end
        D = diag(A);
        realD = zeros(taille,taille);
        for i = 1 : taille
            realD(i,i) = D(i);
        end 
        newA = (L*realD*L');
        dif =+ norm(Abis - newA);
        
    end
    mfprintf(fic,"%.17lf %.14lf %.14lf %d\n",dif/NB_REP,time_ldlt/NB_REP,time_lu/NB_REP,taille);
end
mclose(fic);
