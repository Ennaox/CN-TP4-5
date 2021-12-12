exec LDLT.sci
exec mylu3b.sci
exec mylu1b.sci

TAILLE_MAX = 300;
NB_REP = 5;

[fic, mod] = mopen("data/LDLT.dat", "w");
for taille = 10:10:TAILLE_MAX
    dif = 0;
    disp(string(taille)+"/"+string(TAILLE_MAX));
    time_ldlt = 0;
    time_lu = 0;
    time_lu1b = 0;
    time_lu3b = 0;
    for rep = 1 : NB_REP
        B = tril(rand(taille,taille)+ones(taille,taille));
        A = zeros(taille,taille);
        A = B + B';
        Abis = A;
        tic();
        A = LDLT(A);
        time_ldlt = time_ldlt + toc();
        tic()
        [L1,U1] = lu(Abis);
        time_lu = time_lu + toc();
        tic()
        [L2,U2] = mylu3b(Abis);
        time_lu3b = time_lu3b + toc();
        tic()
        [L3,U3] = mylu1b(Abis);
        time_lu1b = time_lu1b + toc();
        
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
        dif = dif + norm(Abis - newA);
        
    end
    mfprintf(fic,"%.17lf %.14lf %.14lf %.14lf %.14lf %d\n",dif/NB_REP,time_ldlt/NB_REP,time_lu/NB_REP,time_lu1b/NB_REP,time_lu3b/NB_REP,taille);
end
mclose(fic);
