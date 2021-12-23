function [A,b] = gen_poisson1D(n,T0,T1)
    A = zeros(n, n);
    for i = 1:n
        A(i, i) = 2;
        
    end

    for i = 1:n-1
        A(i + 1, i) = -1;
        A(i, i + 1) = -1;
    end
    
    b = zeros(n, 1);
    b(1) = T0;
    b(n) = T1;
endfunction

function [x,iter,err]=jacobi(A,b,e,maxiter)
    err = list();
    taille = size(A,1);
    
    D = diag(diag(A));
    invD = inv(D);
    
    x = zeros(taille,1);
    
    iter = 0;
    
    residual = norm(b - A * x);
    err($+1) = residual;
    while residual > e && iter < maxiter do
        
        x = x + invD * (b - A * x)
        
        iter = iter + 1;
        residual = norm((b - A * x));
        err($+1) = residual;
    end
endfunction

[A,b] = gen_poisson1D(100,-5,5);

[x,iter,err] = jacobi(A,b,1e-10,100000);

disp(norm(A\b - x)/norm(A));

[fic, mod] = mopen("../jacobi.dat", "w");
for i = 1:length(err)
    mfprintf(fic,"%.16lf %d\n",err(i),i);
end
mclose(fic);

