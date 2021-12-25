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

function [x,iter,err]=richardson(A,b,e,maxiter,alpha)
    err = list();
    taille = size(A,1);

    x = zeros(taille,1);
    
    iter = 0;
    
    residual = norm(b - A * x)/norm(b);
    err($+1) = residual;
    while residual > e && iter < maxiter do
        
        x = x + alpha * (b - A * x)
        
        iter = iter + 1;
        residual = norm((b - A * x))/norm(b);
        err($+1) = residual;
    end
endfunction

iter_max = 50000;

[A,b] = gen_poisson1D(100,-5,5);

[x,iter,err] = richardson(A,b,1e-10,iter_max,1/2);
disp(err(length(err)));

[x,iter,err2] = richardson(A,b,1e-10,iter_max,1/3);
disp(err2(length(err2)));

[x,iter,err3] = richardson(A,b,1e-10,iter_max,1/4);
disp(err3(length(err3)));

[x,iter,err4] = richardson(A,b,1e-10,iter_max,1/5);
disp(err4(length(err4)));

[x,iter,err5] = richardson(A,b,1e-10,iter_max,1/6);
disp(err5(length(err5)));

[fic, mod] = mopen("../richardson1_2.dat", "w");
[fic2, mod] = mopen("../richardson1_3.dat", "w");
[fic3, mod] = mopen("../richardson1_4.dat", "w");
[fic4, mod] = mopen("../richardson1_5.dat", "w");
[fic5, mod] = mopen("../richardson1_6.dat", "w");

max_ = max(length(err),length(err2),length(err3),length(err4),length(err5));
for i = 1:max_
    if i<length(err)
        mfprintf(fic,"%d %.16lf\n",i,err(i));
    end
     if i<length(err2)
        mfprintf(fic2,"%d %.16lf\n",i,err2(i));
    end
     if i<length(err3)
        mfprintf(fic3,"%d %.16lf\n",i,err3(i));
    end
     if i<length(err4)
        mfprintf(fic4,"%d %.16lf\n",i,err4(i));
    end 
    if i<length(err5)
        mfprintf(fic5,"%d %.16lf\n",i,err5(i));
    end
end
mclose(fic);
mclose(fic2);
mclose(fic3);
mclose(fic4);
mclose(fic5);

