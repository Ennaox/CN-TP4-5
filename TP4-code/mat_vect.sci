function [v]= my_sparse_matvect(AA,JA,IA,n,m,nnz_,v)
cpt = 1;

for i = 1:n;
    vi = v(i);
    disp(IA(i));
    disp(IA(i+1));
    disp("-----------");
    for j = 1:IA(i+1)-IA(i)
        v(i) = v(i) + AA(cpt) * vi;
        disp(cpt);
        cpt = cpt + 1; 
    end
    disp("-----------");
end
endfunction
