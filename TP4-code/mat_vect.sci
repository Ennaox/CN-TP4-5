function [r]= my_sparse_matvect(AA,JA,IA,n,m,nnz_,v)
cpt = 1;
r = zeros(m,1);
for i = 1:m;
    for j = 1:IA(i+1)-IA(i)
        r(i) = r(i) + AA(cpt) * v(j);
        cpt = cpt + 1; 
    end
end
endfunction
