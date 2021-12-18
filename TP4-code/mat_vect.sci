function [r]= my_sparse_matvect(AA,JA,IA,li,col,nnz_,v)
cpt = 1;
r = zeros(col,1);
for i = 1:col;
    for j = 1:IA(i+1)-IA(i)
        r(i) = r(i) + AA(cpt) * v(JA(cpt));
        cpt = cpt + 1; 
    end
end
endfunction
