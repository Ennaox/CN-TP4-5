function [AA,JA,IA,n,m,nnz_]=myCSR(A)
n = size(A,2);
m = size(A,1);
AA = list();
JA = list();
IA = list(1);

for i = 1:m
    nb_elem = 0;
    for j = 1:n
        if A(i,j) ~= 0
            AA($+1) = A(i,j);
            JA($+1) = j;
            nb_elem = nb_elem + 1;
        end
    end
    IA($+1) = IA($) + nb_elem;
end
nnz_ = IA($) - 1;

endfunction
