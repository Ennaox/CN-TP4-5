function [AA,JA,IA,li,col,nnz_]=myCSR(A)
col = size(A,2);
li = size(A,1);
AA = list();
JA = list();
IA = list(1);

for i = 1:li
    nb_elem = 0;
    for j = 1:col
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
