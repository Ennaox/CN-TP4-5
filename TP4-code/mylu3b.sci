function [L,U] = mylu3b(A)
n = size(A,1);
for k = 1 : n - 1
    for i = k + 1 : n
        A(i, k) = A(i, k)/A(k, k);
    end
    for i = k + 1 : n
        for j = k + 1 : n
            A(i, j) = A(i, j) - A(i, k) * A(k, j);
        end
    end
end
U = triu(A);
L = tril(A);
for k = 1 : n
    L(k,k) = 1;
end
endfunction
