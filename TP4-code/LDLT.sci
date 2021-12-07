function [A] = LDLT(A)
n = size(A,1);
v = matrix(n,1);
for j = 1 : n
    for i = 1:j-1
        v(i) = A(j,i) * A (i,i);
    end
    if j>1
        A(j,j) = A(j,j)-A(j,1:j-1)*v(1:j-1);
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1) * v(1:j-1))*1/A(j,j);
    else    
        A(j+1:n,j) = A(j+1:n,j)*1/A(j,j);
    end
end
endfunction
