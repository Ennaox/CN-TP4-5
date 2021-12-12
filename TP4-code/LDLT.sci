function [A] = LDLT(A)
li = size(A,1);
v = matrix(li,1);
for j = 1 : li
    for i = 1:j-1
        v(i) = A(j,i) * A (i,i);
    end
    if j>1
        A(j,j) = A(j,j)-A(j,1:j-1)*v(1:j-1);
        A(j+1:li,j) = (A(j+1:li,j) - A(j+1:li,1:j-1) * v(1:j-1))*1/A(j,j);
    else    
        A(j+1:li,j) = A(j+1:li,j)*1/A(j,j);
    end
end
endfunction
