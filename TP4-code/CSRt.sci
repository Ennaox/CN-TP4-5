function [AAt,JAt,IAt,colt,lit,nnz_]=myCSRt(AA,JA,IA,col,li,nnz_)
AAt = list();
JAt = list();
IAt = list(1);
l = 0;
colt = li;
lit = col;
for i = 1:li
    cpt = 0;
    for  j = 1:length(AA)
        if JA(j) == i
            AAt($+1) = AA(j); 
            l = l + 1;
            cpt = cpt + 1;
            for k = 1 : length(IA) - 1
                if j >= IA(k) & j <= IA(k+1)
                    JAt(l) = k;
                end
            end
        end
    end
    IAt($+1) = IAt($) + cpt;
end
endfunction
