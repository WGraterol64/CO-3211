function [a, b]= eliminacionGaussiana(a,b, p)
    if p == 0
        [m,n] = size(a);
        for k = 1 : 1 : n-1
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j) ;
                endfor
                b(i) = b(i) - alf*b(k);
            endfor
        endfor
    elseif p == 1
        [m,n] = size(a);
        for k = 1 : 1 : n-1
            for i =k:1:n
              if abs(a(k,k)) < abs(a(i,k))
                for j = k: 1: n
                  t=a(k,j);
                  a(k,j) = a(i,j);
                  a(i,j) = t;
                endfor
                t = b(k);
                b(k) = b(i);
                b(i) = t;
              endif
            endfor
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j);
                endfor
                b(i) = b(i) - alf*b(k);
            endfor
        endfor
    elseif p != 0 & p !=1
       display("Error, definir pivoteo")
       return
    endif
endfunction
 