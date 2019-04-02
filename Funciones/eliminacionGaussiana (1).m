
function [a, b]= eliminacionGaussiana(a,b, p)
    % Sin Pivoteo
    if p == 0
        [m,n] = size(a);
        for k = 1 : 1 : n-1
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j) ;
                end
                b(i) = b(i) - alf*b(k);
            end
        end
    % Con pivoteo
    if p == 1
        [m,n] = size(a);
        for k = 1 : 1 : n-1
		max = abs(a(k,k))
        imax = k
	    for i = k+1 :1:n
		    if abs(a(i,k))>maximo
                maximo = abs(a(i,k))
                imax = i
            end
        end
        for j = k :1:n
            temp = a(k,j)
            a(k,j) = a(imax,j)
            a(imax,j) = temp
        end
        temp = b(k)
        b(k) = b(imax)
        b(imax) = temp
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j);
                end
                b(i) = b(i) - alf*b(k);
            end
        end
     end
 end
