
function x = matrizTriangular(a,b)
    [m,n] = size(a)
    if m ~= n 
        display('Error; Matrices no cuadradas')
        return
    end
    x = NaN*ones(n,1);
    x(n)= b(n)/a(n, n);
    for i = n-1 : -1 : 1
        %if abs(a(n)) < epsilon
         %   disp('Cero en la diagonal')
          %  return
        %end
        
        suma = 0.0;
        
        for j = i+1 : 1 :n
            suma = suma + a(i,j)*x(j);
        end
        
        x(i) = (b(i) - suma)/a(i,i);
    end
end
    
