
function x = matrizTriangularL(a,b)
    [m,n] = size(a);
    if m ~= n 
        display('Error; Matrices no cuadradas')
        return
    end
    x = NaN*ones(n,1);
    x(1)= b(1)/a(1, 1);
    for i = 2 : 1 : n
        %if abs(a(n)) < epsilon
         %   disp('Cero en la diagonal')
          %  return
        %end
        
        suma = 0.0;
        
        for j = 1 : 1 :i-1
            suma = suma + a(i,j)*x(j);
        end
        
        x(i) = (b(i) - suma)/a(i,i);
    end
endfunction
    
