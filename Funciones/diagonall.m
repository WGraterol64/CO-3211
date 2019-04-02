

function x = diagonal(a,b)
    [n,m] = size(a)
    if m ~= n 
        display('Error; Matrices no cuadradas')
        return
    end
    x = NaN*ones(n,1)
    for i = 1 :1: n
        if abs(a(i,i)) < epsilon 
            disp('Cero en la diagonal')
            return
        end
        x(i) = b(i) / a(i,i)
    end
end
