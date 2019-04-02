function A = LUdecomp(A)
    [n,m] = size(A);
    for k = 2:1:n
        % Calculo de L: Guardar los multiplicadores
        for i = k:1:n  
             A(i, k-1) = A(i, k-1)/A(k-1, k-1);
            for j = k:1:n
                % Calculo de U: Hacer la eliminacion sobre U
                A(i, j) = A(i, j) - A(i, k-1)*A(k-1, j);
            endfor
        endfor
    endfor
endfunction