function [r, y, k] = met_inv(A, x)
  % Hallar LU
  [m, n] = size(A);
  T = LUdecomp(A);
  L = tril(T,-1) + eye(n);
  U = triu(T);
  for k = 1:1000
    % Resolver el sistema para la inversa
    y = LUsolve(L, U, x);
    % Actualizar el autovalor
    r = y(1)/x(1);
    % Revisar si parar
    if norm(x-(y/norm(y,Inf)),Inf) < 10e-6
      % Retornar el cambio para regresar el menor autovalor
      r = 1/r;
      return
    endif
    % Hacer el vector viejo el vector nuevo normalizao 
    x = y/norm(y,Inf);
  endfor
  disp("numero maximo de iteraciones alcanzadas")
  r = 1/r;
endfunction
