function [r, y, k] = met_dir_desp(A, x, delta)
  [m, n] = size(A);
  A = A-delta*eye((n)); 
  for k = 1:1000
    % Resolver el sistema
    y = A*x;
    % Actualizar el autovalor
    r = y(1)/x(1);
    % Revisar si parar
    if norm(x-(y/norm(y,Inf)),Inf) < 10e-6
      % Retornar el cambio para regresar el menor autovalor
      r = r+delta;
      return
    endif
    % Hacer el vector viejo el vector nuevo normalizado
    x = y/norm(y,Inf);
  endfor
  disp("numero maximo de iteraciones alcanzadas")
  r = r+delta;
endfunction
