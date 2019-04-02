function [r, y, k] = met_dir(A, x)
  for k = 1:1000
    % Aproximar al nuevo autovector
    y = A*x;
    % Aproximar al nuevo autovalor
    r = y(1)/x(1);
    % Revisar si parar
    if norm(x-(y/norm(y,Inf)),Inf) < 10e-6
      return
    endif
    % Hacer el vector viejo el vector nuevonormalizado 
    x = y/norm(y,Inf);
  endfor
  disp("numero maximo de iteraciones alcanzadas")
endfunction
