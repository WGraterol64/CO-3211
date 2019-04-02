function [xn, k, totaltime] = jacobi(A, b, xb, maxiter, tol)
  n = size(xb);
  xn = NaN*ones(n,1);
  initial = cputime;
  for k = 0:1:maxiter
    display(cputime)
    for i = 1:1:n
      s1 = 0;
      for j = 1:1:i-1
        s1 = s1+A(i,j)*xb(j);
      endfor
      s2 = 0;
      for j = i+1:1:n
        s2 = s2+A(i,j)*xb(j);
      endfor
      xn(i) = (b(i)-s1-s2)/A(i,i);
    endfor
    if norm(xn-xb, Inf) < tol
      display("Solucion Hallada! Numero de iteraciones: ")
      disp(k)
      endtime = cputime
      totaltime = endtime-initial
      display("Tiempo de ejecucion: ")
      disp(totaltime)
      return
    endif
    xb = xn;
  endfor
  display("Numero maximo de iteraciones alcanzadas :(")
  endtime = cputime;
  totaltime = endtime-initial;
  display("Tiempo de ejecucion: ")
  disp(totaltime)
endfunction
