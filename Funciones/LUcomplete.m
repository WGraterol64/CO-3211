function [x, totaltime] = LUcomplete(A,b)
  initial = cputime;
  [m, n] = size(A)
  T = LUdecomp(A);
  L = tril(T,-1) + eye(n)
  U = triu(T)
  x = LUsolve(L, U, b);
  endtime = cputime;
  totaltime = initial - endtime;
  display("Tiempo de ejecucion: ")
  disp(totaltime)
endfunction
