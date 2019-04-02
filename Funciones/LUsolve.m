function x = LUsolve(L,U,b)
  [m,n] = size(U);
  % Para resolver LUx=b se resuelve primero Ly=b y luego Ux=y
  y = NaN*ones(n,1);
  x = NaN*ones(n,1);
  % Resolver Ly=b por sust. hacia adelante
  y(1) = b(1)/L(1,1);
  for i = 2:n
    s = 0;
    for j=1:i-1
      s = s+L(i,j)*y(j);
    endfor
    y(i) = (b(i)-s)/L(i,i);
  endfor
  % Resolver Ux=y por sust. hacia atras
  x(n) = y(n)/U(n,n);
  for i = n-1:-1:1
    s = 0;
    for j=i+1:n
      s = s+U(i,j)*x(j);
    endfor
    x(i) = (y(i)-s)/U(i,i);
  endfor
endfunction
