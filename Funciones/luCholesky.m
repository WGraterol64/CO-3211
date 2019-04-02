
function l = luCholesky(A)
  [m, n] = size(A);
  l = zeros(n,n);
  for k = 1:n
    s = 0;
    for i = 1:k-1
        s = s+l(k,i)*l(k,i);
    endfor
    l(k,k) = sqrt(A(k,k)-s);
    for i = k+1:n
      s = 0;
      for j = 1:k-1
        s = s + l(i,j)*l(k,j);
      endfor
      l(i,k) = (A(i,k)-s)/l(k,k);
    endfor
  endfor
endfunction
