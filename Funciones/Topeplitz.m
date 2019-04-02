function T = Topeplitz(v)
  [m, n] = size(v)
  % v(0) = v'(n+1)
  k = floor(n/2)+1
  T = NaN*ones(k,k);
  for i = 1:1:k
    a(i)=v(k-1+i);
  endfor
  for i = 1:1:k-1
    b(i)=v(i);
  endfor
  for i = 1:1:k
    for j= 1:1:k
      if j<i
        T(i,i-j) = v(k-j);
      else
        T(i,j) = v(k+j-i);
      endif
    endfor
  endfor
endfunction
