function s = neville(e, x, y)
  [n, m] = size(y);
  [l, t] = size(e);
  s = NaN*ones(t,1);
  q = NaN*zeros(m,m);
  for i = 1:m
    q(i, 1) = y(i) ;
  endfor
  for k = 1:t
    for i = 2:m
      for j = 2:i
        q(i,j) = (((e(k)-x(i-j+1)).*q(i, j-1))-((e(k)-x(i)).*q(i-1,j-1)))./(x(i)-x(i-j+1));    
      endfor
    endfor
    s(k) = q(m,m);
  endfor
endfunction
