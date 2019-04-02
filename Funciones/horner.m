function y = horner(coef, x)
  [n,m] = size(x);
  [l,t] = size(coef);
  y = NaN*ones(n,1);
  for k = 1:n
    y(k) = coef(1, 1);
    for i = 2 : l
      y(k) = y(k)*x(k)+coef(i, 1);
    endfor
  endfor
endfunction
