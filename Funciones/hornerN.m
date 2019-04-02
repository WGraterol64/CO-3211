function y = hornerN(coef, x, t)
  [n,m] = size(x);
  [l,s] = size(coef);
  y = NaN*zeros(n,1);
  for k = 1: n
    y(k) = coef(l);
    for i = l-1 : -1: 1
      y(k) = y(k)*(x(k)-t(i))+coef(i);
    end
  endfor
end