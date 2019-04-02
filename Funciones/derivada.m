function dp = derivada(a,x,t)
  [n,m] = size(a)
  [l, s] = size(x)
  dp = NaN*zeros(l,1);
  for i = 1:l
  dp(i) = 0; 
  p=a(n);
    for k=n-1:-1:1
      dp(i) = dp(i) * (x(i)-t(k));
      dp(i) = dp(i) + p;
      p = p * (x(i)-t(k)); 
      p = p + a(k);
    endfor
  endfor
endfunction
