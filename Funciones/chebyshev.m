function c = chebyshev(m,a,b)
  c = NaN*ones(m,1)
  for i =1:m
    c(i) = 0.5*(a+b) + 0.5 * (b-a) * cos(((2*i-1)/(2*m))*pi)
  endfor
endfunction
