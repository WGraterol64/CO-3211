function A = Hilbert(n)
  A = NaN*ones(n,n);
  for i = 1 : 1 : n
    for j = 1 : 1 : n 
      A(i,j) = 1/(i+j-1);
    endfor
  endfor
endfunction