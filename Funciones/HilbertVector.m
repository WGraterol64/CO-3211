function b = HilbertVector(n)
  A = Hilbert(n);
  b = A*ones(n,1);
endfunction