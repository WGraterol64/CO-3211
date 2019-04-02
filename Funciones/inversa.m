function Ai= inversa(A)
  [m,n] = size(A);
  Ai = eye(n);
  for k = 1 : 1 : n-1
     for i = k+1 :1 :n
          alf = A(i,k)/A(k,k);
          for j = k:1:n
              A(i,j) = A(i,j) - alf*A(k,j) ;
          endfor
          for j = 1:1:n
              Ai(i,j) = Ai(i,j) - alf*Ai(k,j) ;
          endfor
     endfor
  endfor
  for k = n : -1 : 2
      for i = k-1 :-1 :1
          alf = A(i,k)/A(k,k);
          for j = n:-1:k
             A(i,j) = A(i,j) - alf*A(k,j) ;
          endfor
          for j = n:-1:1
              Ai(i,j) = Ai(i,j) - alf*Ai(k,j) ;
          endfor
      endfor
  endfor
  display(A)
  display(Ai)
  for i =1:1:n
    p = A(i,i)
      for j = 1: 1: n
        Ai(i,j) = Ai(i,j)/p;
      endfor
      A(i,i) = A(i,i)/p;
  endfor
endfunction
