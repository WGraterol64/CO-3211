function x = tempplaca(Ts, Tiz, Td, Tin)
  A = [[4,-1,-1,0,0,0,0,0];[-1,4,0,-1,0,0,0,0];[-1,0,4,-1,-1,0,0,0];
           [0,-1,-1,4,0,-1,0,0];[0,0,-1,0,4,-1,-1,0];[0,0,0,-1,-1,4,0,-1];
           [0,0,0,0,-1,0,4,-1];[0,0,0,0,0,-1,-1,4]]
  b = [Ts+Tiz;Tiz+Tin;Ts;Tin;Ts;Tin;Ts+Td;Td+Tin]
  L = luCholesky(A);
  x = LUsolve(L, L',b);
endfunction
