% Pregunta 1:

% Parte a:

% Si aplicamos las formulas dadas anteriormente nos quedara el siguiente sistema
% matricial At=b:

A = [[4,-1,-1,0,0,0,0,0];[-1,4,0,-1,0,0,0,0];[-1,0,4,-1,-1,0,0,0];[0,-1,-1,4,0,-1,0,0];
     [0,0,-1,0,4,-1,-1,0];[0,0,0,-1,-1,4,0,-1];[0,0,0,0,-1,0,4,-1];[0,0,0,0,0,-1,-1,4]]
b = [5;15;0;10;0;10;5;15]

% A =
%
%   4  -1  -1   0   0   0   0   0
%  -1   4   0  -1   0   0   0   0
%  -1   0   4  -1  -1   0   0   0
%   0  -1  -1   4   0  -1   0   0
%   0   0  -1   0   4  -1  -1   0
%   0   0   0  -1  -1   4   0  -1
%   0   0   0   0  -1   0   4  -1
%   0   0   0   0   0  -1  -1   4

% Parte 2:

% A es pentadiagonal pues las unicas entradas distintas de 0 estan sobre las dos
% diagonales superiores y las dos inferiores a la diagonal principal.

% La matriz es simetrica pues ella es igual a su traspuesta

% Es definida positivo pues es simetrica y todos sus autovalores son positivos 
% (Ver parte 1 del ejercicio2)

eig(A)
% ans =
%   1.3820
%   2.3820
%   3.3820
%   3.6180
%   4.3820
%   4.6180
%   5.6180
%   6.6180

% Parte 3:

% Descomp. Cholesky

function l = luCholesky(A)
  [m, n] = size(A);
  l = zeros(n,n);
  for k = 1:n
    s = 0;
    for i = 1:k-1
        s = s+l(k,i)*l(k,i);
    endfor
    l(k,k) = sqrt(A(k,k)-s);
    for i = k+1:n
      s = 0;
      for j = 1:k-1
        s = s + l(i,j)*l(k,j);
      endfor
      l(i,k) = (A(i,k)-s)/l(k,k);
    endfor
  endfor
endfunction

% Ejecucion sobre A

L = luCholesky(A)

% L =
%
%   2.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
%  -0.50000   1.93649   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
%  -0.50000  -0.12910   1.93218   0.00000   0.00000   0.00000   0.00000   0.00000
%   0.00000  -0.51640  -0.55205   1.85164   0.00000   0.00000   0.00000   0.00000
%   0.00000   0.00000  -0.51755  -0.15430   1.92570   0.00000   0.00000   0.00000
%   0.00000   0.00000   0.00000  -0.54006  -0.56257   1.84170   0.00000   0.00000
%   0.00000   0.00000   0.00000   0.00000  -0.51929  -0.15862   1.92488   0.00000
%   0.00000   0.00000   0.00000   0.00000   0.00000  -0.54298  -0.56426   1.84032

LT = luCholesky(A)'
% LT =

%   2.00000  -0.50000  -0.50000   0.00000   0.00000   0.00000   0.00000   0.00000
%   0.00000   1.93649  -0.12910  -0.51640   0.00000   0.00000   0.00000   0.00000
%   0.00000   0.00000   1.93218  -0.55205  -0.51755   0.00000   0.00000   0.00000
%   0.00000   0.00000   0.00000   1.85164  -0.15430  -0.54006   0.00000   0.00000
%   0.00000   0.00000   0.00000   0.00000   1.92570  -0.56257  -0.51929   0.00000
%   0.00000   0.00000   0.00000   0.00000   0.00000   1.84170  -0.15862  -0.54298
%   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   1.92488  -0.56426
%   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   1.84032

% Estas matrices obtenidas son la descomp. LL' de A

% Parte 4:

% Se tiene que L es de banda puesto que para todo j>i y i+2>j, L(i,j) = 0
% Se tiene que L es de banda puesto que para todo j<i y i+2<j, L(i,j) = 0

% Parte 5:

% Funcion para resolver el sistema

function x = LUsolve(L,U,b)
  [m,n] = size(U);
  % Para resolver LUx=b se resuelve primero Ly=b y luego Ux=y
  y = NaN*ones(n,1);
  x = NaN*ones(n,1);
  % Resolver Ly=b por sust. hacia adelante
  y(1) = b(1)/L(1,1);
  for i = 2:n
    s = 0;
    for j=1:i-1
      s = s+L(i,j)*y(j);
    endfor
    y(i) = (b(i)-s)/L(i,i);
  endfor
  % Resolver Ux=y por sust. hacia atras
  x(n) = y(n)/U(n,n);
  for i = n-1:-1:1
    s = 0;
    for j=i+1:n
      s = s+U(i,j)*x(j);
    endfor
    x(i) = (y(i)-s)/U(i,i);
  endfor
endfunction

LUsolve(L,LT,b)
% ans =
%   3.6842
%   6.3158
%   3.4211
%   6.5789
%   3.4211
%   6.5789
%   3.6842
%   6.3158