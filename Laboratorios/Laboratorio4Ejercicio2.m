% Ejercicio 2:

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

% Decomp de LU normal:
function [L,U] = LUdecomp(A)
    [n,m] = size(A);
    for k = 2:1:n
        % Calculo de L: Guardar los multiplicadores
        for i = k:1:n  
             A(i, k-1) = A(i, k-1)/A(k-1, k-1);
            for j = k:1:n
                % Calculo de U: Hacer la eliminacion sobre U
                A(i, j) = A(i, j) - A(i, k-1)*A(k-1, j);
            endfor
        endfor
    endfor
    L = tril(A,-1) + eye(n)
    U = triu(A)
endfunction

% Resolucion de sistemas LUx=b
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

% Matrices y vectores de cada sistema de ecuaciones:

A = zeros(3,3)
ba = NaN
B = zeros(4,4)
bb = NaN

% Guardamos la primera matriz y el primer vector:

A(1,1) = 1.012
A(1,2) = -2.132
A(1,3) = 3.1041
A(2,1) = -2.132
A(2,2) = 4.096
A(2,3) = -7.013
A(3,1) = 3.1040
A(3,2) = -7.013
A(3,3) = 0.014
ba = [1.984; -5.049; -3.895]

% Guardamos la segunda matriz con el segundo vector:

B(1,1) = 6
B(1,2) = 2
B(1,3) = 1
B(1,4) = -1
B(2,1) = 2
B(2,2) = 4
B(2,3) = 1
B(2,4) = 0
B(3,1) = 1
B(3,2) = 1
B(3,3) = 4
B(3,4) = -1
B(4,1) = -1
B(4,2) = 0
B(4,3) = -1
B(4,4) = 3
bb = [0; 7; -1; -2]

% Parte 1:

% Sobre a primera matriz no podemos aplicar Cholesky pues ella tiene autovalores 
% negativos: 
eig(A)
% ans =
%   10.687475
%   -0.060824
%   -5.504651

% y se tiene que toda matriz simetrica definida positiva tiene autovalores 
% positivos ( Extraido de: https://es.wikipedia.org/wiki/Matriz_definida_positiva)

% De igual manera, B es simetrica y tiene autovalores positivos:
eig(B)

% ans =
%   1.9261
%   3.2198
%   3.8458
%   8.0083

% Por lo tanto B admite Cholesky.

% Parte 2:
% Aplicamos luCholesky a B

luCholesky(B)
% ans =

%   2.44949   0.00000   0.00000   0.00000
%   0.81650   1.82574   0.00000   0.00000
%   0.40825   0.36515   1.92354   0.00000
%  -0.40825   0.18257  -0.46789   1.60657

% Parte 3:
% Aplicamos LUdecomp a A
LUdecomp(A)

% L =

%   1.00000   0.00000   0.00000
%  -2.10672   1.00000   0.00000
%   3.06719   1.19776   1.00000

% U =

%   1.01200  -2.13200   3.10410
%   0.00000  -0.39553  -0.47353
%   0.00000   0.00000  -8.93970


% Parte 3:

% Para A

[L, U] = LUdecomp(A)
LUsolve(L,U,ba)
% ans =
%   1.00137
%   1.00061
%   0.99994

% Para B:

L = luCholesky(B)
LUsolve(L,L',bb)
% ans =
%  -0.85864
%   2.41885
%  -0.95812
%  -1.27225

% Pregunta 5: 

% Si pues es un teorema que una matriz no singular (Ellas son no singulares)
% que puede ser factorizada como LL' tiene que ser definida positiva.

% Entonces como la matriz B cumple con las hipotesis se tiene que B debe ser
% definida positiva