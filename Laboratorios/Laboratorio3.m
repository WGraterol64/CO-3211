% CO3211-Laboratorio 3:
% Wilfredo Graterol
% 15-10639

% Pregunta 1:
function [a, b]= eliminacionGaussiana(a,b, p)
    if p == 0
        [m,n] = size(a);
        for k = 1 : 1 : n-1
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j) ;
                endfor
                b(i) = b(i) - alf*b(k);
            endfor
        endfor
    elseif p == 1
        [m,n] = size(a);
        for k = 1 : 1 : n-1
            for i =k:1:n
              if abs(a(k,k)) < abs(a(i,k))
                for j = k: 1: n
                  t=a(k,j);
                  a(k,j) = a(i,j);
                  a(i,j) = t;
                endfor
                t = b(k);
                b(k) = b(i);
                b(i) = t;
              endif
            endfor
            for i = k+1 :1 :n
                alf = a(i,k)/a(k,k);
                for j = k:1:n
                    a(i,j) = a(i,j) - alf*a(k,j);
                endfor
                b(i) = b(i) - alf*b(k);
            endfor
        endfor
    elseif p != 0 & p !=1
       display("Error, definir pivoteo")
       return
    endif
endfunction
 
 
function A = LUdecomp(A)
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
endfunction

function x = sustAdel(a,b)
    [m,n] = size(a);
    if m ~= n 
        display('Error; Matrices no cuadradas')
        return
    end
    x = NaN*ones(n,1);
    x(1)= b(1)/a(1, 1);
    for i = 2 : 1 : n
        suma = 0.0;
        
        for j = 1 : 1 :i-1
            suma = suma + a(i,j)*x(j);
        end
        
        x(i) = (b(i) - suma)/a(i,i);
    end
endfunction
    

function x = sustAtras(a,b)
    [m,n] = size(a)
    if m ~= n 
        display('Error; Matrices no cuadradas')
        return
    end
    x = NaN*ones(n,1);
    x(n)= b(n)/a(n, n);
    for i = n-1 : -1 : 1
        suma = 0.0;
        
        for j = i+1 : 1 :n
            suma = suma + a(i,j)*x(j);
        end
        
        x(i) = (b(i) - suma)/a(i,i);
    end
end

% Pregunta 2:

% Matriz de Toeplitz
function T = Toeplitz(v)
  [m, n] = size(v)
  % v(0) = v'(n+1)
  k = floor(n/2)+1
  T = NaN*ones(k,k);
  for i = 1:1:k
    a(i)=v(k-1+i);
  endfor
  for i = 1:1:k-1
    b(i)=v(i);
  endfor
  for i = 1:1:k
    for j= 1:1:k
      if j<i
        T(i,i-j) = v(k-j);
      else
        T(i,j) = v(k+j-i);
      endif
    endfor
  endfor
endfunction

% Parte (a):
% Vectores de prueba

v0 = [1;1/2;1/3;1/5;1/7]
v1 = ones(25,1);
v1(1) = 1;
v1(2) = 1/2;

for i = 1:1:23
  v1(i+2) = 1./(2*i+1);
endfor

v2 = [-3;-2;-1;0.01; 1; 2; 3];

 t0 = Topeplitz(v0')
 t1 = Topeplitz(v1')
 t2 = Topeplitz(v2')
 
% Calculo de autovalores

% eig(t0) = [1.007938; -0.043187; 0.035249]

% eig(t1) = [6.9970e-01;-1.4687e-01; 1.2934e-02; -3.9639e-04; -1.3814e-04;
%           -2.5390e-05; 1.7640e-05; -5.7048e-07;  1.4360e-08; -3.9429e-10;
%            8.9616e-12; -1.4040e-13; 9.5901e-16]

% eig(t2) = [0.01000 + 4.47214i; 0.01000 - 4.47214i; 0.01000 + 0.00000i; 
%            0.01000 + 0.00000i]

% Calculo de determinantes

% det(t0) = -0.0015344
% det(t1) = -1.0659e-82
% det(t2) =  0.0020000

% Calculo de normas

% norm(t0,Inf) = 1.8333
% norm(t1,Inf) =  2.7244
% norm(t2,Inf) =  6.0100

% Parte (b)

% Vectores solucion
b0 = t0*[1;1;1]
b1 = t1*[1;1;1;1;1;1;1;1;1;1;1;1;1]
b2 = t2*[1;1;1;1] 
% b0 = [0.67619; 1.03333; 1.83333]
% b1 = [0.38994; 0.41629; 0.44670; 0.48226; 0.52454; 0.57582; 0.63970;
%       0.72224; 0.83480; 1.00254; 1.30139; 1.76435; 2.72435]
% b2 = [6.0100; 2.0100; -1.9900; -5.9900]

% Eliminacion Gaussiana sin pivoteo

 [A,B] = eliminacionGaussiana(t0,b0,0)
 
% A =
%
%   0.33333   0.20000   0.14286
%   0.00000   0.03333  -0.01429
%   0.00000   0.00000  -0.13810
%
% B =
%
%   0.676190
%   0.019048
%  -0.138095

sustAtras(A,B)

%ans =

%   1.00000
%   1.00000
%   1.00000

[A, B] = eliminacionGaussiana(t1,b1,0)
sustAtras(A,B)
% ans =

%   1.000000
%   1.000000
%   0.999990
%   1.000393
%   0.994237
%   1.043331
%   0.808201
%   1.534312
%   0.036283
%   2.124064
%   0.181483
%   1.338305
%   0.939402

[A, B] = eliminacionGaussiana(t2,b2,0)
% A =
%
% 0.01000     1.00000     2.00000     3.00000
% 0.00000   100.01000   201.00000   302.00000
% 0.00000     0.00000     0.06000     0.08009
% 0.00000     0.00000     0.00000     0.03333

% B =

%     6.010000
%   603.010000
%     0.140087
%     0.033333

sustAtras(A,B)
% ans =
%   1.00000
%   1.00000
%   1.00000
%   1.00000

% Eliminacion Gaussiana con pivoteo
[A, B] = eliminacionGaussiana(t0,b0,1)
% A =
%
%   1.00000   0.50000   0.33333
%   0.00000   0.08333   0.03333
%   0.00000   0.00000   0.01841

% B =

%   1.833333
%   0.116667
%   0.018413

sustAtras(A,B)
% ans =

%  1.00000
%  1.00000
%  1.00000


[A, B] = eliminacionGaussiana(t1,b1,1)
sustAtras(A,B)
% ans =

%   1.00000
%  1.00000
%  0.99999
%  1.00025
%  0.99635
%  1.02731
%  0.87965
%  1.33404
%  0.39944
%  1.69851
%  0.49263
%  1.20924
%  0.96259

[A, B] = eliminacionGaussiana(t2,b2,1)
% A =

%  -3.00000  -2.00000  -1.00000   0.01000
%   0.00000   0.99333   1.99667   3.00003
%   0.00000   0.00000  -0.02681  -0.04698
%   0.00000   0.00000   0.00000  -0.02503

% B =

%  -5.990000
%   5.990033
%  -0.073792
%  -0.025031

sustAtras(A,B)
%ans =

%  1.00000
%  1.00000
%  1.00000
%  1.00000

% Factorizacion LU

LU = LUdecomp(t0)
U = triu(LU)
L = tril(LU, -1) + eye(3)
y = sustAdel(L, b0)
x = sustAtras(U, y)

% x =

%  1.00000
%  1.00000
%  1.00000

LU = LUdecomp(t1)
U = triu(LU)
L = tril(LU, -1) + eye(13)
y = sustAdel(L, b1)
x = sustAtras(U, y)

% x =
%
%   1.00000
%   1.00000
%   0.99999
%   1.00054
%   0.99213
%   1.05928
%   0.73730
%   1.73252
%  -0.32227
%   2.54330
%  -0.12443
%   1.46497
%   0.91668


LU = LUdecomp(t2)
U = triu(LU)
L = tril(LU, -1) + eye(4)
y = sustAdel(L, b2)
x = sustAtras(U, y)
% x =

%  1.00000
%  1.00000
%  1.00000
%  1.00000

% Pregunta 2:

T0 = vander([0.5, 0.6, 0.7, 0.8, 0.9])
T1 = vander([0.5, 0.6, 7, 8, 9, 10, 13])
b0 = T0*[1;1;1;1;1]
b1 = T1*[1;1;1;1;1;1;1]

% Eliminacion Gaussana sin pivoteo

[A, B] = eliminacionGaussiana(T0,b0,0)
sustAtras(A,B)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

[A, B] = eliminacionGaussiana(T1,b1,0)
sustAtras(A,B)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

% Eliminacion Gaussana con pivoteo

[A, B] = eliminacionGaussiana(T0,b0,1)
sustAtras(A,B)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

[A, B] = eliminacionGaussiana(T1,b1,1)
sustAtras(A,B)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

% Factorizacion LU

LU = LUdecomp(T0)
U = triu(LU)
L = tril(LU, -1) + eye(5)
y = sustAdel(L, b0)
x = sustAtras(U, y)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

LU = LUdecomp(T1)
U = triu(LU)
L = tril(LU, -1) + eye(5)
y = sustAdel(L, b1)
x = sustAtras(U, y)

% ans =

%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

% Pregunta 3:

% Calculo de inversa por eliminacion Gauss-Jordan hasta reduccion a la matriz inversa
function Ai= inversa(A)
  [m,n] = size(A);
  Ai = eye(n);
  
  % Eliminacion de parte triangular inferior
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
  
  % Eliminacion de parte triangular superior
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

  % Reduccion de la diagonal a 1's
  for i =1:1:n
    p = A(i,i);
      for j = 1: 1: n
        Ai(i,j) = Ai(i,j)/p;
      endfor
      A(i,i) = A(i,i)/p;
  endfor
endfunction

% Pruebas

t0 = inversa(t0)
t11 = inversa(t1)
t22 = inversa(t2)
T00 = inversa(T0)
T11 = inversa(T1)

% norm(t0*t00-eye(3), Inf) =    1.9984e-15

% norm(t1*t11-eye(13), Inf) =  49.846

% norm(t2*t22-eye(4), Inf) =    9.1251e-12

% norm(T0*T00-eye(5), Inf) =    6.4944e-12

% norm(T1*T11-eye(7), Inf) =  0.0000027990