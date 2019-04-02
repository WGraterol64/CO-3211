% Pregunta1:
a = NaN
b = NaN
x = NaN
A = NaN
B = NaN
n= NaN
% Parte a: Algoritmo de elminacion Gaussiana con y sin pivoteo y sust hacia atras

function x = matrizTriangular(a,b)
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
    
function [a, b]= eliminacionGaussiana(a,b, p)
    % Sin Pivoteo
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
    % Con pivoteo
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
 

% Parte b:

A1=[[10e-8 , -1, 10e10]; [10, 10e-8 - 10e7, 1]; [1, -1, 1]]
b1 = A1*[1;-1;1]
A2=[[4 , -1, 1]; [10e-9, 10e-6, 10e-9]; [10e-16, -10e-16, 1]]
b2 = A2*[1;-1;1]

% El determinante de A1 = 10e17 - 10e11 - 10e2 -10 -1 - 10e-1 +10e-8 + 10e-16
% y el determinante de A2 = 4*10e-6 +10e-9 - 10e-22 + 2*10e-25 y ambos son distintos
% de 0 por lo que los sistemas tienen solucion.

% Parte c:

%  Sistema A1x=B1:
%  Sin pivoteo:
[a1, B1] = eliminacionGaussiana(A1, b1, 0)
x10 = matrizTriangular(a1, B1)
%  x10 =
%  10000000
%         0
%         1
%
% Error:
norm(([1;-1;1]-x10)/norm([1;-1;1], Inf),Inf)
%  ans = 9999999
%         
%  Con pivoteo:
[a1, B1] = eliminacionGaussiana(A1, b1, 1)
x11 = matrizTriangular(a1, B1)
%  x11 =
%
%  1
% -1
%  1
%
%  Error: 
norm(([1;-1;1]-x11)/norm([1;-1;1], Inf),Inf)
%  ans = 0  
%

%  Sistema A2x=B2 :
%  Sin pivoteo:
[a2, B2] = eliminacionGaussiana(A2, b2, 0)
x20 = matrizTriangular(a2, B2)
%  x20 =
%
%  1.00000
% -1.00000
%  1.00000
%  
% Error:
norm(([1;-1;1]-x20)/norm([1;-1;1], Inf),Inf)
%  ans =    1.1102e-16
%
%  Con pivoteo:
[a2, B2] = eliminacionGaussiana(A2, b2, 1)
x21 = matrizTriangular(a2, B2)
%  x21 =
%
%  1.00000
% -1.00000
%  1.00000
%  
% Error:
norm(([1;-1;1]-x21)/norm([1;-1;1], Inf),Inf)
%  ans =    1.1102e-16
%  
max(abs(eig(A1))) 
%  ans = 1.0000e+08
%
min(abs(eig(A1))) 
%  ans = 3.1623e+05
%
max(abs(eig(A2))) 
%  ans= 4.0000
%
min(abs(eig(A2))) 
% ans = 1.0003e-05
%
%  En el segundo sistema de ecuaciones vemos un error pequeno por cualquiera de 
%  dos metodos. Sin embargo en el primer sistema vemos una diferencia gigante
%  esta diferencia es causada porque si hacemos la eliminacion tal como esta la 
%  matriz vamos a empezar dividiendo por cifras muy pequenas y poco precisas lo 
%  que hara que ese error aumente con cada operacion de manera rapida, en cambio 
%  al momento de usar el pivoteo se reduce bastante error pues se evita manipular
%  estos numeros lo menos posible. Si el computador trabajara con aritmetica exacta
%  y con un numero infinito de cifras significativas esto no seria necesario pues
%  no seria necesario redondear o truncar numeros y por lo tanto, no se generaria
%  error pues siempre mantendriamos la definicion exacta de cada numero. En estos
%  computadores inexactos influye el pivoteo pues nos ayuda a disminuir el error
%  y obtener resultados mas fieles a lo esperado
%
% Los autovalores en la primera matriz al ser tan elevados hacen que en la matriz
% inversa sean muy pequenos y esto lo que quiere decir es que la matriz inversa 
% va a tender a ser singular y por lo tanto el numero de condicion sea muy grande.
% Entonces podemos asociar autovalores extremos con mal condicionamiento de la matriz.

% Parte d:

% condA1 
condA1 = norm(A1, Inf)*norm(inv(A1), Inf)
% ans = 1.0000e+11
%
% condA2:
condA2 = norm(A2, Inf)*norm(inv(A2), Inf)
% ans = 599850.04349
%
 norm((b1-A1*x10), Inf)/norm(b1, Inf)
% ans =  0.000100000
 norm((b1-A1*x11), Inf)/norm(b1, Inf) 
% ans = 0
 norm((b1-A2*x20), Inf)/norm(b1, Inf) 
% ans =  1.00000
norm((b1-A2*x21), Inf)/norm(b1, Inf)
% ans =  1.00000

% Entonces nos queda que:

condA1*norm((b1-A1*x10), Inf)/norm(b1, Inf) 
% ans = 10000000
%
condA1*norm((b1-A1*x11), Inf)/norm(b1, Inf)
% ans = 0
%
condA2*norm((b1-A2*x20), Inf)/norm(b1, Inf) 
% ans = 599850.04349
%
condA2*norm((b1-A2*x21), Inf)/norm(b1, Inf) 
% ans = 599850.04349


%Y luego efectivamente para todos los casos se tiene que: 
%Error = 9999999 <= 10000000
%Error = 0 <= 0
%Error =    1.1102e-16 <= 599850.04349
%Error =    1.1102e-16 <= 599850.04349

% Ejercicio 2:

% Parte a:

% Funcion que genera el vector de Hilbert
function A = Hilbert(n)
  A = NaN*ones(n,n);
  for i = 1 : 1 : n
    for j = 1 : 1 : n 
      A(i,j) = 1/(i+j-1);
    endfor
  endfor
endfunction

% Funcion que genera el vector asociado a la matriz de Hilbert
function b = HilbertVector(n)
  A = Hilbert(n);
  b = A*ones(n,1);
endfunction

% Parte a y b:
function HilbertAndError(n)
  for i = 2:1:n
    % Parte a: Genera la matriz con su vector asociado y resuelve el sistema
    % de las dos maneras pedidas.
    A = Hilbert(n)
    b = HilbertVector(n)
    [A0,b0] = eliminacionGaussiana(A,b,1)
    x0 = matrizTriangular(A0,b0)
    x = A\b
    
    % Parte b: Muestra la norma del error en pantalla para cada n
    r = b-A*x0
    dx = x0-x
    norm(r, Inf)
    norm(dx, Inf)
    % por ejemplo para n == 3 se obtienen los siguientes resultados:
    % r = [2.2204e-16; 0.0000e+00; 0.0000e+00]
    % dx = [-6.6613e-16; 2.6645e-15; -2.5535e-15]
    % norm(r,Inf) =  2.2204e-16
    % norm(dx, Inf) = 2.6645e-15
    
    % para n == 5
    % r = [4.4409e-16; 0.0000e+00; 2.2204e-16; 1.1102e-16; 1.1102e-16]
    % dx = [4.3965e-14; -8.7419e-13; 3.8874e-12; -5.9914e-12; 2.9738e-12]
    % norm(r,Inf) =  4.4409e-16
    % norm(dx, Inf) =  5.9914e-12
    
    % para n == 7 
    % r = [4.4409e-16; 0.0000e+00; 2.2204e-16; 0.0000e+00; -2.2204e-16; -1.1102e-16
    % ;1.1102e-16; 0.0000e+00]
    % dx = [4.2433e-12; -2.0031e-10;  2.3706e-09; -1.1859e-08;  2.9935e-08
    % ; -4.0139e-08; 2.7291e-08; -7.4036e-09]
    % norm(r,Inf) =  4.4409e-16
    % norm(dx, Inf) = 0.000000040139

  endfor
endfunction

% Parte c:

% Para una computadora de 32 bits, basta con que el condicionamiento este en el
% orden de 10e8 y esto ocurre a partir de n = 6 (el condicionamiento para la 
% matriz de Hilbert de tamano 6 es 29070279.00408)

% Parte d:

% Como no se conoce el n se da una formula en base un n dao que luego el lector 
% puede probar

% sin embargo aqui hay algunos casos de prueba:
% para n = 3 => cond(Hilbert(3)) = 524.06
% para n = 5 => cond(Hilbert(5)) = 476607.25024
% para n = 7 => cond(Hilbert(7)) = 475367356.36862

for n = 3:2:7
  r = 2: 1 :n
  y = NaN*ones(1,n-1)
  for i = 1:n-1
    y(i) = cond(Hilbert(i+1))
  endfor
  figure(n)
  plot(r, y)
  xlim([2, n])
endfor


% Parte e:

% Podemos ver que a medida que n aumenta el condicionamiento de la matriz aumenta 
% y esto se debe a que cada vez que incrementamos el tamano de n se van agregando
% entradas mas pequenas a la matriz lo que hace que su manipulacion con la aritmetica
% del computador sea mas propensa a errores, el aumento de error es exponencial 
% lo que nos lleva a concluir que la matriz de Hilbert es una matriz muy mal 
% condicionada para n pequenos y es casi imposible de manejar para n grandes

