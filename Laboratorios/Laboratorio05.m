% Laboratorio 05:

% Pregunta 1: 
% Se anexan en el correo las rutinas de jacobi y gauss_seidel
% Esta es la ejecucion de ellos sobre un caso de prueba
% Ambas rutinas fueron disenadas para que retornen no solo el vector 
% solucion, si no tambien el tiempo de ejecucion y numero de iteraciones.
A = hilb(11)+1.025*eye(11);
b = A*ones(11,1);
maxiter = 1500;
tol = 1e-11;
xv = zeros(11,1);

% Ejecucion de Jacobi:
jacobi(A, b, xv, maxiter,tol)

% Output: 
% Solucion Hallada! Numero de iteraciones:
% 1402
% ans =
%
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
%   1.0000
   
% Ejecucion de gauss_seidel
gauss_seidel(A, b, xv, maxiter,tol)

% Output: 
% Solucion Hallada! Numero de iteraciones:
% 17
% ans =
%
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000
%   1.00000

% Pregunta 2:

% En el correo se anexan 3 programas para la resolucion por LU. Uno que descompone
% la matriz en una L y una U (ambas guardadas en la misma matriz). Uno que recibe
% dos matrices y entonces calcula la resolucion por sustitucion hacia adelante y
% sustitucion hacia atras y retorna el vector solucion. Y una rutina principal que
% llama a esas dos y luego retorna el vector solucion y el tiempo de ejecucion.


% Pregunta 3: 

% Cargamos los datos
load("datos20")

% Separamos B en 3 vectores
b1 = NaN*ones(480,1);
b2 = NaN*ones(480,1);
b3 = NaN*ones(480,1);

for i = 1:480
  b1(i,1) = B(i,1);
  b2(i,1) = B(i,2);
  b3(i,1) = B(i,3);
endfor

% Recordar se modificaron las rutinas para que no solo 
% retornen el vector solucion si no tambien el tiempo y numero
% de iteraciones. Esto nos ayudara a responder la pregunta 4

% Resolucion de los 3 sistemas por LU:
[sol1, time1 ] =LUcomplete(A,b1)
[sol2, time2 ] =LUcomplete(A,b2)
[sol3, time3 ] =LUcomplete(A,b3)

% Resolucion de los 3 sistemas por jacobi:
% No se especifica el numero maximo de iteraciones asi
% que se elige 1e5

% Matrices para ejecutar Jacobi

% Creacion de la matriz de Jacobi

D = diag(diag(A));
E = -tril(A,-1);
F = -triu(A,1);

M = D;
N = F+E;
H = inv(M)*N;

x0 = ones(480, 1);
x1 = zeros(480, 1);
tol = 1e-5;

% Pruebas
[sol4, k4, time4 ] = jacobi(A, b1, x0, 1e5,tol)
[sol5, k5, time5 ] = jacobi(A, b2, x0, 1e5,tol)
[sol6, k6, time6 ] = jacobi(A, b3, x0, 1e5,tol)
[sol7, k7, time7 ] = jacobi(A, b1, x1, 1e5,tol)
[sol8, k8, time8 ] = jacobi(A, b2, x1, 1e5,tol)
[sol9, k9, time9 ] = jacobi(A, b3, x1, 1e5,tol)

% Resolucion de los 3 sistemas por gauss_seidel:
% No se especifica el numero maximo de iteraciones asi
% que se elige 1e5

% Creacion de la matriz GS
M = (D-E);
N = F;
GS = inv(M)*N;

% Pruebas
[sol10, k10, time10 ] = gauss_seidel(A, b1, x0, 1e5,tol)
[sol11, k11, time11 ] = gauss_seidel(A, b2, x0, 1e5,tol)
[sol12, k12, time12 ] = gauss_seidel(A, b3, x0, 1e5,tol)
[sol13, k13, time13 ] = gauss_seidel(A, b1, x1, 1e5,tol)
[sol14, k14, time14 ] = gauss_seidel(A, b2, x1, 1e5,tol)
[sol15, k15, time15 ] = gauss_seidel(A, b3, x1, 1e5,tol)

% Pregunta 4:

% En teooria se discutio un teorema que dice "Para cualquier vector inicial
% el metodo iterativo converje si y solo si el radio espectral es menor que 1"
%
% Entonces calculemos el redio espectral de H y GS:
 
max(abs(eig(H)))
% ans =  1.0044

max(abs(eig(GS)))
% ans =  0.91925

% Entonces en base al teorema anterior tendremos que con el metodo de Jacobi nunca
% convergera el sistema de ecuaciones planteado, por lo menos para cualquier
% vector inicial. El metodo convergeria si elegimos a priori un vector muy cercano
% a la solucion para el iterado inicial.

% En el caso de Gauss-Seidel se espera que si bien la sucesion convergera, lo hara
% de forma muy lenta pues su radio espectral esta muy cercano a 1. A medida que el
% el radio espectral se aleja de 0 se requieren de mas iteraciones para alcanzar el 
% resultado para una entrada arbitraria. Entonces si bien terminaremos en un numero
% finito de iteraciones, el numero de iteraciones requeridas tendera a un numero grande.
% por ejemplo para el caso de prueba de la linea 126 tendremos que el numero de 
% iteraciones es k10 =  117 y que el tiempo de ejecucion sera time10 =  1164.2
% lo cual es mejor cuando lo comparamos con el de tiempo de ejecucion de LU.

% Entonces el mejor metodo para calcular estos sistemas de ecuaciones depende de
% cual sea la matriz. Si sabemos que el radio espectral de la matriz es pequeno
% entonces gauss_seidel es mucho mejor que LU pero en general es mejor LU
%pues si bien es costoso, nos asegura llegar a la solucion en un numero exacto de 
% pasos.











