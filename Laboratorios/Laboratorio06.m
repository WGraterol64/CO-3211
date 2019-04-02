% Laboratorio 06:

% Ejercicio 1:
A = [[-1+j, 0, 1/4];[1/4, 1, 1/4];[1, 1, 3]]
B = [[1, 3, -2];[-1, -2, 3];[1, 1, 2]]

% Funcion para dibujar circuclos:
% (Extraida de internet)
function h = circle(x,y,r)
   hold on
   th = 0:pi/50:2*pi;
   xunit = r * cos(th) + x;
   yunit = r * sin(th) + y;
   h = plot(xunit, yunit);
   hold off
end
eia = (A)
% ega =

%  -1.05404 + 0.98885i
%   3.17798 + 0.01412i
%   0.87606 - 0.00297i

eib = eig(B)
% egb =

%  -0.81446 + 1.27238i
%  -0.81446 - 1.27238i
%   2.62892 + 0.00000i

% Radio de los discos
ra1 = 0;
ra2 = 0;
ra3 = 0;

rb1 = 0;
rb2 = 0;
rb3 = 0;

for i = 1:3
  if i ~= 1
      ra1 = ra1 + abs(A(1,i));
      rb1 = rb1 + abs(B(1, i));
  endif
  if i ~= 2
      ra2 = ra2 + abs(A(2,i));
      rb2 = rb2 + abs(B(2, i));
  endif
  if i ~= 3
      ra3 = ra3 + abs(A(3,i));
      rb3 = rb3 + abs(B(3, i));
  endif
endfor

% Graficar:

circle(-1, 1, ra1)
hold on
circle(1, 0, ra2)
hold on
circle(3, 0, ra3)
hold on
circle(1, 0, rb1)
hold on
circle(-2, 0, rb2)
hold on
circle(2, 0, rb3)
hold on
xlim([-8, 8])
ylim([-7, 7])

% Ejercicio 2: 
%
% Los programas pedidos se anexan en el correo bajo los nombres de "met_dir",
% "met_inv" y "met_dir_desp". En el programa del metodo desplazado
% se agrega ademas la entrada de el desplazamiento (No tiene sentido no ponerla
% y tener que modificarla para cada desplazamiento)

% Ejercicio 3:

% Autovalores maximo y minimo de A
[eaM, uaM, itaM] = met_dir(A, [1;1;1]) % eigval = 3.177978 + 0.014124i; Numero max de iteraciones
[eam, uam, itam] = met_inv(A, [1;1;1]) % eigval = 0.8760614 - 0.0029702i; Numero max de iteraciones

eig(A)
% ans =

%  -1.05404 + 0.98885i
%   3.17798 + 0.01412i
%   0.87606 - 0.00297i

% Autovalores maximo y minimo de B
[ebM, ubM, itbM] = met_dir(B, [1;1;1]) % eigval = 2.6289; Numero de iteraciones = 23
[ebm, ubm, itbm] = met_inv(B, [1;1;1]) % eigval = -4.4286; Numero max de iteraciones

eig(B)
% ans =

%  -0.81446 + 1.27238i
%  -0.81446 - 1.27238i
%   2.62892 + 0.00000i

% Se tiene que el unico caso para el que no nos aproximamos a la solucion es en
% la matriz B y esto de puede deber a que ella posee dos autovalores que en magnitud 
% son iguales y que son candidatos a ser los menores. Estos errores se deben a que 
% todos estos metodos iterativos asumen que hay un autovalor mayor o menor que todos. 
% Entonces al no cumplirse esto no podemos asegurar la convergencia de los metodos pues
% el residuo de aplicar las operaciones no tiende a cero. 
% Por lo tanto tampoco se contracide lo que se dio en teoria pues necesitamos que
% se cumplan las hipotesis para ello.

% Ejercicio 4:

% Generacion de las matrices:

A1 = eye(3)
A2 = eye(3)
A8 = eye(3)

A1(1,2) = 1
A2(1,2) = 1
A8(1,2) = 1
A1(3,2) = 1
A1(2,3) = 1
A2(3,2) = 10e-2
A2(2,3) = 10e-2
A8(3,2) = 10e-8 
A8(2,3) = 10e-8

% Vector inicial: 
% Se elige un vector no nulo para poder dividir entre sus coordenadas
x = [1;1;1]

% Calculo de autovalores
[e1, u1, it1] = met_dir(A1, x) % ans: eigval = 2; iteraciones = 1
[e2, u2, it2] = met_dir(A2, x) % ans: eigval = 1.1001; iteraciones = 72
[e8, u8, it8] = met_dir(A8, x) % ans: eigval = 1.0032; iteraciones = 316


% Calculo de autovalores reales con eig para calcular el error relativo

e11 = max(abs(eig(A1))) % = 2
e22 = max(abs(eig(A2))) % = 1.1000
e88 = max(abs(eig(A8))) % = 1.0000

% Error relativo:

error1 = abs(e1-e11)/abs(e11) % = 0
error2 = abs(e2-e22)/abs(e22) % = 0.000094280
error8 = abs(e8-e88)/abs(e88) % = 0.0031645

% Notese que aqui la estabilidad de los resultados dependen de la estabilidad de 
% el sistema Ax=y pues lo que estaremos haciendo sera multiplicar esto por A en
% cada iteracion. Entonces al tener numeros mas pequenos en la matriz cuando hagamos
% estos productos iremos acomulando error y este error dependera de las operaciones
% aritmeticas que hagamos. En particular en el caso de la matriz A8 el multiplicar
% por 10e-8 muchas veces ira generando mucho mas error en cada iteracion que multiplicar
% por 1 y estos errores numericos generaran que el error relativo a la salida aumente. 

% Ejercicio 5:

A = [[-149, -50, 154]; [537, 180, 546]; [-27, -9, -25]]

eigMin = met_inv(A, x)
% Resultado: Retorna la respuesta esperada pues da 1.00001 y el menor 
% autovalor es 1.00000 +  0.00000i

% Ahora para calcular el otro autovalores utilizaremos el metodo desplazado
% para ello haremos los circulos de Gershgorin para saber cual valor de delta tomar

% Calculo de radios
r1 = 0;
r2 = 0;
r3 = 0;

for i = 1:3
  if i ~= 1
      r1 = r1 + abs(A(1,i));
  endif
  if i ~= 2
      r2 = r2 + abs(A(2,i));
  endif
  if i ~= 3
      r3 = r3 + abs(A(3,i));
  endif
endfor

% Graficar:

circle(-149, 0, r1)
hold on
circle(537, 0, r2)
hold on
circle(-27, 0, r3)
hold on
xlim([-1000, 1000])
ylim([-1000, 1000])


% Aqui podemos ver que tendremos los 3 circulos, uno adentro del otro, entonces
% sospechamos que los autovalores estan adentro de ellos y debemos desplazar la matriz.
% Entonces desplazemos la matriz para arriba y calculemos:

met_dir_desp(A,x,100+20j)

% numero maximo de iteraciones alcanzadas
% ans =   2.5000 - 91.1907i

% Esto anterior es una buena aproximacion al autovalor verdadero. Como tenemos
% una raiz compleja sabemos que su conjugado tambien es una raiz, si queremos
% calcularlo de igual manera podemos desplazar la matriz para abajo ahora y luego: 

met_dir_desp(A,x,100-20j)
% numero maximo de iteraciones alcanzadas
% ans =   2.5000 + 91.1907i

% Entonces ya hemos hayado los 3 autovalores.
% Notese que no podiamos utilizar el metodo de la potencia directo pues al igual que
% En la pregunta 3 teniamos dos autovalores con la misma magnitud y entonces fallaban
% los metodos.

% Si queremos ver la prueba corremos
% Si queremos ver la prueba corremos
eigMax = met_dir(A, x) 
% Resultado: llega al numero maximo de iteraciones y no se aproxima a la respuesta
% esperada pues retorna 64.490 y el resultado es  2.50000 + 91.19073i. 


% Para referencia: los autovalores verdaderos son
eig(A)
% ans =
%
%    2.50000 + 91.19073i
%    2.50000 - 91.19073i
%    1.00000 +  0.00000i