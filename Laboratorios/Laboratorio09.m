% Laboratorio 09
% Ejercicio 1:

% Definicion de la funcion f
y =@(x)  (x.*sin(x))./((x.*x) .+ 1);
dy =@(x) (((x.*(1 .+ x.*x)).*cos(x)) - ((x.*x-1).*sin(x)))./((1 .+ x.*x).*(1 .+ x.*x));

% Parte (a)
% Definicion de los intervalos; 
i0 = linspace(-pi*4,4*pi, 15);
i1 = linspace(-pi*4,4*pi, 25);
i2 = linspace(-pi*4,4*pi, 30);
i3 = linspace(-pi*4,4*pi, 50);

% coefiicientes de los polinomios de Newton con 30 y 50 ptos respectivamente:
new2 = newton(i2',y(i2)')
new3 = newton(i3',y(i3)')

% Parte (b):
% coefiicientes de los polinomios de Hermite con 15 y 25 ptos respectivamente:
[h0, z0] = hermite(i0', y(i0'), dy(i0'))
[h1, z1] = hermite(i1', y(i1'), dy(i1'))

% Parte (c):
% Definicion de los puntos de chebyshev (funcion incluida en la entrega):
c0 = chebyshev(15,-4*pi, 4*pi)
c1 = chebyshev(25,-4*pi, 4*pi)
c2 = chebyshev(30,-4*pi, 4*pi)
c3 = chebyshev(50,-4*pi, 4*pi)

% coefiicientes de los polinomios de Newton con 30 y 50 ptos respectivamente:
newc2 = newton(c2,y(c2)')
newc3 = newton(c3,y(c3)')

% coefiicientes de los polinomios de Hermite con 15 y 25 ptos respectivamente:
[hc0, zc0] = hermite(c0, y(c0'), dy(c0'))
[hc1, zc1] = hermite(c1, y(c1'), dy(c1'))

% Parte (d):
% En esta parte haremos las graficas.

% Particion para graficar la funcion f:
ii = linspace(-4*pi,4*pi, 500);

% En el enunciado se pide graficar todo sobre un mismo lienzo asi que aqui las graficas
% se ponen todas juntas (por mas que no sea lo mejor para la visualizacion). 
plot(ii, hornerN(new2, ii', i2'), "-")
hold on
plot(ii, hornerN(new3, ii', i3'), "-")
hold on
plot(ii, hornerN(h0, ii', z0), "-")
hold on
plot(ii, hornerN(h1, ii', z1), "-")
hold on 
plot(ii, hornerN(newc2, ii', c2), "-")
hold on
plot(ii, hornerN(newc3, ii', c3), "-")
hold on
plot(ii, hornerN(hc0, ii', zc0), "-")
hold on
plot(ii, hornerN(hc1, ii', zc1), "-")
hold on 
plot(ii,y(ii), "-")
hold on
plot(i0, y(i0),"*")
hold on
plot(i1, y(i1),"*")
hold on
plot(i2, y(i2),"*")
hold on
plot(i3, y(i3),"*")
hold on
plot(c0, y(c0),"*")
hold on
plot(c1, y(c1),"*")
hold on
plot(c2, y(c2),"*")
hold on
plot(c3, y(c3),"*")
hold on
title("Interpolaciones a f(x) mediante Newton y Hermite con 15, 25, 30, 50 ptos equidistantes y de Chebyshev")
legend("Newton: 30 ptos equidistantes", "Newton: 50 ptos equidistantes", "Hermite: 15 ptos equidistantes",
"Hermite: 25 ptos equidistantes","Newton: 30 ptos Chebyshev", "Newton: 50 ptos Chebyshev", "Hermite: 15 ptos Chebyshev",
"Hermite: 25 ptos Chebyshev", "Grafica de f(x)", "15 Puntos de prueba equidistantes", "25 Puntos de prueba equidistantes",
"30 Puntos de prueba equidistantes", "50 Puntos de prueba equidistantes", "15 Puntos de prueba de Chebyshev",
"25 Puntos de prueba de Chebyshev", "30 Puntos de prueba de Chebyshev", "50 Puntos de prueba de Chebyshev")
xlim([-pi*4,4*pi])
ylim([-0.4, 0.6])

% Parte (e):

% Para realizar interpolaciones por lo general no siempre es buena idea agregar mas
% nodos arbitrariamente. En este caso vimos que las interpolaciones con menor cantidad
% de puntos se aproximaban a la funcion a lo largo del intervalo [-pi*4,4*pi] mientras
% que si se agregaban mas puntos ellas se comportaban casi que iguales en un intervalo pequeno
% pero a medida que se acercaban a los bordes del intervalo [-pi*4,4*pi] la interpolacion
% empezaban a oscilar mucho lo cual hacia la aproximacion menos confiable, esto sucedia para
% todos los metodos y para puntos equidistantes. Cuando utilizabamos nodos de Chevyshev vimos
% que estas oscilaciones se reducian pero ellas se mantenian tambien se mantenian concentradas
% en los bordes del intervalo. Al considerar la derivada en la interpolacion vimos una
% mejora cuando teniamos pocos puntos en cuanto a la presicion, sin embargo hay
% que recordar que al incluir informacion de las derivadas se duplica el grado
% del polinomio resultante, lo cual nos trae el mismo problema que antes. Por lo tanto 
% si queremos incluir informacion sobre las derivadas lo mas recomendable es utilizar pocos
% nodos.

% Ejercicio 2:

tiempo = [0;3;5;8;13]
dist = [0;225;383;623;993]
vel = [75;77;80;74;72]

% Parte (a)
[ha, za] = hermite(tiempo, dist, vel);
t10 = hornerN(ha, 10, za)

% Parte (b)
it = linspace(0,13, 100);
% Para calcular la derivada en base de Newton utilizamos un programa que se consiguio
% en esta direccion: 
% https://math.stackexchange.com/questions/808738/algorithm-to-compute-newton-polynomial-derivative
% Que aqui le pusimos el nombre de "derivada".
dp = derivada(ha, it', za);
for i = 1:500
  if floor(dp(i)) == 80
    minimo = it(i)
    break
  endif
endfor

% El ciclo anteriar nos permite buscar el tiempo buscado que nos da que es 5.
% Entonces el tiempo minimo para alcanzar la velocidad de 80 pies/seg esta alrededor de 5.

% Parte(c):
% Para esta parte sencillamente tomaremos el maximo valor posible de dp en el intervalo
% con los puntos dados (desde 0 hasta 13).
max(dp)

% ans =  119.35
% Esta sera la maxima velocidad que posemos predecir en el intervalo de datos que tenemos 









