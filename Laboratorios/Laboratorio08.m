% Laboratorio 08
% Ejercicio 1:

% Definicion de la funcion f
f =@(x)  (x.*sin(x))./((x.*x) .+ 1);

% Definicion de los intervalos; 
i1 = linspace(-pi*4,4*pi, 10);
i2 = linspace(-pi*4,4*pi, 30);
i3 = linspace(-pi*4,4*pi, 50);

% Parte (a):

% Para esto solo guardaremos los coeficientes del polinomio
% como a(n),a(n-1),...,a(1). Como tenemos permitido utilizar
% la funcion de matlab para calcular la matriz de vandermonde
% la utilizaremos y con ella calcularemos los coeficientes.

a1 = linsolve(vander(i1'),f(i1)')
a2 = linsolve(vander(i2'),f(i2)')
a3 = linsolve(vander(i3'),f(i3)')

% Parte(b):

% Como se pidio, se provee el codigo del polinomio de Lagrange implementado
% con el codigo de neville aparte. Este polinomio se evaluara para las graficas.

% Parte (c):

% En esta parte haremos las 3 graficas.

% Particion para graficar la funcion f:
ii = linspace(-4*pi,4*pi, 500);

% Grafica para la aproximacion con 10 puntos:

% En la grafica 1 vemos que los pos polinomios coinciden perfectamente.
subplot(3,2,[1,2])
plot(ii, horner(a1, ii'), "-")
hold on
plot(ii, neville(ii, i1, f(i1)), "-")
hold on
plot(ii,f(ii), "-")
hold on
plot(i1, f(i1),"*")
hold on
title("Aproximaciones a f(x) con los polinomios de Vandermonde y Lagrange: 10 puntos de prueba")
legend("Polinomio de Vandemonde", "Polinomio de Lagrange", "Grafica de f(x)", "Puntos de prueba")
xlim([-pi*4,4*pi])

% Grafica para la aproximacion con 30 puntos:
subplot(3,2,[3,4])
plot(ii, horner(a2, ii'), "-")
hold on
plot(ii, neville(ii, i2, f(i2)), "-")
hold on
plot(ii,f(ii), "-")
hold on
plot(i2, f(i2),"*")
hold on
title("Aproximaciones a f(x) con los polinomios de Vandermonde y Lagrange: 30 puntos de prueba")
legend("Polinomio de Vandemonde", "Polinomio de Lagrange", "Grafica de f(x)", "Puntos de prueba")
xlim([-pi*4,4*pi])
ylim([-0.4, 0.6])
% Aqui se tuvieron que ajustar el eje y para que se pudiera ver la grafica bien escalada
% y se pudiera ver bien la diferencia entre las graficas de los polinomios y f.

% Grafica para la aproximacion con 50 puntos:
subplot(3,2,[5,6])
plot(ii, horner(a3, ii'), "-")
hold on
plot(ii, neville(ii, i3, f(i3)), "-")
hold on
plot(ii,f(ii), "-")
hold on
plot(i3, f(i3),"*")
hold on
title("Aproximaciones a f(x) con los polinomios de Vandermonde y Lagrange: 50 puntos de prueba")
legend("Polinomio de Vandemonde", "Polinomio de Lagrange", "Grafica de f(x)", "Puntos de prueba")
xlim([-pi*4,4*pi])
ylim([-0.4, 0.6])
% Aqui se tuvieron que ajustar el eje y para que se pudiera ver la grafica bien escalada
% y se pudiera ver bien la diferencia entre las graficas de los polinomios y f.

% Parte (d):

% Evaluacion de la interpolacion con 50 puntos con el polinomio de Vandermonde
h1 = horner(a3, -6.1333)
% h1 = -0.0000000055684
h2 = horner(a3, -1.4142)
% h2 = -5.7081e-34

% Evaluacion de la interpolacion con 50 puntos con el polinomio de Lagrange
l1 = neville(-6.1333, i3, f(i3))
% l1 =  0.035793
l2 = neville(-1.4142, i3, f(i3))
% l2 =  0.46368

% Calculos de errores relativos:
% Para Vandermonde:

ev1= norm(f(-6.1333)-h1, Inf)/norm(f(-6.1333),Inf)
% ev1 = 1.00000
ev2= norm(f(-1.4142)-h2, Inf)/norm(f(-1.4142),Inf)
% ev2 =  1

% Para Lagrange:

el1= norm(f(-6.1333)-l1, Inf)/norm(f(-6.1333),Inf)
% el1 =  2.5092
el2= norm(f(-1.4142)-l2, Inf)/norm(f(-1.4142),Inf)
% el2 =  0.0042002

% Parte (e):

% Antes de empezar el analisis note que:
cond(vander(i1'))
% ans =  9831931952.60850
cond(vander(i2'))
% ans =    3.7574e+32
cond(vander(i3'))
% ans =    2.6639e+55

% Entonces, en las graficas pudimos notar que para los 3 casos de prueba
% (con 10, 30 y 50 puntos) se obtuvieron las mejores aproximaciones con la interpolacion
% utilizando 10 puntos. Las interpolaciones con 30 y 50 puntos se alejaban demasiado de
% la funcion como para ser una aproximacion aceptable, recordemos que al tener mas puntos
% tendremos un polinomio de mas grado y con ello tendremos mas raices y mnas oscilaciones. 
% Entre las dos interpolaciones vemos que la de Lagrange se aproxima mejor a la funcion 
% en general para las 3 aproximaciones. Tambien vimos que interpolacion de Vandermonde 
% despues de 10 puntos se alejaba de la funcion original incluso en los puntos con los que 
% se hacia la interpolacion. Esto se debe al mal condicionamiento de la matriz de Vandermonde,
% como vemos arriba, estas son matrices que son muy mal condicionadas. Entoncer los resultados
% numericos que obtenemos no son los que deberiamos obtener al aplicar el metodo. Recordemos
% que como estamos haciendo una interpolacion de grado n con n+1 puntos deberiamos obtener
% un polinomio unico, es decir, que los polinomios que obtenemos mediante Lagrange y Vandemonde
% deberian ser los mismos, pero aqui vemos esta diferencia justamente por el error numerico
% generado por el condicionamiento de Vandermonde. Estos dos metodos de aproximacion tienen
% la caracteristica de ser muy ineficientes, sin embargo podemos rescatar como ventaja que 
% son sencillos de implementar ambos. Aunque para el caso de Vandermonde tenemos el problema
% del condicionamiento, tambien tenemos el problema de que es costoso agregar otro punto
% para hacer la aproximacion. En cambio para el metodo de Lagrange tenemos que tiene una 
% implementacion mas dificil, da resultados mas confiables y es mas facil agregar otro punto
% para hacer la aproximacion.
