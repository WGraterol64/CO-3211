% Laboratorio 07:

% Modificacion para mostrar 15 decimales de precision
format long e 

% Ejercicio 1:

load("data.mat")

% Parte (a): 

p5 = polyfit(xd,yd,5)
p15 = polyfit(xd,yd,15)
p20 = polyfit(xd,yd,20)

% Calculo de errores cuadraticos medios:

% Para el polinomio de grado 5: 
Ec5 = sum((polyval(p5,xd)-yd).^2)/length(xd)


% Para el polinomio de grado 10: 
Ec15 = sum((polyval(p15,xd)-yd).^2)/length(xd)


% Para el polinomio de grado 15: 
Ec20 = sum((polyval(p20,xd)-yd).^2)/length(xd)

% Resultados
% Ec5 =    8.079756116110327e-02
% Ec15 =    4.435952150670028e-03
% Ec20 =    2.061003331691524e-03

% Parte (b): 

% Esto lo utilizamos para hacer que las curvas sean mas suaves:
x = linspace(min(xd), max(xd), 1000)

% Grafica de valores en data:
plot(xd,yd, "*b")
hold on

% Grafica de p5
plot(x, polyval(p5, x, "-r"))
hold on

% Grafica de p15
plot(x, polyval(p15, x, "-g"))
hold on

% Grafica de p20
plot(x, polyval(p20, x, "-y"))
hold on

% Etiquetas y ajustes:
title("Grafico de ajuste de curvas para las temperaturas en el tiempo")
legend("Valores originales", "Ajuste con polinomio de grado 5", "Ajuste con polinomio de grado 15", "Ajuste con polinomio de grado 20")
xlabel("Tiempo")
ylabel("Temperatura")
xlim([0,8])
ylim([33, 45])

% Parte (c):

% Primero vale notar que ninguno de estos polinomios son confiables pues octave 
% nos dice que el numero de condicion generado por la matriz para calcular los 
% valores de los polinomios es muy grande, sin embargo trabajando con estas aprox.
% Podemos notar que el polinomio con el menor error cuadratico es el de grado 20,
% sin embargo cuando vemos las graficas vemos que el de grado 20 en el intervalo
% [4,6] (donde no hay datos en xd) oscila con unas ondas muy grandes y alejados 
% de los datos de entrada. En cambio el de grado 15 posee mas error cuadratco pero
% ella es una curva mas suave y a lo largo del intervalo [0,8] se comporta mas 
% cercano a los valores de xd-yd que el de grado 20 a lo largo de todo el intervalo 
% [0,8]. El de grado 5 es una peor aproximacion en general pues tiene mayor error 
% cuadratico y se aleja mas de los puntos esperados que los oros dos polinomios.
% Por lo tanto para hacer aproximarnos a los valores exactos de las entradas podemos
% utilizar el polinomio de grado 20 pues esta muy cercano de los valores, pero
% para hacer predicciones (que es lo que nos importa en este laboratorio)
% con ese modelo nos conviene mas utilizar el polinomio de grado 15 pues es el 
% que menos oscilaciones tiene y no se aleja mucho de los valores de entrada.

temp1 = polyval(p15,4.5)
temp2 = polyval(p15,5)
temp3 = polyval(p15,5.5)

% Resultado
% temp1 =    3.687829458276998e+01
% temp2 =    3.676211299827744e+01
% temp3 =    3.669325058325114e+01

% Parte (d)
p15(end) = p15(end)-36.612
tiempos = roots(p15)

% tiempos =
%
%   1.001332606706059e+01 + 0.000000000000000e+00i
%   8.228123448954827e+00 + 5.337061500553587e-01i
%   8.228123448954827e+00 - 5.337061500553587e-01i
%   7.141819442978571e+00 + 8.602717817795251e-01i
%   7.141819442978571e+00 - 8.602717817795251e-01i
%   5.349624757560823e+00 + 1.736151413182860e-01i
%   5.349624757560823e+00 - 1.736151413182860e-01i
%   4.099116961905411e+00 + 0.000000000000000e+00i
%   3.296382033527212e+00 + 1.307877990527400e+00i
%   3.296382033527212e+00 - 1.307877990527400e+00i
%   1.071746686220147e+00 + 7.230452508162215e-01i
%   1.071746686220147e+00 - 7.230452508162215e-01i
%   2.069208913129984e-01 + 4.615115540128579e-01i
%   2.069208913129984e-01 - 4.615115540128579e-01i
%  -1.556139108347021e-01 + 0.000000000000000e+00i

% Entonces obtenendremos una temperadura de 36.612 cuando hayan transcurrido 
% 4.099116961905411e+00 minutos (en el intervalo entre
% [0, 8] minutos.