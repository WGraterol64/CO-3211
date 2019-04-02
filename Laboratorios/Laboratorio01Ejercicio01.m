% Laboratorio 1:
% CO-3211
% Wilfredo Graterol
% 15-10639
% Seccion:1

% Ejercicio 1: 
y=@(x) (1-cos(x))./(x.*x);

% Parte (a): 
x = -3:0.1:3;
plot(x,y(x),'-b')

% Parte (b):
y((1.2).*(10.^-8))

% Observamos que el resultado esta fuera del rango esperado pues 
% en el intervalo [-3,3] el maximo valor que tomaba la funcion era 0.5.
%
% El resultado no es confiable, incluso si no supieramos que el valor
% esta eroneo. Este error es generado por la manera en la que el computador
% maneja numeros muy grandes y lo que le hace a ello operar numeros muy
% pequenos. Al tener una resta que involucre numeros pequenos junto con 
% divisiones entre cosas muy pequenas nos genera errores de redondeo
% en cada operacion que es lo que genera el resultado erroneo.

% Parte (c):
y=@(x) (1/2).*(((sin(x/2))./(x/2)).^2);
y((1.2).*(10.^-8))

% Aqui nos da un resultado esperado, lo que sucede es que al tener menos 
% operaciones que involucren cifras muy pequenas se acomula menos error,
% entonces es esa diferencia en la acomulacion de error lo que hace que 
% esta representacion (el algoritmo para calcular la funcion) sea mejor 
% que la anterior.
