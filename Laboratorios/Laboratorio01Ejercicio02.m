% Ejercicio 2:

x = NaN;
n = NaN;
function tcos = taylorCos(x, n)
  tcos = 0;
  for m =0:n
    tcos = tcos .+ ((-1).^m).*((x.^(2*m))./factorial(2*m));
  endfor
endfunction

% Parte (a)
function y = recursiveCos(x,n)
  y=(-1)^n;
  for i = n:-1:1
    y = (y.*(x.^2)./((2*i)*(2*i-1))) .+ (-1)^(i-1);
  endfor
endfunction

% Parte (b)
x = 8*pi:0.1:14*pi;
n = 51;

plot(x, taylorCos(x,n),'-b')
xlim([8*pi, 14*pi])
ylim([-1, 1])
plot(x, recursiveCos(x,n),'-b')
xlim([8*pi, 14*pi])
ylim([-1, 1])

% Parte (c)

taylorCos(2*pi, n)
recursiveCos(2*pi, n)

% ans =  1.00000
% ans =  1
% Error relativo: 0 en ambos casos

taylorCos(8*pi, n)
recursiveCos(8*pi, n)

% ans =  1.0000
% ans =  1.00000
% Error relativo: 0 en ambos casos
taylorCos(50.265, n)
recursiveCos(50.265, n)

% ans =   -6.7520e+10
% ans =   -6.7520e+10 
% Error relativo: aproximadamente 1 en ambos casos


% Parte (d):

% Estos resultados probablemente se deben a que los algoritmos utilizados
% para realizar estas aproximaciones son inestables y por lo tanto, 
% si bien aproximan a la funcion para valores pequenos, para valores muy 
% grandes posee un error muy alto (recordemos que si bien pi es un numero 
% pequeno en terminos de su magnitud, ella posee mas decimales de los que 
% se pueden representar en un computador y esto da a pie a errores). 
% Podemos decir que el segundo algoritmo es una mejor aproximacion pues 
% se mantiene estable para mas valores como podemos ver en la primera 
% prueba, pero ninguna de las dos aproximaciones es optima. Para mejorar
% estos resultados es necesario cambiar el algoritmo.

% Parte (e)
function tcos = taylorCosM(x, n)
  tcos = 0;
  xx = x - floor(x./(2*pi))*2*pi;
  for m =0:n
    tcos = tcos .+ ((-1).^m).*((xx.^(2*m))./factorial(2*m));
  endfor
endfunction

function y = recursiveCosM(x,n)
  y=(-1)^n;
  xx = x - floor(x./(2*pi))*2*pi;
  for i = n:-1:1
    y = (y.*(xx.^2)./((2*i)*(2*i-1))) .+ (-1)^(i-1);
  endfor
endfunction

taylorCosM(2*pi, n)
recursiveCosM(2*pi, n)

% ans =  1
% ans =  1
% Error relativo: 0 en ambos casos

taylorCosM(8*pi, n)
recursiveCosM(8*pi, n)

% ans =  1
% ans =  1
% Error relativo: 0 en ambos casos

taylorCosM(50.265, n)
recursiveCos(50.265, n)

% ans =  1.00000
% ans =   -6.7520e+10 
% Error relativo: 0 en el primer caso y aproximadamente 1 en el segundo

% En este caso el error fue alterado en la ultima prueba en el programa de
% la serie de taylor, este error se debe a la difrencia producida por el 
% nuevo valor de x que se esta evaluando. Esta prueba mejora los resultados
% para la primera funcion por la manera en la que se esta tomando el valor 
% de x lo que produce un error inicial dado pero que se crezca mucho mas 
% lento que el de la funcion polinomial original.
