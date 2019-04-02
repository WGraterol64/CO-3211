% Ejercicio 3:

x = NaN;
function xprom = promedio(x);
  xprom = 0;
  for i = 1:length(x)
    xprom = xprom .+ x(i);
  endfor
  xprom = xprom ./ length(x)
endfunction

function var = varianza(x);
  xbar = promedio(x);
  var = 0;
  for i = 1:length(x)
    var = var + (x(i)-xbar).^2;
  endfor
  var = var./(length(x)-1)
endfunction

function sum = sumaValores(x);
  sum = 0;
  for i = 1:length(x)
    sum = sum + x(i);
    endfor
endfunction

function sum = sumaCuadrados(x);
    sum = 0;
  for i = 1:length(x)
    sum = sum + x(i).*x(i);
  endfor
endfunction

function var = varianzaMejorada(x);
  xc = sumaCuadrados(x);
  xs = sumaValores(x);
  var = (xc-((xs.^2))./length(x))./(length(x)-1);
endfunction

% Parte (a)
x = [10000000000, 10000000001, 10000000002]
 varianzaA = varianza(x)
 % varianza = 1

 % Parte (b)
 varianzaM = varianzaMejorada(x)
 % varianza = 32768 
 
 % El resultado obtenido es radicalmente distinto al obtenido utilizando
 % el primer algoritmo. Esta diferencia se debe al tamano de los numeros,
 % dado que los numeros que manejamos son muy grandes y les estamos produciendo
 % solo pequenas diferencias entonces el computador no reconoce esta diferencia
 % y la redondea dejandola igual en el caso del segundo algoritmo. En el caso 
 % del primer algoritmo, como estamos realizando las operaciones sobre cantidades
 % grandes y luego es que las estamos juntando estamos acomulando menos error
 % que en el primer caso que es lo que nos esta produciendo un valor mas cercano
 % al real. El resultado exacto de la varianza es 1 (lo que nos da el primer
 % algoritmo) y lo podemos calcular a mano: 
 %note que el promedio de los tres valores es x(2), entonces 
 % (x(1)-x(2))^2 = 1, (x(2)-x(2))^2 = 0, (x(3)-x(2))^2 = 1, luego la suma de todos
 % es 2 que dividido por 2 (n-1) da 1.