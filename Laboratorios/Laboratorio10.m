% Laboratorio 10
% Parte 1:

% Definicion de la funcion:
y =@(x) sin(x.*x);
dy =@(x) 2.*x.*cos(x.*x);

% Definicion de intervalos
ie = linspace(-pi,pi, 17);
ic = chebyshev(17,-pi, pi);

% Evaluaciones de y:

ee = y(ie);
ec = y(ic);
dec = dy(ic);

% Polinomio de Hermite:
[hc, zc] = hermite(ic, ec', dec');
% Evaluaciones de los polinomios:
ii = linspace(-pi,pi, 500);
plot(ii', spline(ie',ee',ii'), "--")
hold on
plot(ii', hornerN(hc, ii', zc), ".-")
hold on 
plot(ii', y(ii'))
hold on
title("Interpolaciones a sin(x^2) mediante Hermite con 17 ptos de Chebyshev y Splines en 17 ptos equidistantes y de Chebyshev")
legend("Splines: 17 ptos equidistantes", "Hermite: 17 ptos de Chebyshev", "Funcion original")

% Observaciones:
% En la grafica se puede notar que la interpolacion hecha por el polinomio de hermite se ajusta mucho mas a la 
% grafica real que la interpolacion hecha con los esplines, en los puntos en los que la funcion oscila los splines
% se tienden a alejar mucho mas de la grafica real de el polinomio de hermite.

% Pregunta 2:
ptoMedE = NaN*zeros(16,1)
ptoMedC = NaN*zeros(16,1)
errorH1 = NaN*zeros(16,1)
errorHI = NaN*zeros(16,1)
errorS1 = NaN*zeros(16,1)
errorSI = NaN*zeros(16,1)
for i =1:16
  ptoMedE(i) = (ie(i)+ie(i+1))/2
  ptoMedC(i) = (ic(i)+ic(i+1))/2
endfor

errorH1 = norm(y(ptoMedC)-hornerN(hc, ptoMedC, zc),1)/norm(y(ptoMedC),1)
errorHI = norm(y(ptoMedC)-hornerN(hc, ptoMedC, zc),Inf)/norm(y(ptoMedC),Inf)
errorS1 = norm(y(ptoMedE)-spline(ie',ee',ptoMedE), 1)/norm(y(ptoMedE),1)
errorSI = norm(y(ptoMedE)-spline(ie',ee',ptoMedE), Inf)/norm(y(ptoMedE),Inf)

% Resultados:
% errorH1 =  0.00000018708
% errorHI =  0.00000092844 < 0.0000005 
% errorS1 =  0.12178
% errorSI =  0.32149 < 5e-1

% Esto corroba nuestras observaciones sobre lo que vimos en las graficas, el error
% que acomula el polinomio de Hermite es muchisimo menor que el de los splines y por
% lo tanto tiene mas cifras significaticas que los splines. En general, para pocos puntos
% la aproximacion con el polinomio interpolante de hermite es una mejor aproximacion
% que los splines, mucho mas aun si los puntos de prueba no son equidistantes si no de Chebyshev.

% Se tiene que el numero de las cifras significativas para el polinomio de hermite
% es (en norma infinito) 6 pues 0.00000092844 < 0.000005 = 5e-6 y para los splines
% es 1 pues  0.32149 < 0.5 = 5e-1
