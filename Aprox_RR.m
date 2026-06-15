% Parámetros de Entrada:
% f = función a calcular el cero. Diseñada en un script .m
% x_min
% x_max
% M
% cota
% Parámetros de Salida:
% R = vector R, conteniendo las raíces aproximadas
function R = Aprox_RR(f,x_min,x_max,M,cota)
if nargin < 5
disp('Debe ingresar f, x_min, x_max, M y cota');
return;
end
R = []; % Inicializo vector de raíces aproximadas
h = (x_max - x_min)/M; %iterador
for k = 2:M-1 %unica forma de la que puedo tener todos los valores de x_k es
%si arranca el k en 2

x_km1 = x_min + (k-2)*h; % x_{k-1}
x_k = x_min + (k-1)*h; % x_k
x_kp1 = x_min + k*h; % x_{k+1}
y_km1 = feval(f, x_km1); %evaluo las funciones en cada coso
y_k = feval(f, x_k);
y_kp1 = feval(f, x_kp1);
% Criterio (i): cambio de signo
if y_k * y_km1 < 0
R = [R, (x_km1 + x_k)/2]; % raíz estimada en el medio

% Criterio (ii): valor pequeño y cambio de concavidad
elseif abs(y_k) < cota && (y_k - y_km1)*(y_kp1 - y_k) < 0
R = [R, x_k]; % raíz estimada en x_k
end
end
end