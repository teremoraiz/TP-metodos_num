function v = derivada_numerica(x, h)
% Entradas:
%   x  : vector de la señal (posición x(t))  
%   h  : paso de tiempo Δt [s]
% Salida:
%   v  : derivada numérica (velocidad v(t))  [mismo tamaño que x]

x = x(:);     % uso la columna
N = length(x);
v = zeros(N, 1);
v(1)     = (x(2) - x(1)) / h; % diferencia hacia adelante 
v(2:N-1) = (x(3:N) - x(1:N-2)) / (2*h);%diferencias centradas  
v(N)     = (x(N) - x(N-1)) / h;%diferencia hacia atrás

end