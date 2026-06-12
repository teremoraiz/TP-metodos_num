%% DATOS DEL SISTEMA %%
p=985;         %% kg/m3
A=30;          %% m2
a=0.01;        %% m2
rh=1.2*10^3;   %% N*s/m5
mp=10;         %% kg
k=2.57*10^3;   %% N/m
g=9.81;        %% m/s2
b=15;          %%N*s/m
h0=4;          %%m
h=10^-3;       %% resolucion
t= 0:h:10;
n=length(t);

%% --EJERCICIO 1-- %%
Xt= Runge_Kutta_o4(@edo1,A,mp,b,a,rh,p,g,h0,k, 2,[0,0], t, h);
figure (1)
plot(t, Xt(1,:), 'Color', [0.8 0.2 0.6], 'LineWidth', 1.5) ;
xlabel('Tiempo [s]')
ylabel('x(t) [m]')
title('Ejercicio 1 - Posición traslacional del sistema')
legend('x(t)')


%% -- EJERCICIO 2-- %%
xf_num=Xt(1,end) 
xf_analitico = (p*g*a*A*h0) / (A*k + p*g*a^2);
error_abs= abs(xf_analitico - xf_num);
fprintf('\n--- Ejercicio 2 ---\n')
fprintf('Posición final numérica:   %.6f m\n', xf_num)
fprintf('Posición final analítica:  %.6f m\n', xf_analitico)
fprintf('Error absoluto:            %.2e m\n', error_abs)


%% --EJERCICIO 3-- %%
Xt3=Runge_Kutta_o4(@edo3,A,mp,b,a,rh,p,g,h0,k, 3,[p*g*h0; 0; 0], t, h);
RTA=[0 -1/k 0; 0 0 -1]*Xt3;

figure(2)
subplot(2,1,1)
plot(t, RTA(1,:), 'Color', [0.7 0.1 0.9], 'LineWidth', 1.5)
xlabel('Tiempo [s]')
ylabel('x(t) [m]')
title('Ejercicio 3 - Posición traslacional (modelo de estados)')
legend('x(t)')

subplot(2,1,2)
plot(t, RTA(2,:), 'Color', [1.0 0.4 0.8], 'LineWidth', 1.5)
xlabel('Tiempo [s]')
ylabel('v(t) [m/s]')
title('Ejercicio 3 - Velocidad traslacional (modelo de estados)')
legend('v(t)')

%% --EJERCICIO 6--%%
fprintf('\n--- Ejercicio 6 ---\n')
[media_x, desvio_x, media_v, desvio_v] = Secante(t, RTA, n, xf_num)


