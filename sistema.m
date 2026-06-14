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

%% --EJERCICO 4--%%
v_num = derivada_numerica(Xt(1,:), h);
v_ref = Xt(2,:)';% Comparación con RK4 del ejercicio 1
fprintf('Error máximo = %.4e m/s\n', max(abs(v_num - v_ref)));
figure(3)
plot(t, v_ref, '--k', 'LineWidth', 2,   'DisplayName', 'v(t) RK4 (ref)')
hold on
plot(t, v_num,  'Color', [1 0.41 0.71], 'LineWidth', 1.2, 'DisplayName', 'v(t) dif. centradas')
xlabel('Tiempo [s]'); ylabel('v(t) [m/s]')
title('Ejercicio 4 - Velocidad por derivada numérica')
legend('Location','best'); grid on

%% --EJERCICIO 5-- %%
t0 = t(1);
tF = t(end);
x_t = Xt(1,:);   % posicion x(t) -> Ejercicio 1
v_t = RTA(2,:);  % velocidad v(t) -> Ejercicio 3
% interpolamos los vectores para tener una funcion porque la necesitamos
% como parametro
fx_x = Spline_Cubica(t, x_t);
fv_v = Spline_Cubica(t, v_t);
%Funciones |x(t)|^2 y |v(t)|^2 evaluadas para pasarlas a tropezoidal y simpson ---
fx2 = @(tq) Eval_Spline(t, fx_x, tq).^2;
fv2 = @(tq) Eval_Spline(t, fv_v, tq).^2;
M = 1000;  % cantidad de intervalos
Ix_trap = Regla_Trapezoidal_Compuesta(fx2, t0, tF, M);
Iv_trap = Regla_Trapezoidal_Compuesta(fv2, t0, tF, M);
Ix_simp = Regla_Simpson_Compuesta(fx2, t0, tF, M);
Iv_simp = Regla_Simpson_Compuesta(fv2, t0, tF, M);
%ahora calculamos las potencias
Px_trap = Ix_trap / (tF - t0);
Pv_trap = Iv_trap / (tF - t0);
Px_simp = Ix_simp / (tF - t0);
Pv_simp = Iv_simp / (tF - t0);

fprintf('\n--- Ejercicio 5 ---\n')
fprintf('Px (Trapezoidal) = %.6f\n', Px_trap)
fprintf('Px (Simpson)     = %.6f\n', Px_simp)
fprintf('Pv (Trapezoidal) = %.6f\n', Pv_trap)
fprintf('Pv (Simpson)     = %.6f\n', Pv_simp)

%% --EJERCICIO 6--%%
%fprintf('\n--- Ejercicio 6 ---\n')
%[media_x, desvio_x, media_v, desvio_v] = Secante(t, RTA, n, xf_num)


