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
T= 0:h:10;
n=length(T);

%% --EJERCICIO 1-- %%
fprintf('\n--- Ejercicio 1 ---\n')
Xt= Runge_Kutta_o4(@edo1,A,mp,b,a,rh,p,g,h0,k, 2,[0,0], T, h);
figure (1)
plot(T, Xt(1,:), 'Color', [0.8 0.2 0.6], 'LineWidth', 1.5) ;
xlabel('Tiempo [s]')
ylabel('x(t) [m]')
title('Ejercicio 1 - Posición traslacional del sistema')
legend('x(t)')


%% -- EJERCICIO 2-- %%
fprintf('\n--- Ejercicio 2 ---\n')
xf_analitico = (p*g*a*A*h0) / (A*k + p*g*a^2);
error_abs= abs(xf_analitico - Xt(1,end));
fprintf('Posición final numérica:   %.6f m\n', Xt(1,end))
fprintf('Posición final analítica:  %.6f m\n', xf_analitico)
fprintf('Error absoluto:            %.2e m\n', error_abs)


%% --EJERCICIO 3-- %%
fprintf('\n--- Ejercicio 3 ---\n')
Xt3=Runge_Kutta_o4(@edo3,A,mp,b,a,rh,p,g,h0,k, 3,[p*g*h0; 0; 0], T, h);
RTA=[0 -1/k 0; 0 0 -1]*Xt3;

figure(2)
subplot(2,1,1)
plot(T, RTA(1,:), 'Color', [0.7 0.1 0.9], 'LineWidth', 1.5)
xlabel('Tiempo [s]')
ylabel('x(t) [m]')
title('Ejercicio 3 - Posición traslacional (modelo de estados)')
legend('x(t)')

subplot(2,1,2)
plot(T, RTA(2,:), 'Color', [1.0 0.4 0.8], 'LineWidth', 1.5)
xlabel('Tiempo [s]')
ylabel('v(t) [m/s]')
title('Ejercicio 3 - Velocidad traslacional (modelo de estados)')
legend('v(t)')

%% Interpolacion de vectores x y v por spline cubica %%
x_t = Xt(1,:);   % posicion x(t)
v_t = Xt(2,:);  % velocidad v(t)
fx= Spline_Cubica(T, x_t);
fv= Spline_Cubica(T, v_t);
fx_ev = @(tq) Eval_Spline(T, fx, tq);
fv_ev = @(tq) Eval_Spline(T, fv, tq);

%% --EJERCICO 4--%%
v_num = zeros(n, 1);

% Derivada punto por punto usando Richardson
for i = 1:n
    [D, err, relerr, n] = Extrapolacion_Richardson_O2n(fx_ev, T(i), 1e-5, 1e-8);
    v_num(i) = D(end, end); 
end

% 3. Comparación con RK4 del ejercicio 1
fprintf('Error máximo = %.4e m/s\n', max(abs(v_num - (v_t)')));

% 4. Gráficos
figure(3)
plot(T, v_t, '--k', 'LineWidth', 2,   'DisplayName', 'v(t) RK4 (ref)')
hold on
plot(T, v_num,  'Color', [1 0.41 0.71], 'LineWidth', 1.2, 'DisplayName', 'v(t) dif. centradas')
xlabel('Tiempo [s]'); ylabel('v(t) [m/s]')
title('Ejercicio 4 - Velocidad por derivada numérica')
legend('Location','best'); grid on

%% --EJERCICIO 5-- %%
fprintf('\n--- Ejercicio 5 ---\n')
t0 = T(1);
tF = T(end);

%Funciones |x(t)|^2 y |v(t)|^2 
fx5 = @(tq) fx_ev(tq).^2; 
fv5 = @(tq) fv_ev(tq).^2;
M = 1000;  % cantidad de intervalos
%Integrales
Ix_trap = Regla_Trapezoidal_Compuesta(fx5, t0, tF, M);
Iv_trap = Regla_Trapezoidal_Compuesta(fv5, t0, tF, M);
Ix_simp = Regla_Simpson_Compuesta(fx5, t0, tF, M);
Iv_simp = Regla_Simpson_Compuesta(fv5, t0, tF, M);
%Potencias
Px_trap = Ix_trap / (tF - t0);
Pv_trap = Iv_trap / (tF - t0);
Px_simp = Ix_simp / (tF - t0);
Pv_simp = Iv_simp / (tF - t0);

fprintf('Px (Trapezoidal) = %.6f\n', Px_trap)
fprintf('Px (Simpson)     = %.6f\n', Px_simp)
fprintf('Pv (Trapezoidal) = %.6f\n', Pv_trap)
fprintf('Pv (Simpson)     = %.6f\n', Pv_simp)

%% --EJERCICIO 6--%%
%valor medio y desvio estandar calculado con x(t)
xt_c= x_t - x_t(end);              %centramos en cero
fx6 = Spline_Cubica(T, xt_c );
fx6_ev = @(tq) Eval_Spline(T, fx6, tq);
cota = 1e-6; 

Rx= Aprox_RR(fx6_ev,t0, tF, M, cota );  %guardamos las raices de las pos
N= length(Rx);

for i= 1: N-1
    Tce_x(i)= Rx(i+1) - Rx(i);     %guardamos los tiempos consecutivos 
end


T_medio_expx = sum(Tce_x) / length(Tce_x)        
T_desvio_expx= sqrt(sum((Tce_x-T_medio_expx).^2)/(length(Tce_x)-1))


%valor medio y desvio estandar calculado con v(t)
Rv= Aprox_RR(fv_ev,t0, tF, M, cota );  %guardamos las raices de la vel
N= length(Rv);

for i= 1: N-1
    Tce_v(i)= Rv(i+1) - Rv(i);     %guardamos los tiempos consecutivos 
end

T_medio_expv = sum(Tce_v) / length(Tce_v)        %formula de consigna
T_desvio_expv= sqrt(sum((Tce_v-T_medio_expv).^2)/(length(Tce_v)-1))


%calculo teorico valor medio
Wn= sqrt(((A*k)+(p*g*(a*a)))/(A*mp));
eps= (1/(2*Wn))*((b+(a*a)*rh)/mp);
Wd= Wn*sqrt(1-(eps*eps));

T_medio_teo= pi/Wd 
errorTx= abs(T_medio_teo - T_medio_expx)

errorTv=abs(T_medio_teo - T_medio_expv)

