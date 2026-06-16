function Mat_S=Spline_Cubica(X,Y)
%X: variables independientes en forma de vector.
%Y: variables dependientes en forma de vector
%ambos vectores o fila o columna y del mismo tamanio.
% devuelve matriz de Nx4 con los coeficientes de la cubica (Ak,Bk,Ck,Dk)
% Sk(x) = Ak(x - Xk)^3 + Bk(x - Xk)^2 + Ck(x - Xk) + Dk

N=length(X)-1;
hk=diff(X); %distancia en el eje horizontal entre un punto y el siguiente
dk=diff(Y)./hk; %pendiente de la recta que une un punto con el siguiente
m0=0; mN=0;

a = hk(2:N-1);                        % diagonal inferior
b = 2*(hk(1:N-1) + hk(2:N));         % diagonal principal (el doble de la suma de intervalos vecinos)
c = hk(2:N-1);                        % diagonal superior

uk = 6 * diff(dk);                    % término independiente
B = uk;
B(1)   = uk(1)   - hk(1)*m0;         % aplica condición de borde izquierda
B(N-1) = uk(N-1) - hk(N-1)*mN;       % aplica condición de borde derecha


mk=zeros(N+1,1); 
mk(1)=m0; mk(N+1)=mN;
N1=N-1;

%eliminacion gaussiana
for k=2:N1
    Piv=a(k-1)/b(k-1);
    b(k)=b(k)-Piv.*c(k-1);
    B(k)=B(k)-Piv*B(k-1);
end
X2=zeros(N1,1);
X2(N1)=B(N1)/b(N1); % al terminar este for, el sistema quedó triangular superior

for k=N1-1:-1:1
    X2(k)=(B(k)-c(k)*X2(k+1))/b(k);
end
mk(2:N)=X2;
Mat_S=zeros(N,4);

for k=0:N-1
    Mat_S(k+1,1)=(mk(k+2)-mk(k+1))/(6*hk(k+1));
    Mat_S(k+1,2)=mk(k+1)/2;
    Mat_S(k+1,3)=dk(k+1)-hk(k+1)*(2*mk(k+1)+mk(k+2))/6;
    Mat_S(k+1,4)=Y(k+1);
end
end