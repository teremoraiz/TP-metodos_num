function S = Regla_Trapezoidal_Compuesta(f,a,b,M)
h = (b-a)/M;
Sum = 0;
x0 = a;
for k=1:(M-1)
    Sum = Sum + feval(f,x0+k*h);
end
S = h*(feval(f,a) + feval(f,b))/2 + h*Sum;
end