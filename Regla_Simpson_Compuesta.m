function S = Regla_Simpson_Compuesta(f,a,b,M)
h = (b-a)/(2*M);
x = a:h:b;
Sum = 0;
for k=1:M
    Sum = Sum + feval(f,x(2*k-1)) + 4*feval(f,x(2*k)) + feval(f,x(2*k+1));
end
S = Sum*h/3;
end