%% --declaro la edo de orden 2 del ejercicio 1-- %%
function dydt = edo1 (t,Xt,A,mp,b,a,rh,p,g,h0,k)
        
        c1=p*g*a*h0/mp;
        c2=-(b+(a^2*rh))/mp;
        c3=- ((A*k)+(p*g*a^2))/(A*mp);
        dydt=zeros(2,1);

        dydt(1)=Xt(2);
        dydt(2)=c1+c2*Xt(2)+c3*Xt(1);
       
end
