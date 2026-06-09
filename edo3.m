%% -declaro matriz de ME de ejercicio 3-- %%
function X= edo3 (t,Xt,A,mp,b,a,rh,p,g,h0,k)
    
    mat_a=[0 0 (p*g*a)/A; 0 0 k; -a/mp -1/mp -(b+(a^2*rh))/mp];
    X=mat_a*Xt;

end