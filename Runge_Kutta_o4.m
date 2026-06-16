function X = Runge_Kutta_o4 (f,A,mp,b,a,rh,p,g,h0,k, m,X0, t, h)
    
    n=length (t);
    Xt=zeros(m,n);
    Xt(:,1) = X0;

   
    for i= 1:n-1
        k1=feval(f,t(i),Xt(:,i),A,mp,b,a,rh,p,g,h0,k);
        k2=feval(f,t(i)+h/2,Xt(:,i)+h/2*k1,A,mp,b,a,rh,p,g,h0,k);
        k3=feval(f,t(i)+h/2,Xt(:,i)+h/2*k2,A,mp,b,a,rh,p,g,h0,k);
        k4=feval(f,t(i)+h,Xt(:,i)+h*k3,A,mp,b,a,rh,p,g,h0,k);

        Xt(:,i+1)=Xt(:,i)+ h/6* (k1+2*k2+2*k3+k4);  
    end
    X= Xt;
    
   
end
  
%%%% poner en el informe

