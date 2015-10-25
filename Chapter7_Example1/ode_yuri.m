% This is a self-written ode45 solver with fixed steps
function [t,Y]=ode_yuri(t0,tf,x0,dt)
h=dt;                                                   % stepsize                                      % step numbers, 100 is the final time                                   
t=t0:dt:tf;                                             % time                         
y=x0;       
Y=[];
% initial value                         
    for clock=t           
      %ODE solver (with fixed step size)---------- 
       k1=polysys(clock,y);                      %dx=closeloopsys(t,x) is the function
       k2=polysys(clock+h/2,y+h*k1/2);
       k3=polysys(clock+h/2,y+h*k2/2);
       k4=polysys(clock+h,y+h*k3);
       y=y+h*(k1+2*k2+2*k3+k4)/6;      
       Y=[Y y];
    end 
Y=Y';
end