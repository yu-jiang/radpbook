function [t,Y]=ode_yu_NF(t0,tf,x0,dt)
h=dt;                                                         % stepsize                                      % step numbers, 100 is the final time                                   
t=t0:dt:tf;                                                      % initial time                         
y=x0;       
Y=[];
    for clock=t  
        %clock
       y=y+motor2DNF(clock,y)*h;      
       Y=[Y y];
    end 
Y=Y';
end