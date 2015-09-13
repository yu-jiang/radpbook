function [t,Y]=ode_yu_lrn(t0,tf,x0,dt,Kadp)
h=dt;                                                         % stepsize                                      % step numbers, 100 is the final time                                   
t=t0:dt:tf;                                                      % initial time                         
y=x0;       
Y=[];
    for clock=t  
        clock;
       y=y+motor_adp(clock,y,dt,Kadp)*h;      
       Y=[Y y];
    end 
Y=Y';
end

function dX=motor_adp(t,X,dt, Kadp)
para;

%global Kadp%
%Kadp=Ko;

x=X(1:3); 
hx=x+randn(3,1)*ww; % measurement noise

u=-Kadp*hx;

w=randn*sqrt(dt);

%M=[c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];
M=c*u;

v=M*w; % control dependent noise

dx=A*x+B*u+B*v./dt; %6
dIr=hx'*Q*hx+u'*R*u;%kron(x,x); %1
dIxu=kron(x,R*v)./dt; % 1 %12;
dIuu=kron(v,v)./dt;         % 4
dX=[dx;   %   3
    dIr; %    1
    dIxu; %   3
    dIuu];%   1
end