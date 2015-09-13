function [t,Y]=ode_nl(t0,tf,x0)
dt=0.005;;                                                         % stepsize                                      % step numbers, 100 is the final time
t=t0:dt:tf;                                                      % initial time
y=x0;
Y=[];
for clock=t
    %clock
    y=y+motor(clock,y,dt)*dt;
    Y=[Y y];
end
Y=Y';
end


function dx=motor(t,x,dt)
%global K dt

%%% Start Para
para;
%%%   End Para

x=x(:);

u=-Ko*x;

w=randn*sqrt(dt);

M=c*u;

v=M*w; % control dependent noise

dx=A*x+B*u+B*v./dt; %6

end