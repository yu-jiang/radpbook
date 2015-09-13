function [t,Y]=ode2(t0,tf,x0)
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
m0=1;
m=3;  % mass
b=10; % viscosity
c  = 0.1;
tau = 0.05;
A0=[    0    1          0;
       0  -b/m0        1/m0;
       0     0     -1/tau];
A=[    0    1          0;
       0  -b/m        1/m;
       0     0     -1/tau];
B=[    0;
       0;
   1/tau];

Q=diag([400,1,.001]);
%Q=diag([400,0,0]);


R=0.01;

[Ko,Po,Eo]=lqr(A0,B,Q,R);
%%%   End Para

x=x(:);

%hx=x+randn(3,1)*ww;

u=-Ko*x;


w=randn*sqrt(dt);

M=c*u;

v=M*w; % control dependent noise

dx=A*x+B*u+B*v./dt; %6

end