% ADPlearning
% start of one trial
%clear all
function Kadp=ADPlearning(Kadp,iter_max)
clc
disp('[Simulating the online learning process]')
%global Kadp
%para;
%Kadp=Ko;

disp(['[0-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
for num_iter=1:iter_max
    
    x_save=[];t_save=[];
    N=400; %length of the window, should be at least greater than xn^2
    NN=4;  %max iteration times
    T=.001;
    dt=0.0001;
    
    %X=[0,-.25,0,0,0,0,zeros(1,12+1+4)]; % initial conditions
    X=[-.25,0,0,zeros(1,5)];
    
    Dxx=[];
    Iq=[];
    Ixu=[];
    Iuu=[];
    
    for i=1:N
        [t,X]=ode_yu_lrn((i-1)*T,i*T,X(end,:)',dt,Kadp);
        Dxx=[Dxx; kron(X(end,1:3),X(end,1:3))-kron(X(1,1:3),X(1,1:3))];
        Iq=[Iq; X(end,3+1)-X(1,3+1)];
        Ixu=[Ixu; X(end,3+1+1:3+1+3)-X(1,3+1+1:3+1+3)];
        Iuu=[Iuu; X(end,end)-X(1,end)];
        x_save=[x_save;X];
        t_save=[t_save;t'];
    end
    
    %%
    Dxx=Dxx(:,[1,2,3,5,6,9]);%[1:6,8:12,15:18,22:24,29:30,36]); % Only the distinct columns left % ADP learning
    
    %% Learning
    Y=-Iq;
    X1=[Dxx,-2*Ixu, Iuu];
    pp=inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
    P=[pp(1)  pp(2)     pp(3);
        0     pp(4)     pp(5);
        0         0     pp(6)];
    
    P=(P+P')/2;
    Kadp=[pp(7:9)'];
    disp(['[',num2str(num_iter), '-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
    
end
end


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
m=3;
A=[    0    1          0;
       0  -b/m        1/m;
       0     0     -1/tau];
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