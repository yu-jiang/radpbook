function dX=motor2D_lrn(t,X)
para;

global Kadp

x=X(1:6);

u=-Kadp*x;%+100*sin(10*t)+100*sin(0.1*t);

w=randn(2,1)*sqrt(dt);

M=[c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];

v=M*w; % control dependent noise

dx=A*x+B*u+B*v./dt; %6
dIr=x'*Q1*x+u'*R*u;%kron(x,x); %1
dIxu=kron(x,R*v)./dt; % 1 %12;
dIuu=kron(v,v)./dt;         % 4
dX=[dx;   %   6
    dIr; %  1 36
    dIxu; %   12
    dIuu];%   4
end