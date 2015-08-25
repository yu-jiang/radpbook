function dX=motor2DAL(t,X)
global K dt 

para;

x=X(1:6);

% The gain can also be obtained using SimuVF30.m
K30=[410.20   37.53  100.13  -17.50    1.65   -0.06
    -176.24  313.48    3.99   97.01   -0.06    1.62];
u=-K30*x;

w=randn(2,1)*sqrt(dt);

%M=[c1*u(1) c2*u(1); -c2*u(2) c1*u(2)];

M=[c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];

v=M*w; % control dependent noise


dx=A*x+B*u+B*v./dt; %6
dIq=x'*(Q1+K'*R*K)*x; %1
dxxi=kron(x',v'*R)'./dt; %12;
duu=kron(v,v)./dt;         % 4
dX=[dx;  % \dot x
    dIq; % reward
    dxxi;
    duu];


end