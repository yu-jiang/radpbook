function dX=motor2DNF(t,X)
%global K dt 

para;

x=X(1:6);

u=-K*x;

w=randn(2,1)*sqrt(dt);

%M=[c1*u(1) c2*u(1); -c2*u(2) c1*u(2)];

M=[c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];

v=M*w; % control dependent noise


dx=A0*x+B*u+B*v./dt; %6
dIq=x'*(Q+K'*R*K)*x; %1
dxxi=kron(x',v'*R)'./dt; %12;
duu=kron(v,v)./dt;         % 4
dX=[dx;  % \dot x
    dIq; % reward
    dxxi;
    duu];


end