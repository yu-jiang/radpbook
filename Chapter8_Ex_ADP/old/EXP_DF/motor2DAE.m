function dX=motor2DAE(t,X)
global K dt 

para;

x=X(1:6);


K30=[831.94 12.53 96.08 1.55 1.61 0.02;
    15.83 330.97 1.62 71.57 0.02 1.36];
u=-K30*x;

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