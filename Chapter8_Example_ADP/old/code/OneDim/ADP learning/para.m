m=3;  % mass
m0=1;
b=10; % viscosity
c  = 0.1;
tau = 0.05;
ww=0.005;% measurment noise


A=[    0    1          0;
    0  -b/m        1/m;
    0     0     -1/tau];
A0=[    0    1          0;
    0  -b/m0        1/m0;
    0     0     -1/tau];
B=[    0;
    0;
    1/tau];

Q=diag([400,1,.001]);

%dt=0.001;


R=0.01;

[Ko,Po,Eo]=lqr(A0,B,Q,R);
% [K1,P1,E1]=lqr(A,B,Q,R);