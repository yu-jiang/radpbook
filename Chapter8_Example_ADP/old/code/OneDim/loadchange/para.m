m=1;  % mass
b=10; % viscosity
c  = 0.1;
tau = 0.05;
ww=0;

A=[    0    1          0;
       0  -b/m        1/m;
       0     0     -1/tau];
B=[    0;
       0;
   1/tau];

Q=diag([400,1,.001]);



R=0.01;

[Ko,Po,Eo]=lqr(A,B,Q,R);
%[Ko,Po,Eo]=lqr(A0,B,Q,R);