%para.m

m   = 1;  % mass 
b   = 10; % viscosity
c1  = 0.15/2;
c2  = 0.05/2;
dt  = 0.005;
tau = 0.05;

A =[0 0    1    0;
    0 0    0    1;
    0 0 -b/m    0;
    0 0    0 -b/m];
A1=A;
A1(3,1)=200;


B=[0   0;
   0   0;
   1/m 0;
   0 1/m];
B1=[1/tau 0;
    0     1/tau];

A12=[A1        B;
    zeros(2,4) -B1];
B12=[zeros(4,2);
    B1];


%Q=10*eye(4);
Q0=diag([100,100,.1,.1]);
Q1=diag([1000,100,.1,.1]);
R1=0.01*eye(2);
%clc
K0=lqr(A,B,Q0,R1);

Q1=Q1*10;
R1=R1*10;

K1=lqr(A,B,Q1,R1);


Q2=0.1*eye(2);
R2=0.01*eye(2);
% 
% Q2=Q2*10;
% R2=R2*10;

K2=lqr(zeros(2),B1,Q2,R2);

%K1
%K
%eig(A1-B*K1)'
%eig(A12-B12*K)'


%R=0.01*eye(2);

K=K2*[K1,eye(2)];


%disp(['check if it is very close to zero:', num2str(norm(P*A+A'*P-P*B*inv(R+D1'*P*D1+D2'*P*D2)*B'*P+Q))])





