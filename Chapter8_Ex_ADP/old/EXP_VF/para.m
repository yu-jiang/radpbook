%para.m

tau=0.05;
m1=2;
m2=2;
kai=.7;
d11=13*kai;
d12=-18*kai;
d21=18*kai;
d22=13*kai;

c1=.15/2;
c2=.05/2;
dt=0.005;

in=0;

A0=[0 0 1        0      0      0;
    0 0 0        1      0      0;
    0 0 0        0     1/m1    0;
    0 0 0        0      0      1/m2;
    0 0 0        0      -1/tau 0;
    0 0 0        0      0      -1/tau];


%% simulating the velocity force field
A=A0+[0 0 0     0      0 0;
    0 0 0     0      0 0;
    0 0 d11/m1 d12/m1      0 0;
    0 0 d21/m2 d22/m2      0 0;
    0 0 0     0      0 0;
    0 0 0     0      0 0];


% A=A0+[0 0 0     0      0 0;
%       0 0 0     0      0 0;
%       150/m1    0 0     0      0 0;
%       0 0 0     0      0 0;
%       0 0 0     0      0 0;
%       0 0 0     0      0 0];

B=[0 0;
    0 0;
    0 0;
    0 0;
    1/tau 0;
    0 1/tau];

%%
Q0=[500 0; 0 1000];R=diag([.01,.01]);
theta=0*15/180*pi;%0.9453;
Qc = Q0;
Q = blkdiag(Qc,0.01*Qc,0.00005*Qc);


theta1=-36/180*pi;%0.9453;
TM1=[cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
%Qc1 = 250*log(1+30)*[cos(theta1)*cos(theta1) cos(theta1)*sin(theta1); ...
%    cos(theta1)*sin(theta1) sin(theta1)*sin(theta1)]+Q0;
Qc1 = 1500*[cos(theta1)*cos(theta1) cos(theta1)*sin(theta1); ...
    cos(theta1)*sin(theta1) sin(theta1)*sin(theta1)]+Q0;
Q1 = blkdiag(Qc1,0.01*Qc1,0.00005*Qc1);


%K0=[1000 0 250 0 10 0;0 1000 0 250 0 10];
%K0(1:2,1:2)=K0(1:2,1:2)*1.1;
K=lqr(A0,B,Q,R);
K0=K*3;
%K0(1:2,1:2)=K(1:2,1:2)*3;
%K0(1:2,1:2)=K(1:2,1:2)*2;
Ko=lqr(A,B,Q1,R);
% K0=K;
% K0(1,1)=K(1,1)+500;%230;
% K0(2,2)=K(2,2)+0;


%K=K1
% eig(A-B*K);

