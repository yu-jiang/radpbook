%para.m

m   = 2;%1;  % mass 
b   = 10; % viscosity
c1  = 0.15/1.5;
c2  = 0.05/1.5;
dt  = 0.005;
tau = 0.05;

A =[0 0    1    0;
    0 0    0    1;
    0 0 -b/m    0;
    0 0    0 -b/m];

A0=[0 0    1    0;
    0 0    0    1;
    0 0 -b/m    0;
    0 0    0 -b/m];


B=[0   0;
   0   0;
   1/m 0;
   0 1/m];
B1=[1/tau 0;
    0     1/tau];

A12=[A        B;
    zeros(2,4) -B1];
B12=[zeros(4,2);
    B1];


%Q=10*eye(4);
Q1=diag([1000*4,1000*4,10,10]); R1=0.1*eye(2);

%Q1=diag([1000,1000,25,25]); R1=0.01*eye(2); % after learning
Q1(1,2)=0;
Q1(2,1)=Q1(1,2);


K1=lqr(A,B,Q1,R1);

Q2=1*eye(2);
R2=0.1*eye(2);

K2=lqr(zeros(2),B1,Q2,R2);

R=0.01*eye(2);

K=K2*[K1,eye(2)];


%disp(['check if it is very close to zero:', num2str(norm(P*A+A'*P-P*B*inv(R+D1'*P*D1+D2'*P*D2)*B'*P+Q))])





