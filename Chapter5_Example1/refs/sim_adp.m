%simulation_mm
clc
clear all
global A1 A2 B1 B2 K1 K2 B22 B12 B11 B21 Ef1 Ef2 angle10 angle20 ts
ts=15;
D1=1;
T1=2;
B12=.055;
B11=.01;
B22=0.01;
w0=100*pi;
H1=2;
angle10=2;
angle20=1;
D2=1;
T2=1;
B12=B12;
H2=1;
A1=[0     1         0;
    0  -D1/2/H1     w0/2/H1;
    0    0         -1/T1];
B1=[0;
    0;
    -1/T1];
Ef1=2;
Ef2=3;
K1=-[1 1 10];
Q=300*eye(3);
%K1=lqr(A1,B1,Q,1);
A2=[0     1         0;
    0  -D2/2/H2     w0/2/H2;
    0    0         -1/T2];
B2=[0;
    0;
    -1/T2];
K2=lqr(A2,B2,0.5*[1  0 0;0 1 0;0 0 1],10);

x0=[0.1;0;0];z0=[2;-1;0];
%xsave=[];
%zsave=[];


%%%%%%
%%%%%%% For ADP Learning
h=0.01;dt=h;                                                         % stepsize
n=1000; floor(50/h);                                      % step numbers, 100 is the final time                                   
t=0;                                                      % initial time                         
y(:,1)=[x0;z0;zeros(9+3,1)];                                                   % initial value                         
Ixx=[];
Ixu=[];
dlxx=[];


for j=1:30
    j1=(j-1)*10+1;
    j2=j*10;
for ii=j1:j2          
  %ODE solver (with fixed step size)---------- 
   t(ii+1)=h*ii;
   k1=sm_sys_adp(t(ii),y(:,ii));                      %dx=closeloopsys(t,x) is the function
   k2=sm_sys_adp(t(ii)+h/2,y(:,ii)+h*k1/2);
   k3=sm_sys_adp(t(ii)+h/2,y(:,ii)+h*k2/2);
   k4=sm_sys_adp(t(ii)+h,y(:,ii)+h*k3);
   y(:,ii+1)=y(:,ii)+h*(k1+2*k2+2*k3+k4)/6;      
end

if j>=20
Ixx=[Ixx,
     y(7:15,j2)'-y(7:15,j1)'];
Ixu=[Ixu,
     y(16:18,j2)'-y(16:18,j1)'];

dlxx=[dlxx;
      kron(y(1:3,j2),y(1:3,j2))'-kron(y(1:3,j1),y(1:3,j1))'];
end
end
dlxx=[dlxx(:,[1,2,3,5,6,9])]
%%


for k=1:10
  Qk=Q+K1'*K1;
  thetak=[dlxx, -2*Ixx*kron(eye(3),K1')-2*Ixu];
  psik=-Ixx*Qk(:);
  pv=inv(thetak'*thetak)*thetak'*psik;
  P=[pv(1)   pv(2)/2 pv(3)/2;
     pv(2)/2 pv(4)   pv(5)/2;
     pv(3)/2 pv(5)/2 pv(6)];
  L=pv(7:9);
  K1=L'
end



%% post learning
y=y(1:6,:);
h=0.05;dt=h; 
for ii=length(t):n          
  %ODE solver (with fixed step size)---------- 
   t(ii+1)=t(j2)+h*(ii-j2);
   k1=sm_sys_nl(t(ii),y(:,ii));                      %dx=closeloopsys(t,x) is the function
   k2=sm_sys_nl(t(ii)+h/2,y(:,ii)+h*k1/2);
   k3=sm_sys_nl(t(ii)+h/2,y(:,ii)+h*k2/2);
   k4=sm_sys_nl(t(ii)+h,y(:,ii)+h*k3);
   y(:,ii+1)=y(:,ii)+h*(k1+2*k2+2*k3+k4)/6;      
end

xsave=y(1:3,:);
zsave=y(4:6,:);


figure(2)
plot(t,xsave(1,:))


figure(3)
plot(t,zsave(1,:))


%
delta=xsave(1,:);%zeros(1,1000)+1;
delta0=delta;%+.5*sin(.1*[1:length(delta)])
delta1=zsave(1,:);%-0.5;
delta10=delta1;%+.5*sin(.15*[1:length(delta)])+0.01*randn(1,length(delta));


for i=1:length(delta)
    i
make_generators(t(i),angle10,delta(i)+angle10,angle20,delta1(i)+angle20,1*(t(i)>=ts),1*(t(i)>=2)+1*(t(i)>=3))
filename=['c:\withadp\test',num2str(i)];
saveas(gcf,filename,'jpg');
end





   