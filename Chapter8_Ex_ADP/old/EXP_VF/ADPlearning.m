% ADPlearning
% start of one trial
%clear all
global Kadp
warning off
para;
%Kadp=K0;
x_save=[];t_save=[];
[xn,un]=size(B);%size of x
N=400; %length of the window, should be at least greater than xn^2
NN=4;  %max iteration times
T=.001;
dt=0.0005;%0.0001; % period to data recording

%X=[0.25,-.25,0,0,0,0,zeros(1,12+1+4)];
X=[0,-.25,0,0,0,0,zeros(1,12+1+4)];
%for it=0:0

Dxx=[];
Iq=[];
% Ixx=[];
Ixu=[];
Iuu=[];

for i=1:N
    [t,X]=ode_yu_lrn((i-1)*T,i*T,X(end,:)',dt);
    Dxx=[Dxx; kron(X(end,1:6),X(end,1:6))-kron(X(1,1:6),X(1,1:6))];
    %Ixx=[Ixx; X(end,6+1:6+36)-X(1,6+1:6+36)];
    Iq=[Iq; X(end,6+1)-X(1,6+1)];
   % Ixu=[Ixu; X(end,6+36+1:6+36+12)-X(1,6+36+1:6+36+12)];
    Ixu=[Ixu; X(end,6+1+1:6+1+12)-X(1,6+1+1:6+1+12)];
    Iuu=[Iuu; X(end,end-3:end)-X(1,end-3:end)];
    x_save=[x_save;X];
    t_save=[t_save;t'];
end


    Dxx=Dxx(:,[1:6,8:12,15:18,22:24,29:30,36]); % Only the distinct columns left % ADP learning
    Iuu=Iuu(:,[1,2,4]);
%% Learning    
para
    Q=K'*R*K+Q;
    %Y=-Ixx*Q(:);
    Y=-Iq;
    % obtain the left-hand side of the equation
    %X1=[Dxx,-2*Ixx*kron(eye(6),K'*R)-2*Ixu*kron(eye(6),R),-Iuu];
    %X1=[Dxx,-2*Ixx*kron(eye(6),K'*R)-2*Ixu*kron(eye(6),R)];
    X1=[Dxx,-2*Ixu, Iuu];
    %pp=inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
    pp=X1\Y;
    P=[pp(1)  pp(2)     pp(3)     pp(4)   pp(5)    pp(6) ;
        0     pp(7)     pp(8)     pp(9)   pp(10)   pp(11) ;
        0         0     pp(12)    pp(13)  pp(14)   pp(15) ;
        0         0        0      pp(16)  pp(17)   pp(18) ;
        0         0        0          0   pp(19)   pp(20) ;
        0         0        0          0   0        pp(21)];
    P=(P+P')/2
%     Kadp=[pp(22:2:32)';pp(23:2:33)']
%% simulate the rest part of the movement
para
[t,X]=ode_yu_lrn(t(end),1,X(end,:)',dt);
x_save=[x_save;X]; t_save=[t_save;t'];
Kadp=[pp(22:2:32)';pp(23:2:33)']


% end of one trial
