function RADP_MultiMachine()
%This is the Matlab code used for the paper:
%
% Yu Jiang and Zhong-Ping Jiang, "Robust Adaptive Dynamic Programming for
% Large-Scale Systems with an Application to Multimachine Power Systems,"
% IEEE Transactions on Circuits and Systems II: Express Briefs, vol. 59, no.
% 10, pp. 693-697, 2012.
%
% The code is free for everyone to use. Please cite the above paper in your
% publication if you do use the code.
%
% Please contact yu.jiang@nyu.edu if you find
% any bugs or have any suggestions on improving the code. Thanks!
%
clc
clear all
close all
global Kadp
warning off
mm_para %load the parameters

% simulate the system until 4s
% operating in steady-state from 0s to 1s
disp('Simulating the system on steady state to 1s...')
[t0,y0] = ode45(@Local_mmsys_online_radp,[0,1],zeros(15*(Nm-1),1));
% add an impulse disturbance at 1s
disp('Adding an impulse disturbance at 1s and simulating to 4s...')
[t1,y1]=ode45(@Local_mmsys_online_radp,[1,4],y0(end,:)'-[kron(pm(2:end),[0,0,1]),zeros(1,12*(Nm-1))]');


y=y1(end,:);
N=20; %# of learning intervals

Ixx=zeros(N,9,Nm-1);Ixu=zeros(N,3,Nm-1);Dxx=zeros(N,6,Nm-1);

yy=[y0;y1];tt=[t0;t1];
%%
for i=0:N-1
    disp(['simulating the ', num2str(i+1),'-th interval...', num2str(N-i), 'left']);
    % simulate the trajectories for learning
    [t,y]=ode45(@Local_mmsys_online_radp,[4+i/N,4+(i+1)/N],y(end,:));
    for j=1:Nm-1
        Ixx(i+1,:,j)=y(end,imxx(j):imxx(j)+8)-y(1,imxx(j):imxx(j)+8);%y(1,imuu(j):imuu(j)+2)];
        Ixu(i+1,:,j)=y(end,imuu(j):imuu(j)+2)-y(1,imuu(j):imuu(j)+2);
        Dxx(i+1,:,j)=[y(end,1+(j-1)*3)^2-y(1,1+(j-1)*3)^2
            y(end,1+(j-1)*3)*y(end,2+(j-1)*3)-y(1,1+(j-1)*3)*y(1,2+(j-1)*3)
            y(end,1+(j-1)*3)*y(end,3+(j-1)*3)-y(1,1+(j-1)*3)*y(1,3+(j-1)*3)
            y(end,2+(j-1)*3)^2-y(1,2+(j-1)*3)^2
            y(end,2+(j-1)*3)*y(end,3+(j-1)*3)-y(1,2+(j-1)*3)*y(1,3+(j-1)*3)
            y(end,3+(j-1)*3)^2-y(1,3+(j-1)*3)^2]';
    end
    yy=[yy;y];  tt=[tt;t];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% off-Policy learning for all Generators       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
K=zeros(1,3,Nm-1);
for i=1:Nm-1
    K(:,:,i)=[10   50   0]; % initial feedback gain matrix
end


for j=1:Nm-1
    Kprev=[100 100 100];it=0;
    while (norm(K(:,:,j)-Kprev)>0.01)
        it=it+1;
        Kprev=K(:,:,j);
        Qk=1000*eye(3)+K(:,:,j)'*K(:,:,j);
        Theta=[Dxx(:,:,j) -2*Ixx(:,:,j)*kron(eye(3),K(:,:,j)')-2*Ixu(:,:,j)*kron(eye(3),1)];
        Psi=-Ixx(:,:,j)*Qk(:);
        % pv=inv(Theta'*Theta)*Theta'*Psi;
        pv=pinv(Theta)*Psi;
        K(:,:,j)=pv(end-2:end)';
    end
    Kadp(:,:,j)=K(:,:,j);
    disp(['The' num2str(j+1) '-th machine stopped learning after' num2str(it) 'iterations'])
end

%%
[t,y]=ode45(@Local_mmsys_online_radp,[5,15],y(end,:)); % simulate the post-learning performance
y=[yy;y];t=[tt;t];

% plot the angles
myfigure
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local_mmsys_online_radp -- Local function to simulat the dynamics of the
% multimachine power system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxx = Local_mmsys_online_radp(t,xx)
global Kadp
mm_para %load the parameters

x=[0;0;0];
for i=2:Nm
    id=(i-2)*3+1:(i-2)*3+3;
    x(:,i)=xx(id);
end


K=zeros(1,3,Nm-1);
for i=1:Nm-1
    K(:,:,i)=[10   50   0];
end

% calculate angular differences in matrix form
dlt=zeros(Nm);d=zeros(1,Nm);
for i=1:Nm
    dlt(i,:)=x(1,i)+dlt0(i)-x(1,:)-dlt0';
    d(i)= E(i)*(E.*(x(2,i)-x(2,:)))*(BX(i,:).*cos(dlt(i,:))-GX(i,:).*sin(dlt(i,:)))';
end


u=0;
for i=2:Nm
    if t>=4 & t<=5 % learning stage is between 4s and 5s
        u(i)=-K(:,:,i-1)*x(:,i)+0.001*sin(100*t);
    else
        u(i)=-K(:,:,i-1)*x(:,i);
    end
end

if t>5 % update the control policies and start the post-learning stage
    for i=2:Nm
        u(i)=-Kadp(:,:,i-1)*x(:,i);
    end
end

for i=2:Nm
    dx(:,i-1)=A(:,:,i)*x(:,i)+B(:,:,i)*(u(i)-d(i));
end

dIxxu=zeros(12,Nm-1);
if t>=4 & t<=5
    for i=2:Nm
     dIxxu(:,i-1)=[kron(x(:,i),x(:,i));kron(x(:,i),u(i)-d(i))];
    end
end
dxx = [dx(:);
       dIxxu(:)];
end

