function [] = GADP_FT_main()
global K

%% Initialize Parameters
% x1  x2  x1x1 x1x2 x2x2  x1^3  x1x1x2 x1x2x2   x2^3
Params.F = [0   0    0    1    0     -1     0      -1       0;
    1   2    0    0    0      0     0       0       0];
Params.Q = diag([1,1,0,0,0,0,0,0,0]);
Params.R = 0.001;

% Inital control gains, thees values are taken from literature
Params.K0 = [10.283 -13.769; -10.7 -3.805];
% Extend the gains to match the new basis
K = [Params.K0 zeros(2,7)];

Params.Gl = [0 1; 1 1];       % Lower Bound of G
Params.Gu = [0 0.5; 0.5 0.5]; % Upper Bound of G

% Initialize the values function
p_old = InitialValueControlPair(Params);

Kold = K;
x0 =[1 -2]; % Initial Condition
X=x0;


% Compute the objective function for SOSs in Poilcy Iteration
x1min =-0.1;
x1max = 0.1;
x2min =-0.1;
x2max = 0.1;
c = ComputeObjectiveFcn(x1min, x1max, x2min, x2max);

clear x1 x2
psave=[];
ksave= K;
Xsave=[];
tsave=[];

% ji=[0 199 199 199 199 199 199 199 199 199]
% jisum=[0 200 400 600 800 1000 1200 1400 1600];% 100+100+100];

T=0.02; %0.005;
for j=1:7
    % Online Simulation for Data Collection
    Phi=[]; Xi=[]; Theta=[]; % Matricies to collect online data
    for i=0:200-1
        [t,X] = ode45(@jetsysonline, ...
            [i*T,(i+1)*T]+(j-1)*200*T,...
            [X(end,1:2) zeros(1,35+9)]);
        Phi = [Phi;X(end,2+1:2+34+9)];
        Xi = [Xi;X(end,end)];
        Theta = [Theta; bphi(X(end,1),X(end,2))'-bphi(X(1,1),X(1,2))'];
        Xsave = [Xsave; X(:,1:2)];
        tsave = [tsave; t(:)];
    end
    
    % Online Policy Iteration
    [p, K] = LocalOnlinePolicyIteratoin(Theta, Xi, Phi, p_old, c);
    
    psave = [psave;p(:)'];
    %ksave=[ksave;K];
    % P0=P;
    p_old = p;
    if j==1
        % P1=P;
    end
end


Knew = K;

[t,x]=ode45(@jetsys0,[t(end) 30], Xsave(end,:));
tsave=[tsave;t(:)];
Xsave=[Xsave; x];

x00 = Xsave(end,:) +[5,10];
[t,x]=ode45(@jetsys0,[30 35], x00);
tsave=[tsave;t(:)];
Xsave=[Xsave; x];

%% Draw the figures
%close all
K=[K0 zeros(2,7)];
[t,x]=ode45(@jetsys0,[30 35], x00);

%[t1,x1]=ode45(@jetsys0,[15 40], x0+[20,10]);
% t=[t;t1];
% x=[x;x1];
%% figure(1)
subplot(2,2,1)

plot(tsave,Xsave(:,1), 'b-','Linewidth', 2)
axis([0 35 -3 7])
legend('With GADP-based controller')
xlabel('time (sec)')
ylabel('x_1')

subplot(2,2,2)
plot(tsave,Xsave(:,1),'b-',t,x(:,1),'r-.', 'Linewidth', 2)
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)')


axis([29.8 30.5 -3 7])


subplot(2,2,3)
plot(tsave,Xsave(:,2),'b-','Linewidth', 2)
axis([0 35 -4 15])
legend('With GADP-based controller')
xlabel('time (sec)')
ylabel('x_2')


subplot(2,2,4)
plot(tsave,Xsave(:,2),'b-',t,x(:,2),'r-.', 'Linewidth', 2)
axis([29.8 30.5 -4 15])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)')



%subplot(313)
%plot(t,x*K0(:),'b-')%,t,x(:,2),'r:', 'Linewidth', 2)
%state_annotations
%saveas(gcf,'Ex2_state.eps', 'psc2');
%saveas(gcf,'Ex2_state.pdf');
%saveas(gcf,'Ex2_state.jpg');


%%
figure(2)
x1=-10:1:10;
x2=-10:1:10;
vn=zeros(length(x1),length(x2));
v1=zeros(length(x1),length(x2));
vs=[];
us=[];
un=[];
kn=vn;
k1=v1;
for i=1:length(x1)
    for j=1:length(x2)
        vn(i,j)=phi(x1(i),x2(j))'*P*phi(x1(i),x2(j));
        v1(i,j)=phi(x1(i),x2(j))'*P1*phi(x1(i),x2(j));
        k1(i,j)=norm([Kold(1,:)*sigma(x1(i),x2(j)),Kold(2,:)*sigma(x1(i),x2(j))]);
        kn(i,j)=norm([Knew(1,:)*sigma(x1(i),x2(j)),Knew(2,:)*sigma(x1(i),x2(j))]);
    end
end
surf(x1,x2,vn')
hold on
surf(x1,x2,v1')
hold off
xlabel('x_1', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
% Create axes
view(gca,[-40.5 14]);
annotation(gcf,'textarrow',[0.216071428571429 0.174535137214669],...
    [0.845238095238095 0.731440045897881],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_0(x1,x2,0)'});
annotation(gcf,'textarrow',[0.132142857142857 0.159949345986154],...
    [0.140476190476191 0.257882960413087],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_7(x1,x2,0)'});
%saveas(gcf,'Ex2_cost.eps', 'psc2');
%saveas(gcf,'Ex2_cost.pdf');


%%
figure(3)
surf(x1,x2,kn')
hold on
surf(x1,x2,k1')
hold off
xlabel('x_1')
ylabel('x_2')
annotation(gcf,'textarrow',[0.216071428571429 0.174535137214669],...
    [0.845238095238095 0.731440045897881],'TextEdgeColor','none','FontSize',12,...
    'String',{'|u_1|^2'});
annotation(gcf,'textarrow',[0.132142857142857 0.159949345986154],...
    [0.140476190476191 0.257882960413087],'TextEdgeColor','none','FontSize',12,...
    'String',{'|u_7|^2'});
%saveas(gcf,'Ex2_control_surf.eps', 'psc2');
%export_fig Ex2_control -pdf -transparent
%%
% figure(4)
% plot(1:8, psave(:,1),'r-o',1:8, psave(:,3), '-*', 1:8, psave(:,5), '-^', 'linewidth',2)
% axis([0.5,8.5,-2, 17])
% legend('p_{i,1} (x_1^2)', 'p_{i,2} (x_2^2)', 'p_{i,3} (x_1^2x_2)')
% xlabel('Iteration', 'FontSize', 12)
% set(gcf,'PaperPositionMode','auto')
%saveas(gcf,'Ex2_psave.pdf') %, 'psc2');
%saveas(gcf,'Ex2_psave.eps', 'psc2');
%
%
%
% %%
% % figure(2)
% % export_fig Ex2_cost -pdf -transparent
figure(3)
%export_fig Ex2_control -pdf -transparent
%figure(4)
%export_fig Ex2_psave -pdf -transparent
end

function [p,K] = LocalOnlinePolicyIteratoin(Theta, Xi, Phi, p0, c)
cvx_begin sdp
variable p(12,1)
variable P(5,5) symmetric
variable L(9,9) symmetric
% Obj: min integral{V(x)} on the set Omega
% The objective is equivalently converted to
% min c'*[x]_{1,5}
(p - p0) == [P(1,1) P(2,1)+P(1,2) P(2,2) P(1,3)+P(3,1) P(1,4)+P(4,1)+P(2,3)+P(3,2) ...
    P(1,5)+P(5,1)+P(2,4)+P(4,2) P(2,5)+P(5,2) P(3,3) P(3,4)+P(4,3) ...
    P(3,5)+P(5,3)+P(4,4) P(4,5)+P(5,4) P(5,5)]';
minimize(c'*p)

% 1) Equality constraint:
% Given p (V(x)), L and K can be uniquelly determined
LandK = (Phi'*Phi)\(Phi'*(-Xi-Theta*p));

%

%     W = [2*p(1) p(2) 3*p(4) 2*p(5) p(6) 4*p(8) 3*p(9) 2*p(10) p(11);
%         p(2)   2*p(3) p(5) 2*p(6) 3*p(7) p(9) 2*p(10) 3*p(11) 4*p(12)];
P <= 0;

l = LandK(1:25);
K = [LandK(26:34)'; LandK(35:43)']

% 2) SOS contraint:
% l'*[x] = -dV/dx (f+gu) - r(x,u) is SOS
l == [L(1,1);
    L(1,2)+L(2,1);
    L(2,2);
    L(1,3)+L(3,1);
    L(1,4)+L(4,1)+L(2,3)+L(3,2);
    L(1,5)+L(5,1)+L(2,4)+L(4,2);
    L(2,5)+L(5,2);
    L(1,6)+L(6,1)+L(3,3);
    L(1,7)+L(7,1)+L(2,6)+L(6,2)+L(3,4)+L(4,3);
    L(1,8)+L(8,1)+L(2,7)+L(7,2)+L(3,5)+L(5,3)+L(4,4);
    L(1,9)+L(9,1)+L(2,8)+L(8,2)+L(5,4)+L(4,5);
    L(2,9)+L(9,2)+L(5,5);
    L(3,6)+L(6,3);
    L(3,7)+L(7,3)+L(4,6)+L(6,4);
    L(3,8)+L(8,3)+L(4,7)+L(7,4)+L(5,6)+L(6,5);
    L(3,9)+L(9,3)+L(4,8)+L(8,4)+L(5,7)+L(7,5);
    L(4,9)+L(9,4)+L(5,8)+L(8,5);
    L(5,9)+L(9,5);
    L(6,6);
    L(6,7)+L(7,6);
    L(6,8)+L(8,6)+L(7,7);
    L(6,9)+L(9,6)+L(7,8)+L(8,7);
    L(7,9)+L(9,7)+L(8,8);
    L(9,8)+L(8,9);
    L(9,9)];
L>=0;

cvx_end
end
