function SimResults = GADP_FT_main()
% Demo #2 for Global Adaptive Dynamic Programming for Continuous-time
% Nonlinear Systems, by Yu Jiang and Zhong-Ping
% Jiang, IEEE Transasctions on Automatic Control, 2015
% 
% This paper can be found at 
% 1. http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7063901
% 2. http://arxiv.org/pdf/1401.0020.pdf
%
% System requirements:
% - MATLAB (Manually Tested in MATLAB R2014b)
% - MATLAB Symbolic Toolbox
% - CVX (free to download at http://cvxr.com/cvx/)
%
% Copyright 2015 Yu Jiang 
% 
% Contact: yu.jiang@nyu.edu (Yu Jiang)

% Initialize Parameters
            % x1  x2  x1x1 x1x2 x2x2  x1^3  x1x1x2 x1x2x2   x2^3
Params.F = [0   0    0    1    0     -1     0      -1       0;
            1   2    0    0    0      0     0       0       0];
Params.G = [0 0.7; 0.6 0.7];
Params.Q = diag([1,1,0,0,0,0,0,0,0]);
Params.R = 0.001;

% Inital control gains, thees values are taken from literature
Params.K0 = [10.283 -13.769; -10.7 -3.805];

% Extend the gains to match the new basis
K = [Params.K0 zeros(2,7)];
K_old = K;

Params.Gl = [0 1; 1 1];       % Lower Bound of G
Params.Gu = [0 0.5; 0.5 0.5]; % Upper Bound of G

% Initialize the value function
p = LocalInitialValueControlPair(Params);
p_old = ones(size(p))*100; 

x0 =[1 -2]; % Initial Condition
X = x0;

% Policy Iteraton Parameters
IterMax = 7;        %Max iterations
T = 0.02;           %Length of interval for data collection
NumIntervals = 200; %Number of intervals for one interation 
tol_conv = 0.0001;   %Convergence criterion

% Compute the objective function for SOSs in Poilcy Iteration
x1min = -0.1;
x1max =  0.1;
x2min = -0.1;
x2max =  0.1;
c = LocalComputeObjectiveFcn(x1min, x1max, x2min, x2max);

% This varibles go to the output
psave=[];
ksave= K;
Xsave=[];
tsave=[];

for j = 1:IterMax
    % Online Simulation for Data Collection
    [Phi, Xi, Theta, Xsave, tsave, t, X] = LocalOnlineDataCollection(T, ...
        NumIntervals, X, Xsave, tsave, j, Params, K);
    
    % Online Policy Iteration
    if norm(p - p_old)>tol_conv
        p_old = p;
        [p, K] = LocalOnlinePolicyIteratoin(Theta, Xi, Phi, p, c);
        numAIter = j;
    end
    
    psave = [psave;p_old'];    %#ok<AGROW>
    ksave = [ksave;K];         %#ok<AGROW>
end
% Note: Till this point, all the online simulation is finished.

SimResults.psave = psave;
SimResults.ksave = ksave;
SimResults.xsave = Xsave;
SimResults.tsave = tsave;

% Generate figures
SimResults.hFigs = LocalPostProcess(Params, ...
    t, Xsave, tsave, p, psave(1,:), K_old,K, numAIter);
end


%% -------------------------------------------------------------- 
% LocalOnlineDateCollection:
% Local function for simulation and online data collection
% ------------------------------------------------------------------------
function [Phi, Xi, Theta, Xsave, tsave, t, X] = LocalOnlineDataCollection(T, NumIntervals, X, Xsave, tsave, j, Params, K)
Phi=[]; Xi=[]; Theta=[]; % Matricies to collect online data
for i = 0:NumIntervals - 1
    [t,X] = ode45(@(t,x) FTSystemWrapper(t,x,Params,K), ...
        [i,(i+1)]*T+(j-1)*NumIntervals*T,...
        [X(end,1:2) zeros(1,35+9)]);
    Phi = [Phi;X(end,2+1:2+34+9)];
    Xi = [Xi;X(end,end)];
    Theta = [Theta; BasisQuadPhiX(X(end,1),X(end,2))'-BasisQuadPhiX(X(1,1),X(1,2))'];
    Xsave = [Xsave; X(:,1:2)];
    tsave = [tsave; t(:)];
end
end

%% -------------------------------------------------------------- 
% LocalOnlinePolicyIteratoin:
% Local function for implementing online SOS policy iteration.
% Note: This does not require the system dynamics
% ------------------------------------------------------------------------
function [p,K] = LocalOnlinePolicyIteratoin(Theta, Xi, Phi, p0, c)
cvx_begin sdp
% cvx_precision best
variable p(12,1)
variable P(5,5) symmetric
variable L(9,9) symmetric
% Obj: min integral{V(x)} on the set Omega
% The objective is equivalently converted to
% min c'*[x]_{1,5}
minimize(c'*p)

% 1) Equality constraint:
% Given p (V(x)), L and K can be uniquelly determined
LandK = (Phi'*Phi)\(Phi'*(-Xi-Theta*p));

l = LandK(1:25);
K = [LandK(26:34)'; LandK(35:43)'];

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

% 3) SOS constraint:
% V(x) <= V_old(x)
(p - p0) == [P(1,1) 
    P(2,1) + P(1,2) 
    P(2,2) 
    P(1,3) + P(3,1) 
    P(1,4) + P(4,1) + P(2,3) + P(3,2)
    P(1,5) + P(5,1) + P(2,4) + P(4,2)
    P(2,5) + P(5,2) 
    P(3,3)
    P(3,4) + P(4,3)
    P(3,5) + P(5,3)+P(4,4)
    P(4,5) + P(5,4)
    P(5,5)];
P <= 0;

cvx_end
end

function p = LocalInitialValueControlPair(Params)
K = [Params.K0 zeros(2,7)];

cvx_begin sdp
cvx_precision best
variable P(5,5) symmetric
variable L(9,9) symmetric
variable L1(9,9) symmetric
variable L2(9,9) symmetric

p = [P(1,1) P(2,1)+P(1,2) P(2,2) P(1,3)+P(3,1) P(1,4)+P(4,1)+P(2,3)+P(3,2) P(1,5)+P(5,1)+P(2,4)+P(4,2) P(2,5)+P(5,2) P(3,3) P(3,4)+P(4,3) P(3,5)+P(5,3)+P(4,4) P(4,5)+P(5,4) P(5,5)]';
W = [2*p(1) p(2) 3*p(4) 2*p(5) p(6) 4*p(8) 3*p(9) 2*p(10) p(11);
    p(2)   2*p(3) p(5) 2*p(6) 3*p(7) p(9) 2*p(10) 3*p(11) 4*p(12)];
P >= 0;

% Gl
H1 = (1/2*W'*(Params.F+Params.Gl*K)+1/2*(Params.F+Params.Gl*K)'*W)+Params.Q+K'*Params.R*K;
L1(1,1)==H1(1,1);
L1(1,2)+L1(2,1)==H1(1,2)+H1(2,1);
L1(2,2)==H1(2,2);
L1(1,3)+L1(3,1)==H1(1,3)+H1(3,1);
L1(1,4)+L1(4,1)+L1(2,3)+L1(3,2)==H1(1,4)+H1(4,1)+H1(2,3)+H1(3,2);
L1(1,5)+L1(5,1)+L1(2,4)+L1(4,2)==H1(1,5)+H1(5,1)+H1(2,4)+H1(4,2);
L1(2,5)+L1(5,2)==H1(2,5)+H1(5,2);
L1(1,6)+L1(6,1)+L1(3,3)== H1(1,6)+H1(6,1)+H1(3,3);
L1(1,7)+L1(7,1)+L1(2,6)+L1(6,2)+L1(3,4)+L1(4,3)==H1(1,7)+H1(7,1)+H1(2,6)+H1(6,2)+H1(3,4)+H1(4,3);
L1(1,8)+L1(8,1)+L1(2,7)+L1(7,2)+L1(3,5)+L1(5,3)+L1(4,4)==H1(1,8)+H1(8,1)+H1(2,7)+H1(7,2)+H1(3,5)+H1(5,3)+H1(4,4);
L1(1,9)+L1(9,1)+L1(2,8)+L1(8,2)+L1(5,4)+L1(4,5)==H1(1,9)+H1(9,1)+H1(2,8)+H1(8,2)+H1(5,4)+H1(4,5);
L1(2,9)+L1(9,2)+L1(5,5)==H1(2,9)+H1(9,2)+H1(5,5);
L1(3,6)+L1(6,3)==H1(3,6)+H1(6,3);
L1(3,7)+L1(7,3)+L1(4,6)+L1(6,4)==H1(3,7)+H1(7,3)+H1(4,6)+H1(6,4);
L1(3,8)+L1(8,3)+L1(4,7)+L1(7,4)+L1(5,6)+L1(6,5)==H1(3,8)+H1(8,3)+H1(4,7)+H1(7,4)+H1(5,6)+H1(6,5);
L1(3,9)+L1(9,3)+L1(4,8)+L1(8,4)+L1(5,7)+L1(7,5)==H1(3,9)+H1(9,3)+H1(4,8)+H1(8,4)+H1(5,7)+H1(7,5);
L1(4,9)+L1(9,4)+L1(5,8)+L1(8,5)==H1(4,9)+H1(9,4)+H1(5,8)+H1(8,5);
L1(5,9)+L1(9,5)==H1(5,9)+H1(9,5);
L1(6,6)==H1(6,6);
L1(6,7)+L1(7,6)==H1(6,7)+H1(7,6);
L1(6,8)+L1(8,6)+L1(7,7)==H1(6,8)+H1(8,6)+H1(7,7);
L1(6,9)+L1(9,6)+L1(7,8)+L1(8,7)==H1(6,9)+H1(9,6)+H1(7,8)+H1(8,7);
L1(7,9)+L1(9,7)+L1(8,8)==H1(7,9)+H1(9,7)+H1(8,8);
L1(9,8)+L1(8,9)==H1(9,8)+H1(8,9);
L1(9,9)==H1(9,9);
L1<=0;

% Gu
H2 = (1/2*W'*(Params.F+Params.Gu*K)+1/2*(Params.F+Params.Gu*K)'*W)+Params.Q+K'*Params.R*K;
L2(1,1)==H2(1,1);
L2(1,2)+L2(2,1)==H2(1,2)+H2(2,1);
L2(2,2)==H2(2,2);
L2(1,3)+L2(3,1)==H2(1,3)+H2(3,1);
L2(1,4)+L2(4,1)+L2(2,3)+L2(3,2)==H2(1,4)+H2(4,1)+H2(2,3)+H2(3,2);
L2(1,5)+L2(5,1)+L2(2,4)+L2(4,2)==H2(1,5)+H2(5,1)+H2(2,4)+H2(4,2);
L2(2,5)+L2(5,2)==H2(2,5)+H2(5,2);
L2(1,6)+L2(6,1)+L2(3,3)== H2(1,6)+H2(6,1)+H2(3,3);
L2(1,7)+L2(7,1)+L2(2,6)+L2(6,2)+L2(3,4)+L2(4,3)==H2(1,7)+H2(7,1)+H2(2,6)+H2(6,2)+H2(3,4)+H2(4,3);
L2(1,8)+L2(8,1)+L2(2,7)+L2(7,2)+L2(3,5)+L2(5,3)+L2(4,4)==H2(1,8)+H2(8,1)+H2(2,7)+H2(7,2)+H2(3,5)+H2(5,3)+H2(4,4);
L2(1,9)+L2(9,1)+L2(2,8)+L2(8,2)+L2(5,4)+L2(4,5)==H2(1,9)+H2(9,1)+H2(2,8)+H2(8,2)+H2(5,4)+H2(4,5);
L2(2,9)+L2(9,2)+L2(5,5)==H2(2,9)+H2(9,2)+H2(5,5);
L2(3,6)+L2(6,3)==H2(3,6)+H2(6,3);
L2(3,7)+L2(7,3)+L2(4,6)+L2(6,4)==H2(3,7)+H2(7,3)+H2(4,6)+H2(6,4);
L2(3,8)+L2(8,3)+L2(4,7)+L2(7,4)+L2(5,6)+L2(6,5)==H2(3,8)+H2(8,3)+H2(4,7)+H2(7,4)+H2(5,6)+H2(6,5);
L2(3,9)+L2(9,3)+L2(4,8)+L2(8,4)+L2(5,7)+L2(7,5)==H2(3,9)+H2(9,3)+H2(4,8)+H2(8,4)+H2(5,7)+H2(7,5);
L2(4,9)+L2(9,4)+L2(5,8)+L2(8,5)==H2(4,9)+H2(9,4)+H2(5,8)+H2(8,5);
L2(5,9)+L2(9,5)==H2(5,9)+H2(9,5);
L2(6,6)==H2(6,6);
L2(6,7)+L2(7,6)==H2(6,7)+H2(7,6);
L2(6,8)+L2(8,6)+L2(7,7)==H2(6,8)+H2(8,6)+H2(7,7);
L2(6,9)+L2(9,6)+L2(7,8)+L2(8,7)==H2(6,9)+H2(9,6)+H2(7,8)+H2(8,7);
L2(7,9)+L2(9,7)+L2(8,8)==H2(7,9)+H2(9,7)+H2(8,8);
L2(9,8)+L2(8,9)==H2(9,8)+H2(8,9);
L2(9,9)==H2(9,9);
L2<=0;

cvx_end
end

%% -------------------------------------------------------------- 
% LocalComputeObjectiveFcn: - Compute Coefficients for the objective
% functioin in SOSp PI
% ------------------------------------------------------------------------
function c = LocalComputeObjectiveFcn(x1min, x1max, x2min, x2max)
% In the SOS-base Policy Iteration, we need to solve an SOSp in each
% iteraton. The objective of the SOSp is
%
% min integration{V(x)}_Omega
%
% where Omega is a compact set, an area of interested of the system
% performance.

syms x1 x2
v_basis_fcn = BasisQuadPhiX(x1,x2);
c = double(int(int(v_basis_fcn,x1min,x1max),x2min,x2max));
end


%% --------------------------------------------------------------
% LocalPostProcess - Post process data and plot the resuts
% ------------------------------------------------------------------------
function hFigs = LocalPostProcess(Params, t, Xsave, tsave, p, p1, K_old,Knew, numAIter)
close all

[t,x] = ode45(@(t,x)FTSys(x, Params, Knew),[t(end) 30], Xsave(end,:));
tsave = [tsave;t(:)];
Xsave = [Xsave; x];

x00   = Xsave(end,:) + [5,10];         % coordinate transform
[t,x] = ode45(@(t,x)FTSys(x, Params, Knew),[30 35], x00);
tsave = [tsave;t(:)];
Xsave = [Xsave; x];

[t,x] = ode45(@(t,x)FTSys(x, Params, K_old),[30 35], x00);

% Figure 1. Plot state trajectories
h1 = figure(1);
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

% Figure 2. Plot and compare the value functions
h2 = figure(2);
x1 = (-10:1:10)/100;
x2 = (-10:1:10)/100;
vn = zeros(length(x1),length(x2));
v1 = zeros(length(x1),length(x2));
vs = []; us = []; un = [];
kn = vn;
k1 = v1;
for i=1:length(x1)
    for j=1:length(x2)
        vn(i,j) = p(:)'*BasisQuadPhiX(x1(i),x2(j));
        v1(i,j) = p1(:)'*BasisQuadPhiX(x1(i),x2(j));
        k1(i,j) = norm([K_old(1,:)*BasisPlanMono(x1(i),x2(j)), ...
                  K_old(2,:)*BasisPlanMono(x1(i),x2(j))]);
        kn(i,j) = norm([Knew(1,:)*BasisPlanMono(x1(i),x2(j)), ...
                  Knew(2,:)*BasisPlanMono(x1(i),x2(j))]);
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
    [0.845238095238095 0.731440045897881], ...
    'TextEdgeColor','none','FontSize',12,...
    'String',{'V_0(x1,x2)'});
annotation(gcf,'textarrow',[0.132142857142857 0.159949345986154],...
    [0.140476190476191 0.257882960413087], ...
    'TextEdgeColor','none','FontSize',12,...
    'String',{sprintf('V_%d(x1,x2)',numAIter)});

% Figure 3. Plot the control curve
h3 = figure(3);
surf(x1,x2,kn')
hold on
surf(x1,x2,k1')
hold off
xlabel('x_1')
ylabel('x_2')
annotation(gcf,'textarrow',[0.216071428571429 0.174535137214669],...
    [0.845238095238095 0.731440045897881], ...
    'TextEdgeColor','none','FontSize',12,...
    'String',{'|u_1|^2'});
annotation(gcf,'textarrow',[0.132142857142857 0.159949345986154],...
    [0.140476190476191 0.257882960413087], ...
    'TextEdgeColor','none','FontSize',12,...
    'String',{sprintf('|u_%d|^2',numAIter)});
%export_fig Ex2_control -pdf -transparent
hFigs = [h1;h2;h3];
end

%% Basic Minomial Functions are defined in separated files
% Basis for system dynamics and controller
% dx = F*sigma + G*u
% u = K*sigma;
% Basis for V(x)


