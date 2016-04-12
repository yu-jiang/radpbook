function SimResults = Ch4Ex4_main()
%% Global Adaptive Dynamic Programming for an inverted pendulum
%
% System requirements:
%
% - MATLAB (Manually Tested in MATLAB R2014b)
%
% - MATLAB Symbolic Toolbox
%
% - CVX (free to download at http://cvxr.com/cvx/)
%
%    You can download CVX along with the examples, unzip tools.zip to 
%    the root level folder of the examples. Then, run the following
%
%     >> run('.\tools\cvx-w64\cvx\cvx_setup.m')
%
%     >> run('.\tools\cvx-w64\cvx\cvx_startup.m')
%
%
% Copyright 2015 Yu Jiang 
% 
% Contact: yu.jiang@nyu.edu (Yu Jiang)

SimResults = [];
k = 1;   
m = 1;   % Mass of the pendulum
l = 1;   % Length of the pendulum
g = 9.8; % Gravity accelaration rate
Params.F = [0    1     0  ;   0   -k*l/m g ];
Params.G = [0; 1/m];
Params.Q = diag([10 10 0]);   % Weighting matrix
Params.R = 1;                 % Weighting matrix
Params.Noise = 1;
x = [-1.5 1];       % Initial Contidition;
Params.x0 = x;

K = [-10 -1 -15];   % Initial Feedback gains
Psave = [];
Ksave = K;
xsave = [];
tsave = [];

T = 0.01;            % Length of time interval for data collection
P_old = 1000*eye(2); % Initializing the previous Value function
P = 0.9*P_old;       % Initialize the current Value function
c = [1 0 1];         % Coefficient for SOS policy iteration

NumDataIntv  = 50;
Iter = 0;
Tol = 0.1;
while norm(P_old - P)> Tol
    P_old = P;
    % Online Simulation for Data collection
    Theta = []; Xi = []; Phi = [];  % Data matrices
    for IterIntv = 0:NumDataIntv - 1
        [t,x] = ode45(@(t,x) LocalInvertedPendulumSysWrapper(t,x,K,Params), ...
            ([IterIntv, IterIntv + 1] + NumDataIntv * Iter) * T,...
            [x(end,1:2) zeros(1,10)]);
        Theta = [Theta; x(end,1)^2-x(1,1)^2 x(end,1)*x(end,2)-x(1,1)*x(1,2) x(end,2)^2-x(1,2)^2];
        Phi = [Phi; x(end,3:8) -2*x(end,9:11)];
        Xi = [Xi; x(end,12)]; 
        xsave = [xsave; x(:,1:2)];
        tsave = [tsave;t];
    end
    
    % Online SOSp based policy iteration
    [P,K] = LocalOnlinePolicyIteration(Theta,Xi,Phi,P_old,c);
    
    % Save results
    Psave = [Psave; P(:)'];
    Ksave = [Ksave; K(:)'];
    Iter = Iter + 1;
end
Params.Noise = 0;

% Post learning simulation
[t,x] = ode45(@(t,x) LocalInvertedPendulumSysWrapper(t,x,K,Params), ...
    [tsave(end) 5],[x(end,1:2) zeros(1,10)]');
xsave = [xsave; x(:,1:2)];
tsave = [tsave;t];

% Create Simulation output data
SimResults.xsave = xsave;
SimResults.tsave = tsave;
SimResults.Psave = Psave;
SimResults.Ksave = Ksave;
SimResults.P = P;
SimResults.Iter = Iter;

% Plot results
LocalPostProcessData(Params,SimResults);
end

%% LocalOnlinePolicyIteration: Implement online SOS policy iteration
function [P,K] = LocalOnlinePolicyIteration(Theta,Xi,Phi,P_old,c)
cvx_begin sdp
variables p(3,1)
variables bt r1 r2  r3 r4  r5 r6
% Objective of the SDP
minimize(c*p)

% 1) Equality constraint to calculate L and K
LandK = pinv(Phi)*(Xi + Theta*p);

% 2) Inequality constraint
L = LandK(1:6);
dFGQR = [L(1)   L(2)/2 L(4)/2;
    L(2)/2 L(3)   L(5)/2;
    L(4)/2 L(5)/2 L(6)];
O = [0 0 bt;
    0 0 0;
    bt 0 0];
bt>=0;
gamma1 = [r1 0 r2/2; 0  0  0; r2/2 0 0];
gamma2 = [r3 0 0; 0  0  0; 0 0 r4];
gamma3 = [0 0 r5/2; 0  0  0; r5/2 0 r6];
r1 >= 0;
r3 >= 0;
r5 >= 0;
r1+r2 >= 0;
r3+r4 >= 0;
r5+r6 >= 0;
dFGQR + O + gamma1 + gamma2 + gamma3 <= 0;

% 3) Inequality constraint
P = [p(1) p(2)/2;
    p(2)/2 p(3)];
P_old - P >= 0;
cvx_end;

K = LandK(7:9)';
end

%% invertedPendulumSys: Inverted Pendulum System Dynamics
function dx = LocalInvertedPendulumSys(x,K,Params)
u = K * LocalSigma(x);
dx = Params.F * LocalSigma(x) + Params.G * u;
end

%% LocalInvertedPendulumSysWrapper: 
% Inverted Pendulum System Dynamics with external states for learning purpose
function dx = LocalInvertedPendulumSysWrapper(t,x,K,Params)
 x1 = x(1);
 x2 = x(2);
 sgm  =  [x1;x2;sin(x1)];

 e  = LocalExplNoise(t) * Params.Noise;
 u  = K * sgm+e; 
 dx = Params.F * sgm + Params.G * u;
 
 dZ  =  [x1*x1 x2*x1 x2*x2 x1*sin(x1) x2*sin(x1) sin(x1)^2]';
 deZ = sgm*e;
 
 Qk  = Params.Q + K'*Params.R*K;
 dQ  = sgm' * Qk * sgm; 
 dx  = [dx;  %2 
        dZ;  %6
        deZ; %3
        dQ;  %1
        ]; %12
end

%% LocalSigma: Basic function for the plant
function y = LocalSigma(x)
 y = [x(1) x(2) sin(x(1))]';
end

%% LocalExplNoise: Generate exploration noise
function e = LocalExplNoise(t)
e = sin(10*t) + sin(3*t) + sin(50*t) - sin(10*t) + sin(0.7*t) - sin(100*t);
end

%% LocalPostProcessData - Plot results after all simulation is finished
function LocalPostProcessData(Params,SimResults)
tsave = SimResults.tsave;
xsave = SimResults.xsave;
Ksave = SimResults.Ksave;
Psave = SimResults.Psave;
P = SimResults.P;
x1 = -2:0.2:2;
x2 = -5:0.5:5;
vn = zeros(length(x1),length(x2));
v1 = zeros(length(x1),length(x2));
P1 = [Psave(1,1) Psave(1,2);Psave(1,3) Psave(1,4)] ;
for i=1:length(x1)
    for j=1:length(x2)
        vn(i,j)=[x1(i) x2(j)]*P*[x1(i) x2(j)]';
        v1(i,j)=[x1(i) x2(j)]*P1*[x1(i) x2(j)]';
    end
end
figure(1)
surf(x1,x2,vn')
hold on
surf(x1,x2,v1')
hold off
xlabel('x_1')
ylabel('x_2')

K = Ksave(1,:);[t1,y1]=ode45(@(t,x) LocalInvertedPendulumSys(x,K,Params), ...
    [0 5],Params.x0);
%%
figure(2)
subplot(211)
plot(tsave,xsave(:,1),'b-',t1,y1(:,1),'r-.','Linewidth',2)
axis([0 5 -2 2])
ylabel('x_1')
xlabel('time (sec)')
legend('With GADP-based controller','With initial controller')
subplot(212)
plot(tsave,xsave(:,2),'b-',t1,y1(:,2),'r-.','Linewidth',2)
axis([0 5 -4 6])
legend('With GADP-based controller','With initial controller')
ylabel('x_2')
xlabel('time (sec)')
end