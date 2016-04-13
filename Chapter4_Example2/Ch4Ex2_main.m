function SimResults = Ch4Ex2_main()
%% GADP for a scalar polynomial system
% Demo #1 for Global Adaptive Dynamic Programming for Continuous-time
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

%% Setting parameters
% Parameters for Simulation the dynamic system
SysParams.Q = 0.01*[1 0 0;0 1 0;0 0 0];  % The cost function will be 
                                         % ([x]_{1,3})'*Q*[x]_{1,3}+u^2
SysParams.noiseFlag = true;   % Inicator for noise on/off
SysParams.K = [0.1 0 0.01];   % Initial feeback gain matrix
SysParams.F = [0 0.01 0];     % System dynamics dx = F*[x]_1,3 + u, 
                              % F is not used in Online Policy Iteration
SysParamsInit = SysParams;    % Make a copy of initial parameters
                            
% Parameters for learning
Params.iter_tol = 0.007;
Params.T = 0.1;          % Length of each time interval for data collection
Params.MaxIter = 10;     % Maximum iteration numbers 
Params.IterInterval = 5; % Length of each time interval for learning
Params.xinit = 2;        % Initial condition for the actual state x


% Other params      
x = [Params.xinit;zeros(9,1)]'; % Initial condition for the augmented x
P = eye(2)*10;                % Initialize V_0: ([x]_{1,2})'*P*[x]_{1,2}
Pold = -100*eye(2);           % Initialize V_{-1}

% Simulation results to export
SimResults.Ksave = SysParams.K;  % Keep track of the feedback gains
SimResults.Psave = P(:)';        % Keep track of the value function
SimResults.Usave = [];           % Actual control signal during simulation
SimResults.Xsave=[];             % Actual x(t) during simulation
SimResults.Tsave=[];             % Sample time points during simulation


% caculate the weights for V
c = LocalComputeObjective(-1, 1);

%% Start online simulation
for i = 0:Params.MaxIter-1
   % Data collection
   Theta = [];Sigma = []; Xi = [];  % Data matrices for online learning
    for j = 0:Params.IterInterval/Params.T-1
        [t,x] = ode45(@(t,x) SystemWrapper(t,x,SysParams), ...
            [j,j+1]*Params.T + i*Params.IterInterval, ...
            [x(end,1) zeros(1,9)]);
        Theta = [Theta; 
                (x(end,1)^2-x(1,1)^2) (x(end,1)^3-x(1,1)^3) (x(end,1)^4-x(1,1)^4)];
        Sigma = [Sigma; 
                 x(end,2:6)  -x(end,7:9)];
        Xi = [Xi; x(end,10)];
        SimResults.Tsave = [SimResults.Tsave;t];
        SimResults.Xsave = [SimResults.Xsave;x(:,1)];
        for k=1:length(t)
            SimResults.Usave = [SimResults.Usave; 
               LocalComputeControlSignal(x(k,1),SysParams.K,...
               SysParams.noiseFlag,t(k))];
        end
    end
    
    % SOS-based Online Policy iteration
    if norm(P(:)-Pold(:)) > Params.iter_tol 
        % Calling local function for policy improvement and policy
        % evaluation. Notice that we do not pass in the systen dynamics,
        % i.e, the F vector. Because online learning does not rely on it.
        [Pn,K] = LocalOnlinePI(Sigma,Xi,Theta,... % Data matrices collected online
                               c,  ... % weights
                               P);     % previous value function
        Pold = P; P = Pn; % Save the current and the old P
        SysParams.noiseFlag = true;
        SysParams.K = K;
        SimResults.Ksave = [SimResults.Ksave; K(:)'];
        SimResults.Psave = [SimResults.Psave; P(:)'];
        % Qlave=[Qlave;dQl(:)'];

    else
        SysParams.noiseFlag = false;
        disp('Convergence has been attained ...update is not necessary')
    end
end

%% Post-process and plot results
SimResults.hFigs = LocalPostProcess(SysParams, SysParamsInit, ...
   SimResults,Params, P);
end

%% LocalPostProcess: Process results and generate figures
function hFigs = LocalPostProcess(SysParams, SysParamsInit,SimResults, Params, P)
% Figure 1: 
% Comparison of x between GADP and unlearned system
hFig1 = figure(1);
[t0,y0] = ode45(@(t,x) SystemWrapper(t,x,SysParamsInit),...
    [0 50], ...
    [Params.xinit, zeros(1,9)]);
y0 = y0(:,1); % Only need the first column for the actual x
for i=1:length(t0)
 u0(i) = SysParamsInit.K * y0(i).^[1 2 3]'; % Unleared controller
end

plot(SimResults.Tsave, SimResults.Xsave, 'b-', ...  % Learned
     t0,y0, 'r-.', ...                              % Unlearned
     'linewidth', 2)                         
axis([0 50 -.5 2])
myLegend = legend('With GADP-based controller', 'With initial controller');
set(myLegend, 'Fontsize', 12);
xlabel('time (sec)', 'FontSize', 12)
%ylabel('x', 'FontSize', 12)

% Create textarrows
annotation(hFig1,'textarrow', ...
    [0.281132075471698 0.226415094339623],...
    [0.845386533665835 0.800498753117207],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'1st iteration'});
annotation(hFig1,'textarrow',...
    [0.443396226415094 0.44188679245283],...
    [0.244389027431421 0.309127182044887],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'4th iteration'});
annotation(hFig1,'textarrow',...
    [0.369811320754717 0.372452830188679],...
    [0.448877805486284 0.386334164588527],'TextEdgeColor','none',...
    'FontSize',12,...
    'String',{'3rd iteration'});
annotation(hFig1,'textarrow',...
    [0.284905660377358 0.286037735849057],...
    [0.321695760598504 0.416408977556109],'TextEdgeColor','none', ...
    'FontSize',12,...
    'String',{'2nd iteration'});

% Figure 2
hFig2 = figure(2);
plot(SimResults.Tsave, SimResults.Usave,'Linewidth',2)
myLegend = legend('u');
set(myLegend, 'FontSize', 12);
xlabel('time (sec)', 'FontSize', 12)
axis([0 50 -0.5 2]);
% Create textarrows
annotation(hFig2 ,'textarrow',[0.239092495636998 0.2107082880569],...
    [0.221748400852878 0.320754616656652],'TextEdgeColor','none', ...
    'FontSize',12,...
    'String',{'1st iteration'});
annotation(hFig2 ,'textarrow',[0.396160558464223 0.361981626000198],...
    [0.176972281449893 0.27332776800004],'TextEdgeColor','none', ...
    'FontSize',12,...
    'String',{'3rd iteration'});
annotation(hFig2 ,'textarrow',[0.331588132635253 0.2999993414337],...
    [0.439232409381663 0.335385523398326],'TextEdgeColor','none', ...
    'FontSize',12,...
    'String',{'2nd iteration'});
annotation(hFig2 ,'textarrow',[0.471204188481675 0.443631993150911],...
    [0.37953091684435 0.304862789720793],'TextEdgeColor','none', ...
    'FontSize',12,...
    'String',{'4th iteration'});



% Figures 3: 
% Comparing the initial, the improved, and the ideal value fcn 
syms v(y) 
F = SysParams.F;
% Solve the HJB analytically
vsx = dsolve(diff(v)*(F(1)*y+F(2)*y^2+F(3)*y^3) + 0.01*(y^2+y^4)-1/4*(diff(v))^2==0, v(0)==0);
vsx = vsx(1);
x = -2.5:.01:2.5;
vn = []; % V_n
v1 = []; % V_1   
vs = []; % Optimal V*
us = []; % Optimal u*
u1 = []; % Initial u_1
un = []; % Improved u_n
P1 = [SimResults.Psave(2,1) SimResults.Psave(2,2);
      SimResults.Psave(2,3) SimResults.Psave(2,4)] ;

for y = x 
vn = [vn [y y^2]*P*[y ;y^2]];
v1 = [v1 [y y^2]*P1*[y ;y^2]];
vs = [vs eval(vsx)];
u1 = [u1 -1/2*SysParamsInit.K*[y;y^2;y^3]];
un = [un -1/2*SysParams.K'*[y;y^2;y^3]];
us = [us -1/2*eval(diff(vsx))];
%(y*(y*(101*y^2 + 100)^(1/2) + 101*y^2 + 100))/(50*(101*y^2 + 100)^(1/2))
end

hFig3 = figure(3);
plot(x,v1,'g:',x,vn,'r-.',x,vs,'b','linewidth',2)
myLegend = legend('V_1 : Initial cost', 'V_4: Improved cost', ...
    'V^o: Optimal cost');
set(myLegend, 'FontSize', 12);
xlabel('x', 'FontSize', 12)

% Figures 4: 
% Comparing the initial, the improved, and the ideal control input 
hFig4 = figure(4);
plot(x,u1,'g:',x,un,'r-.',x,us,'b','linewidth',2)
myLegend = legend('u_1: Initial control policy', ...
    'u_4: Improved control policy', ...
    'u^o: Optimal control policy');
set(myLegend, 'FontSize', 12);
xlabel('x', 'FontSize', 12)

% Export all figure handles
hFigs = [hFig1; hFig2; hFig3; hFig4];

end

%% LocalOnlinePI
% Local function to implement the online ADP method.
% Note1: CVX solver is required (Download: http://cvxr.com/cvx/)
% Note2: The learning process does not depend on the system dynamics (F)
function [Pn,K] = LocalOnlinePI(Sigma,Xi,Theta,c,P)
cvx_begin sdp
variable pv(3,1)
variable dQl(3,3) symmetric
% SDP Objective function
minimize(c(1)*pv(1)+c(3)*pv(3))
% 1) Equality constraint
LnK = -inv(Sigma'*Sigma)*Sigma'*(Xi + Theta*pv(:));
% 2) SOS constraint:
%    L*[x]_2,6 is SOS
% i.e., there exists Ql>0, such that
%    L*[x]_2,6 = ([x]_1,3)'*Ql*[x]_1,3
LnK(1) == dQl(1,1);                        %#ok<*EQEFF> CVX Syntax
LnK(2) == dQl(1,2) + dQl(2,1);
LnK(3) == dQl(1,3) + dQl(3,1) + dQl(2,2);
LnK(4) == dQl(3,2) + dQl(2,3);
LnK(5) == dQl(3,3);
dQl >= 0;                                  %#ok<*VUNUS> CVX Syntax
% 3) SOS constraint
%    pv_old*[x]_1,3 - pv*[x]_1,3 is SOS
% This implies that
% V_old >= V_new (i.e.,the value function is reduced)
Pn = [ pv(1) 1/2*(pv(2)); 1/2*(pv(2)) pv(3)];
Pn <= P;
K = LnK(6:8);
cvx_end
end


 
%% LocalComputeObjective
% Compute the objective function for the SOSp in
% Policy Iterations. The objective depends on the interval [x_min, x_max].
% This requires MATLAB Symbolic Toolbox.
function c = LocalComputeObjective(x_min, x_max)
syms z
c = double(int(z.^[2,3,4], x_min,x_max));
end

%% LocalComputeControlSignal
% Compute the control input
function u = LocalComputeControlSignal(x,K,noiseFlag,t)
 u = -1/2*K(:)'*x.^[1 2 3]'+ ExplorationNoise(t)*noiseFlag;
 u = abs(u);
end