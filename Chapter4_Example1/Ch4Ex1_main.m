function Ch4Ex1_main()
%% SOS-based Policy Iteration for a car suspension systems
% The function is tested in MATLAB R2014b
% Copyright 2015 Yu Jiang
% Contact Yu Jiang (yu.jiang@nyu.edu)

% System requirements:
% - MATLAB (Manually Tested in MATLAB R2014b)
% - MATLAB Symbolic Toolbox
% - SDPT3-4.0
% - SISOTOOLS (free to download at http://www.cds.caltech.edu/sostools/)

% You can download tools.zip and run setuptools.m in the folder.

syms x1 x2 x3 x4 real
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m
kn = 0.1*ks;

% State matrices
A = [ 0 1 0 0;
    [-ks -bs ks bs]/mb ; ...
    0 0 0 1;
    [ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 10000/mb ; 0 0 ; [kt -10000]/mw];
B = B(:,2);

% LQR feedback gains for the linearized system
% Klqr = lqr(A,B,eye(4),1);
% f(x)
f = A*[x1;x2;x3;x4] + [0;-kn*(x1-x3)^3/mb;0;kn*(x1-x3)^3/mw];

% Polynomial Weighting functions
q0 = 100*x1^2+x2^2+x3^2+x4^2;
rx = 1;%+(x1^2+x2^2+x3^2+x4^2);
vars = [x1;x2;x3;x4];

% Initialize the SOSp
prog = sosprogram(vars);

% The Lyapunov function V(x)
[prog,V] = sospolyvar(prog,monomials([x1;x2;x3;x4],2:4),'wscoeff');

% Objective of the SOSp
myObj = int(int(int(int(V,-.5,.5),-10,10),-.5,.5),-10,10);

% Add Inequality constraint to assure V is positive definite
prog = sosineq(prog,V-0.0001*(x1^2+x2^2+x3^2+x4^2));

% Add Inequality constraint to assure stability and performance
expr = -[diff(V,x1) diff(V,x2) diff(V,x3) diff(V,x4)]*rx*f-rx*q0;
prog = sosineq(prog,expr);

% Solve the SOSp
prog = sossolve(prog);

% Obtain the Initial Lyapunov function
V0 = sosgetsol(prog,V);
V_old = V0;

% Initializing the old contol policy
u_prev = zeros(size(x1));

% Iteration
for i=1:10
    clear prog V
    %------------------------------ SOSp Start ----------------------------
    prog = sosprogram(vars);
    [prog,V] = sospolyvar(prog,monomials([x1;x2;x3;x4],2:4),'wscoeff');
    prog = sosineq(prog,V_old - V);
    prog = sosineq(prog, V);
    u = -1/2*B'*[diff(V_old,x1) diff(V_old,x2) diff(V_old,x3) diff(V_old,x4)].';
    qfcn =rx*q0 + u'*u;
    expr = -[diff(V,x1) diff(V,x2) diff(V,x3) diff(V,x4)]*(rx*f+B*u)-qfcn;
    prog = sosineq(prog,expr);
    prog = sossetobj(prog, myObj);
    prog = sossolve(prog);
    %------------------------------ SOSp End ------------------------------
    V_ = sosgetsol(prog,V);
    V_old = V_;
end

% Save the improved value function
Vnew = V_;

%% Post-processing results
x0 = [0,0,0,0];  % Iniial Condition
tIntv = [0 3];
[t1,y1] = ode23s(@(t,x) LocalSuspSys(t,x,u), [0 3], x0);
[t,y] = ode23s(@(t,x) LocalSuspSys(t,x,0), tIntv, x0);

%% Plot Results
figure(1)
subplot(411);
plot(t1,y1(:,1),t,y(:,1), 'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_1','FontSize',12)
hl = legend('Improved performance', 'Uncontrolled performance');
set(hl, 'FontSize', 12');
subplot(412);
plot(t1,y1(:,2),t,y(:,2),'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_2','FontSize',12)
hl = legend('Improved performance', 'Uncontrolled performance');
set(hl, 'FontSize', 12');
subplot(413);
plot(t1,y1(:,3),t,y(:,3), 'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_3','FontSize',12)
hl = legend('Improved performance', 'Uncontrolled performance');
set(hl, 'FontSize', 12');
subplot(414);
plot(t1,y1(:,4),t,y(:,4),'r:','linewidth',2);
xlabel('time (sec)', 'FontSize',12)
ylabel('x_4','FontSize',12)
hl = legend('Improved performance', 'Uncontrolled performance');
set(hl, 'FontSize', 12');


figure(2)
xx1 = -.4:.04:.4;
xx2 = -5:0.5:5;
vn = zeros(length(xx1),length(xx2));
v1 = zeros(length(xx1),length(xx2));
un=vn;
ulqr = un;
kn=vn;
k1=v1;
x3=0;
x4=0;
for i=1:length(xx1)
    x1 = xx1(i);
    for j=1:length(xx2)
        x2 = xx2(j);
        vn(i,j)=eval(Vnew);
        v1(i,j)=eval(V0);
    end
end
surf(xx1,xx2,vn')
hold on
surf(xx1,xx2,v1')
hold off
xlabel('x_1', 'FontSize', 12)
ylabel('x_2', 'FontSize', 12)
view(gca,[-30.5 28]);
% Create textarrow
annotation(gcf,'textarrow',[0.210714285714286 0.174535137214669],...
    [0.895238095238095 0.631440045897884],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_0(x1,x2,0,0)'});

% Create textarrow
annotation(gcf,'textarrow',[0.139285714285714 0.186735060271868],...
    [0.183333333333333 0.386454388984516],'TextEdgeColor','none','FontSize',12,...
    'String',{'V_{10}(x1,x2,0,0)'});
% export_fig Ex3_cost -pdf -transparent

end


%% LocalSuspSys
% Dynamics of the nonlinear suspension system
function dx = LocalSuspSys(t,x,u)
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

[x1,x2,x3,x4] = deal(x(1),x(2),x(3),x(4));

% State matrices
A = [ 0 1 0 0;
    [-ks -bs ks bs]/mb ; ...
    0 0 0 1;
    [ks bs -ks-kt -bs]/mw];
B = [ 0; 10000/mb ; 0 ; -10000/mw];
B1 = [ 0; 0 ; 0 ; kt/mw];

if ~isdouble(u)
    u = eval(u);
end
    
if t <= 0.001
    r = 10;
else
    r = 0;
end

dx = A*x + B*u + B1*r;
end
