% SOSDEMO10p --- Set containment
% Section 3.10 of SOSTOOLS User's Manual
% 
clear; echo on;
pvar x1 x2; 
vartable = [x1 x2];

eps = 1e-6;
% =============================================
% This is the problem data
p = x1^2+x2^2;
gamma = 1;
g0 = [2 0]*[x1;x2];
theta = 1;

% =============================================
% Initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% The multiplier
Zmon = monomials(vartable,0:4);
[prog,s] = sospolymatrixvar(prog,Zmon,[1 1]);

% =============================================
% Term to be added to g0
Zmon = monomials(vartable,2:3);
[prog,g1] = sospolymatrixvar(prog,Zmon,[1 1]);

% =============================================
% The expression to satisfy the set containment
Sc = [theta^2-s*(gamma-p) g0+g1; g0+g1 1];

prog = sosmatrixineq(prog,Sc-eps*eye(2));

option.solver = 'sdpt3';
option.params.vers = 2;
option.params.gaptol = 1e-7;
prog = sossolve(prog,option);

s = sosgetsol(prog,s);
g1 = sosgetsol(prog,g1);

% Display function g1 after removing small coefficients
cleanpoly(g1,1e-4)

% =============================================
% Program is feasible, { x |((g0+g1) + theta)(theta - (g0+g1)) >=0 } contains { x | p <= gamma }
echo off;

