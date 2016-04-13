% SOSDEMO3 --- Bound on Global Extremum
% Section 3.3 of SOSTOOLS User's Manual
% 

clear; echo on;
pvar x1 x2 gam;
vartable = [x1, x2];
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% Declare decision variable gam too
prog = sosdecvar(prog,gam);

% =============================================
% Next, define SOSP constraints

% Constraint : r(x)*(f(x) - gam) >= 0
% f(x) is the Goldstein-Price function
f1 = x1+x2+1;
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2;
f3 = 2*x1-3*x2;
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2;

f = (1+f1^2*f2)*(30+f3^2*f4);

prog = sosineq(prog,(f-gam));

% =============================================
% Set objective : maximize gam
prog = sossetobj(prog,-gam);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Finally, get solution
SOLgamma = sosgetsol(prog,gam)
echo off





