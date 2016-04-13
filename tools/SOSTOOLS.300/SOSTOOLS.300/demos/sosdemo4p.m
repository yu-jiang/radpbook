% SOSDEMO4 --- Matrix Copositivity
% Section 3.4 of SOSTOOLS User's Manual
% 

clear; echo on;
pvar x1 x2 x3 x4 x5;
vartable = [x1; x2; x3; x4; x5];

% The matrix under consideration
J = [1 -1  1  1 -1;
    -1  1 -1  1  1;
     1 -1  1 -1  1;
     1  1 -1  1 -1;
    -1  1  1 -1  1];

% =============================================
% First, initialize the sum of squares program

prog = sosprogram(vartable);     % No decision variables.

% =============================================
% Next, define SOSP constraints

% Constraint : r(x)*J(x) - p(x) = 0
J = [x1^2 x2^2 x3^2 x4^2 x5^2]*J*[x1^2; x2^2; x3^2; x4^2; x5^2];
r = x1^2 + x2^2 + x3^2 + x4^2 + x5^2;

prog = sosineq(prog,r*J);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Program is feasible. The matrix J is copositive.

echo off





