% SOSDEMO1 --- Sum of Squares Test
% Section 3.1 of SOSTOOLS User's Manual
% 

clear; echo on;
pvar x1 x2;
vartable = [x1, x2];
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);   % No decision variables.

% =============================================
% Next, define the inequality

% p(x1,x2) >=  0
p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;
prog = sosineq(prog,p);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Program is feasible, thus p(x1,x2) is an SOS.
echo off;