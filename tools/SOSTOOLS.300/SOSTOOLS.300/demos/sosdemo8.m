% SOSDEMO8 --- Bounds in Probability
% Section 3.8 of SOSTOOLS User's Manual
 
clear; echo on;
syms x a b c;

% The probability adds up to one.
m0 = 1 ;

% Mean
m1  = 1 ;

% Variance
sig = 1/2 ;

% E(x^2)
m2 = sig^2+m1^2; 

% Support of the random variable
R = [0,5];

% Event whose probability we want to bound
E = [4,5];

% =============================================
% Constructing and solving the SOS program
prog = sosprogram([x],[a,b,c]);

P = a + b*x + c*x^2 ;

% Nonnegative on the support
prog = sosineq(prog, P ,R);

% Greater than one on the event
prog = sosineq(prog,P-1,E);

% The bound
bnd =  a * m0 + b * m1 + c * m2 ;

% Objective: minimize the bound
prog = sossetobj(prog, bnd) ;

prog = sossolve(prog);

% =============================================
% Get solution
BND = sosgetsol(prog,bnd,16)
PP = sosgetsol(prog,P);

echo off;
