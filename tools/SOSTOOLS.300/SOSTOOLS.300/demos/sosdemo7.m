% SOSDEMO7 --- Chebyshev polynomials
% Section 3.7 of SOSTOOLS User's Manual

clear; echo on;
 
ndeg = 8;   % Degree of Chebyshev polynomial

syms x gam;

% =============================================
% First, initialize the sum of squares program
prog = sosprogram([x],[gam]);


% Create the polynomial P
Z = monomials(x,[0:ndeg-1]);
[prog,P1] = sospolyvar(prog,Z);
P = P1 + gam * x^ndeg;           % The leading coeff of P is gam

% Imposing the inequalities
prog = sosineq(prog, 1 - P, [-1, 1]);
prog = sosineq(prog, 1 + P, [-1, 1]);

% And setting objective
prog = sossetobj(prog, -gam);

% Then solve the program
prog = sossolve(prog);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog, P)
GAM = sosgetsol(prog, gam)

echo off
