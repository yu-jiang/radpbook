% SOSDEMO2 --- Lyapunov Function Search 
% Section 3.2 of SOSTOOLS User's Manual
% 

clear; echo on;
pvar x1 x2 x3;
vars = [x1; x2; x3];

% Constructing the vector field dx/dt = f
f = [(-x1^3-x1*x3^2)*(x3^2+1);
    (-x2-x1^2*x2)*(x3^2+1);
    (-x3+3*x1^2*x3)*(x3^2+1)-3*x3];

% =============================================
% First, initialize the sum of squares program
%the following was used to test the sosprogram function by passing the decision variables as an argument GV
%syms dec1 dec2;
%decvartable = [dec1; dec2];
%prog = sosprogram(vars,decvartable);
%which is replacing the previous command

prog = sosprogram(vars);

% =============================================
% The Lyapunov function V(x): 
[prog,V] = sospolyvar(prog,[x1^2; x2^2; x3^2],'wscoeff');

% =============================================
% Next, define SOSP constraints

% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
prog = sosineq(prog,V-(x1^2+x2^2+x3^2));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2)+diff(V,x3)*f(3));
prog = sosineq(prog,expr);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,V)

echo off;



