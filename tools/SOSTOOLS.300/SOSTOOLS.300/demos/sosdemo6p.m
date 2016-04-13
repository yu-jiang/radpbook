% SOSDEMO6 --- MAX CUT
% Section 3.6 of SOSTOOLS User's Manual
% 

clear; echo on;
pvar x1 x2 x3 x4 x5;
vartable = [x1; x2; x3; x4; x5];

% Number of cuts
f = 2.5 - 0.5*x1*x2 - 0.5*x2*x3 - 0.5*x3*x4 - 0.5*x4*x5 - 0.5*x5*x1;

% Boolean constraints
bc{1} = x1^2 - 1 ;
bc{2} = x2^2 - 1 ;
bc{3} = x3^2 - 1 ;
bc{4} = x4^2 - 1 ;
bc{5} = x5^2 - 1 ;

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% Then define SOSP variables

% -- p1(x) -- : sum of squares
% Monomial vector: 5 independent variables, degree <= 1
Z = monomials(vartable,[0 1]); 
[prog,p{1}] = sossosvar(prog,Z);

% -- p2(x) ... p6(x) : polynomials
% Monomial vector: 5 independent variables, degree <= 2
Z = monomials(vartable,0:2);
for i = 1:5
    [prog,p{1+i}] = sospolyvar(prog,Z);
end;

% =============================================
% Next, define SOSP constraints

% Constraint : p1(x)*(gamma - f(x)) +  p2(x)*bc1(x)
%               + ... + p6(x)*bc5(x) - (gamma-f(x))^2 >= 0
gamma = 4;

expr = p{1}*(gamma-f);
for i = 2:6
    expr = expr + p{i}*bc{i-1};
end;
expr = expr - (gamma-f)^2;

prog = sosineq(prog,expr);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Program is feasible, thus 4 is an upper bound for the cut.
echo off