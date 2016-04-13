% SOSDEMO9p --- Matrix SOS decomposition
% Section 3.9 of SOSTOOLS User's Manual
% 
clear; echo on;
pvar x1 x2 x3;
vartable = [x1, x2, x3];
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);   % No decision variables.

% =============================================
% Consider the following candidate sum of squares matrix P(x)

P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2); 
     x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2];  

% Test if P(x1,x2,x3) is an SOS matrix and return H so that P = H.'*H
[Q,Z,H] = findsos(P);

% Verify that P - H'*H = 0 and P- kron(I,Z)'*Q*kron(I,Z)= 0 to within 
% numerical tolerance.
tol = 1e-8;
cleanpoly(P - H'*H ,tol)
cleanpoly(P- kron(eye(2),Z)'*Q*kron(eye(2),Z), tol)

% Verify that Q >=0
min(eig(Q))

% =============================================
% Program is feasible, thus P(x1,x2,x3) is an SOS matrix.

echo off;
