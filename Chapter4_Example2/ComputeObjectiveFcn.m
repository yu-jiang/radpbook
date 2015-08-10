function c = ComputeObjectiveFcn(x1min, x1max, x2min, x2max)
% In the SOS-base Policy Iteration, we need to solve an SOSp in each
% iteraton. The objective of the SOSp is 
%
% min integration{V(x)}_Omega
%
% where Omega is a compact set, an area of interested of the system
% performance.

syms x1 x2
v_basis_fcn=[ x1*x1;
    x1*x2;
    x2*x2;
    x1*x1*x1;
    x1*x1*x2;
    x1*x2*x2;
    x2*x2*x2;
    x1*x1*x1*x1;
    x1*x1*x1*x2;
    x1*x1*x2*x2;
    x1*x2*x2*x2;
    x2*x2*x2*x2;];
c = double(int(int(v_basis_fcn,x1min,x1max),x2min,x2max));
end