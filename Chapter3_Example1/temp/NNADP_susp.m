classdef NNADP_susp < handle	
	properties
		PsiPsi = []
		PsiU = []
		Phi = []
		CostQ = []
		PhiLength % = numel(Phi_fun(zeros(1,4)));
		PsiLength % = numel(Psi_fun(zeros(1,4)));
		r = 1     % Weight on u, r*u^2		
		N = 100   % Number of intervals 
		IterMax = 10  % Number of iterations
		X = [1,.3,1,-.2,zeros(1,PsiLength^2+PsiLength+1)];
	end
	
	methods
	end	
end



% ==================== Functions for system dynamics ==================== %
function dx = susp_sys(x,u)
%% SUSP_SYS decribes the dynamics of the suspension system
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

% Coefficients
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m
kn = ks/10;  % N/m

% System Dynamics
dx1 = x2;
dx2 = -(ks*(x1-x3)+kn*(x1-x3)^3)/mb-(bs*(x2-x4)-u)/mb;
dx3 = x4;
dx4 =  (ks*(x1-x3)+kn*(x1-x3)^3)/mw +(bs*(x2-x4)-kt*x3-u)/mw;

% Combine the output
dx = [dx1;
	  dx2;
	  dx3;
	  dx4];
end

function dX = adpSysWrapper(t,X)
%% ADPSYSWRAPPER augments the systems dynamics function by adding integrators

% System dynamics part
x = X(1:4);
u = sum(0.1*sin([1 3 7 11 13 15]*t));

dx = susp_sys(x,u);       % dx as the first 1-4 states of the wrapper

% Augmented part

psi = Psi_fun(x);

PsiPsi = kron(psi,psi);   % Export \psi\otimes\psi as the 4+[1,24^2] states
Psiu   = psi*u;           % Export \psi*u as the 4+[1,24^2]+ [1,24] states
Q = x'*x;                 % Export q(x) as the last state

dX = [dx; 
	  PsiPsi; 
	  Psiu;
	  Q];
end


% =========================== Utility functions ========================= %
function y = Phi_fun(x)
%% PHI_FUN is the function to approximate the value function V(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
%
y  = [x1^2
	  x2^2
	  x3^2
	  x4^2
	  x1*x2
	  x1*x3
	  x1*x4
	  x2*x3
	  x2*x4
	  x3*x4
	  x1^4
	  x2^4
	  x3^4
	  x4^4
	  x1^2*x2^2
	  x1^2*x3^2
	  x1^2*x4^2
	  x2^2*x3^2
	  x2^2*x4^2
	  x3^2*x4^2	  
	]';
% Notice that the output will be a ROW vector
end


function y = Psi_fun(x)
%% PSI_FUN are the basic functions to approximate u(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
%
y  = [x1                %  #1
      x2                %  #2
      x3                %  #3
      x4                %  #4
      x1*x1*x1          %  #5
      x1*x1*x2          %  #6
      x1*x1*x3          %  #7
      x1*x1*x4          %  #8
      x1*x2*x2          %  #9
      x1*x2*x3          %  #10
      x1*x2*x4          %  #11
      x1*x3*x3          %  #12
      x1*x3*x4          %  #13
      x1*x4*x4          %  #14
      x2*x2*x2          %  #15
      x2*x2*x3          %  #16
      x2*x2*x4          %  #17
      x2*x3*x3          %  #18
      x2*x3*x4          %  #19
      x2*x4*x4          %  #20
      x3*x3*x3          %  #21
      x3*x3*x4          %  #22
      x3*x4*x4          %  #23
      x4*x4*x4          %  #24
	];
end

