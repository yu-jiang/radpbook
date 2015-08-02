function dx = susp_sys_linear(x,u, A, B)

% Dynamics of the syspension system
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

% Coefficients
mdlCoeff

% System Dynamics
dx1 = x2;
dx2 = -(ks*(x1-x3)+kn*(x1-x3)^3)/mb - (bs*(x2-x4)-10000*u)/mb;
dx3 = x4;
dx4 =  (ks*(x1-x3)+kn*(x1-x3)^3)/mw + (bs*(x2-x4)-kt*x3-10000*u)/mw;


% Combine the output
dx = [dx1;
	  dx2;
	  dx3;
	  dx4]-(A*x+B*u);
end