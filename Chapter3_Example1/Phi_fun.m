function y = Phi_fun(x)
% Basis function for V(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

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