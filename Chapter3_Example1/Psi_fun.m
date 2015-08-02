function y = Psi_fun(x)
% Derivative of Phi

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

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