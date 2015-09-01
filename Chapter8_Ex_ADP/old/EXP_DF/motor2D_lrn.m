function dX = armDynamics(t,X)

para;

global Kadp

x = X(1:6);

u = -Kadp*x; 

w = randn(2,1)*sqrt(dt);

M = [c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];

v = M*w; % Noise from CNS

% External integrators for learning purpose
dX = [A*x+B*u+B*v./dt;  % dim = 6 
      x'*Q1*x+u'*R*u;   % dim = 1
      kron(x,R*v)./dt;  % dim = 12
      kron(v,v)./dt];   % dim = 4;
end