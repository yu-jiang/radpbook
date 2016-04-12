function dX = Motor2D(X,K,isForceField)

% Parameters
m   = 2;%1;  % mass 
b   = 10; % viscosity
c1  = 0.15/2;
c2  = 0.05/2;
dt  = 0.005;
tau = 0.05;
A = [ 0 0    1    0;
      0 0    0    1;
      0 0 -b/m    0;
      0 0    0 -b/m];
B = [0   0;
     0   0;
     1/m 0;
     0 1/m];
B1 = [1/tau 0;
     0  1/tau];

x = X(1:6);
z = X(7:8);

dx = x;
x = x(:);
u = -K*x;

w = randn(2,1)*sqrt(dt);
M = [c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];

VF=[-10.1 -11.2; -11.2 11.1];

v = M*w; % control dependent noise
f = z;

dx(1:4) = A * x(1:4)+B*(x(5:6)+f);
dx(5:6) = B1 * (-x(5:6)+u+v./dt);

if isForceField
	dz = -1/0.01*(z-VF*[x(3);x(4)]);
else
	dz=[0;0];
end

dX=[dx;dz];
end