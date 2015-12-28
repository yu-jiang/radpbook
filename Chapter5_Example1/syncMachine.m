function dxx = syncMachine(t,xx, pmgr)

xx = xx(:);

% states of the physical machines
x = xx(1:3);
e = xx(4:6);

% Power deviation
Pd = pmgr.B11*pmgr.Ef1*pmgr.Ef2*(sin(x(1) + ...
	pmgr.angle10 - e(1) - pmgr.angle20) - ...
	sin(pmgr.angle10 - pmgr.angle20));

if t<3
    u = pmgr.K1*x  + 0.1*sin(100*t);
else
    u = -pmgr.KM*x;
end

dx = pmgr.A1*x + pmgr.B1*u - [0; Pd*pmgr.w0/2/pmgr.H1;0];

% System decomposition
A = pmgr.A1(1:2,1:2);
B = pmgr.A1(1:2,3);
F = pmgr.A1(3,3);
G = pmgr.B1(3);

xi = xx(20);
dxi = (pmgr.K0*(A-B*pmgr.K0)-F*pmgr.K0)*x(1:2) + (F + pmgr.K0*B)*xi + G*u;

% Dynamics of the following machine
de = (pmgr.A2 + pmgr.B2*pmgr.K2)*e + [0;Pd*pmgr.w0/2/pmgr.H2;0];

% Prepare the integrators for learning purpose
dIxx = kron(x(1:2),x(1:2));
dIxu = x(1:2)*(u);
dIzz = x(3)^2;
dIxz = x(1:2)*x(3);
dIzu = x(3)*(u);
dIdx = x(1:2)*(-Pd);
dIdz = x(3)*(-Pd);
dIxiu = xi*u;
dIxixi = xi*xi;
dIxix = xi*x(1:2);
xc = xx(25:27);
ec = xx(28:30);

% 
Pdc = pmgr.B11*pmgr.Ef1*pmgr.Ef2*(sin(xc(1) + pmgr.angle10 - ...
	ec(1) - pmgr.angle20) - sin(pmgr.angle10 - pmgr.angle20));
dxc = (pmgr.A1+pmgr.B1*pmgr.K1)*xc-[0;Pdc*pmgr.w0/2/pmgr.H1;0];
dec = (pmgr.A2+pmgr.B2*pmgr.K2)*ec+[0;Pdc*pmgr.w0/2/pmgr.H2;0];


dxx = [dx;
	   de;
	   dIxx;
	   dIxu;
	   dIxz;
	   dIzz;
	   dIzu;
	   dIdx;
	   dIdz;
	   dxi;
	   dIxiu;
	   dIxixi;
	   dIxix;
	   dxc;
	   dec];

%% Old implementation
% x = xx(1:3);
% z = xx(4:6);
% 
% if t >= pmgr.oscTstart
% 	u = sum(sin(pmgr.w*10*t))/5 - pmgr.K1*x;
% else
% 	u = -pmgr.K1*x;
% end
% 
% dIxx = kron(x,x);
% dIxu = kron(x,u);
% 
% if t >= pmgr.connectTstart
% 	dx = (pmgr.A1 - pmgr.B1*pmgr.K1)*x + ...
% 		[0;0;1]*(pmgr.B11*pmgr.Ef1*cos(x(1) - pmgr.angle10)*x(2) + pmgr.B12*pmgr.Ef1*pmgr.Ef2*cos(x(1)-z(1) - pmgr.angle10 + pmgr.angle20)*(x(2)-z(2)));
% 	dz = (pmgr.A2 - pmgr.B2*pmgr.K2)*z + ...
% 		[0;0;1]*(pmgr.B22*pmgr.Ef2*cos(z(1) - pmgr.angle20)*z(2) + pmgr.B12*pmgr.Ef1*pmgr.Ef2*cos(z(1)-x(1) + pmgr.angle10 - pmgr.angle20)*(z(2)-x(2)));
% else
% 	dx = pmgr.A1*x + pmgr.B1*u;
% 	dz = zeros(3,1);
% end
% 
% dxx = [dx;
% 	   dz;
% 	   dIxx;
% 	   dIxu];

end