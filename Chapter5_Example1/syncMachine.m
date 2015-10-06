function dxx = syncMachine(t,xx, pmgr)

x = xx(1:3);
z = xx(4:6);

if t >= pmgr.oscTstart
	u = sum(sin(pmgr.w*10*t))/5 - pmgr.K1*x;
else
	u = -pmgr.K1*x;
end

dIxx = kron(x,x);
dIxu = kron(x,u);

if t >= pmgr.connectTstart
	dx = (pmgr.A1 - pmgr.B1*pmgr.K1)*x + ...
		[0;0;1]*(pmgr.B11*pmgr.Ef1*cos(x(1) - pmgr.angle10)*x(2) + pmgr.B12*pmgr.Ef1*pmgr.Ef2*cos(x(1)-z(1) - pmgr.angle10 + pmgr.angle20)*(x(2)-z(2)));
	dz = (pmgr.A2 - pmgr.B2*pmgr.K2)*z + ...
		[0;0;1]*(pmgr.B22*pmgr.Ef2*cos(z(1) - pmgr.angle20)*z(2) + pmgr.B12*pmgr.Ef1*pmgr.Ef2*cos(z(1)-x(1) + pmgr.angle10 - pmgr.angle20)*(z(2)-x(2)));
else
	dx = pmgr.A1*x + pmgr.B1*u;
	dz = zeros(3,1);
end

dxx = [dx;
	   dz;
	   dIxx;
	   dIxu];

end