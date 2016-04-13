function [t,Y] = LocalSDESolver(t0,tf,x0,dt, K, isForceField)
h=dt;                                                            
t=t0:dt:tf;                                                      
y=x0;
Y=[];
for clock = t
	%clock
	y = y + Motor2D(y,K,isForceField)*h;
	Y = [Y y];
end
Y=Y';
end